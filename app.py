"""
TBAnalytica — Flask Web UI
Provides a browser-based interface to the TB variant analysis pipeline.
Runs analyses in background threads and reports progress via polling.
"""

import io
import logging
import os
import sys
import threading
import uuid
import zipfile
from datetime import datetime
from pathlib import Path

from flask import (
    Flask, render_template, request, jsonify, send_file, abort, url_for,
)

# ── Project path setup ─────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from schema import (
    TBVariant, Mutation, ResistanceGene, ResistanceLevel,
    DrugSensitivity, ClinicalReport,
    DataSource, QUALITY_GATE_THRESHOLD,
)
from modules.variant_db import (
    load_known_variants, get_variant_by_id, add_variant,
    load_resistance_db,
)
from modules.api_calls import (
    fetch_ncbi_gene_sequence, search_ncbi_gene,
    fetch_ncbi_metadata, blast_sequence,
    fetch_uniprot_sequence, fetch_alphafold_prediction,
    fetch_all_gene_data, generate_quality_summary,
    build_quality_score, TB_REVIEWED_UNIPROT,
)
from modules.protein_compare import find_closest_protein_match
from modules.gene_analysis import full_gene_analysis, find_all_mutations
from modules.similarity_score import (
    full_similarity_analysis, calculate_risk_score,
    interpret_and_recommend, predict_resistance_level,
)
from modules.treatment import (
    get_treatment_for_known_variant, generate_treatment_recommendation,
    infer_treatment_from_similarity, format_treatment_summary,
    check_drug_interactions,
)
from modules.report_generator import generate_doctor_report, generate_patient_report
from modules.charts import (
    generate_all_charts, resistance_profile_chart, mutation_frequency_chart,
    comparison_heatmap, treatment_availability_chart, risk_score_gauge,
)

# ── Flask app ──────────────────────────────────────────────────────────────

app = Flask(__name__)
app.secret_key = os.urandom(24)

LOG_DIR = PROJECT_ROOT / "logs"
LOG_DIR.mkdir(exist_ok=True)
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(exist_ok=True)
CHARTS_DIR = OUTPUT_DIR / "charts"
CHARTS_DIR.mkdir(exist_ok=True)

log = logging.getLogger("TBAnalytica.web")

# ── Job store ──────────────────────────────────────────────────────────────

jobs: dict[str, dict] = {}

KNOWN_STEPS = [
    "Loading variant from database",
    "Checking data quality",
    "Calculating risk score",
    "Generating treatment recommendation",
    "Comparing to known variants",
    "Generating PDF reports",
    "Generating charts",
]

NEW_STEPS = [
    "Fetching sequence data from NCBI",
    "Building variant profile",
    "Fetching protein data from UniProt",
    "Running protein comparison",
    "Analyzing gene mutations",
    "Generating similarity scores",
    "Generating treatment recommendation",
    "Generating PDF reports",
    "Generating charts",
    "Saving variant to database",
]

COMPARE_STEPS = [
    "Loading first variant",
    "Fetching second variant data",
    "Running comparison analysis",
    "Calculating risk scores",
    "Checking data quality",
]

DEMO_STEPS = [
    "Loading known variants database",
    "Analysing known MDR variant",
    "Building new unknown variant",
    "Running risk and similarity analysis",
    "Checking data quality",
    "Generating reports and charts",
]

# ── Lineage mapping ───────────────────────────────────────────────────────

LINEAGE_ALIASES = {
    "beijing": "Lineage 2 (Beijing)",
    "cas": "Lineage 3",
    "cas/delhi": "Lineage 3",
    "central asian": "Lineage 3",
    "indo-oceanic": "Lineage 1 (Indo-Oceanic)",
    "indo oceanic": "Lineage 1 (Indo-Oceanic)",
    "euro-american": "Lineage 4 (Euro-American)",
    "euro american": "Lineage 4 (Euro-American)",
    "lineage 1": "Lineage 1",
    "lineage 2": "Lineage 2",
    "lineage 3": "Lineage 3",
    "lineage 4": "Lineage 4",
}


def _resolve_variant_input(user_input: str) -> TBVariant | None:
    variant = get_variant_by_id(user_input)
    if variant:
        return variant
    search_lineage = LINEAGE_ALIASES.get(user_input.lower(), user_input)
    known = load_known_variants()
    for v in known:
        if v.lineage and search_lineage.lower() in v.lineage.lower():
            return v
    for v in known:
        if user_input.lower() in v.variant_id.lower():
            return v
        if v.name and user_input.lower() in v.name.lower():
            return v
    return None


def identify_variant(input_str: str) -> dict:
    """Identify whether user input matches a known variant.

    Multi-step check: local DB → NCBI metadata cross-ref → BLAST similarity.
    Returns dict with keys: type, variant, source, ncbi_metadata, matched_on.
    """
    base: dict = {
        "type": "new", "variant": None, "source": None,
        "ncbi_metadata": {}, "matched_on": "",
    }

    # ── Step 1: Direct local DB match ────────────────────────────────────
    local_match = _resolve_variant_input(input_str)
    if local_match:
        log.info("identify_variant: '%s' → local DB match %s", input_str, local_match.variant_id)
        return {**base, "type": "known", "variant": local_match,
                "source": "local_db", "matched_on": f"variant ID / name '{input_str}'"}

    # ── Step 2: Fetch NCBI metadata + cross-reference ────────────────────
    ncbi_meta: dict = {}
    try:
        ncbi_meta = fetch_ncbi_metadata(input_str)
    except Exception as exc:
        log.debug("identify_variant: NCBI metadata fetch failed: %s", exc)

    if ncbi_meta.get("organism") or ncbi_meta.get("strain") or ncbi_meta.get("title"):
        log.info(
            "identify_variant NCBI metadata: organism=%r  strain=%r  gene=%r  title=%r",
            ncbi_meta.get("organism", ""), ncbi_meta.get("strain", ""),
            ncbi_meta.get("gene", ""), ncbi_meta.get("title", ""),
        )

        keywords: list[str] = []
        strain = ncbi_meta.get("strain", "")
        organism = ncbi_meta.get("organism", "")
        title = ncbi_meta.get("title", "")

        if strain:
            keywords.append(strain)
        if organism:
            keywords.append(organism)
        if title:
            keywords.append(title)
        for hint, lineage_key in [
            ("H37Rv", "H37Rv"), ("Beijing", "Beijing"), ("CAS", "CAS"),
            ("Central Asian", "CAS"), ("Indo-Oceanic", "Indo-Oceanic"),
            ("Euro-American", "Euro-American"), ("Haarlem", "Euro-American"),
            ("LAM", "Euro-American"),
        ]:
            if hint.lower() in (strain + " " + organism + " " + title).lower():
                keywords.append(lineage_key)

        known = load_known_variants()
        for kw in keywords:
            if not kw:
                continue
            resolved = _resolve_variant_input(kw)
            if resolved:
                log.info(
                    "identify_variant: '%s' → NCBI keyword '%s' matched %s",
                    input_str, kw, resolved.variant_id,
                )
                return {**base, "type": "known", "variant": resolved,
                        "source": "ncbi_metadata", "ncbi_metadata": ncbi_meta,
                        "matched_on": f"NCBI metadata keyword '{kw}'"}

            kw_lower = kw.lower()
            for v in known:
                if v.name and kw_lower in v.name.lower():
                    log.info(
                        "identify_variant: '%s' → NCBI keyword '%s' matched name of %s",
                        input_str, kw, v.variant_id,
                    )
                    return {**base, "type": "known", "variant": v,
                            "source": "ncbi_metadata", "ncbi_metadata": ncbi_meta,
                            "matched_on": f"NCBI strain/name '{kw}'"}

    # ── Step 3: BLAST sequence similarity ────────────────────────────────
    sequence = ncbi_meta.get("sequence", "")
    if sequence and len(sequence) >= 30:
        log.info("identify_variant: running BLAST for '%s' (%d aa/bp)", input_str, len(sequence))
        try:
            blast_hits = blast_sequence(sequence, max_results=5)
            known = load_known_variants()
            for hit_data, hit_qs in blast_hits:
                hit_title = hit_data.get("title", "")
                pct_id = hit_data.get("percent_identity", 0)
                if pct_id < 95:
                    continue
                hit_lower = hit_title.lower()
                for v in known:
                    if v.name and v.name.lower() in hit_lower:
                        log.info(
                            "identify_variant: BLAST hit '%s' (%.1f%%) matched %s",
                            hit_title, pct_id, v.variant_id,
                        )
                        return {
                            **base, "type": "known", "variant": v,
                            "source": "blast_match", "ncbi_metadata": ncbi_meta,
                            "matched_on": f"BLAST hit '{hit_title}' ({pct_id:.1f}% identity)",
                        }
                    for rg in v.resistance_genes:
                        if rg.gene_name.lower() in hit_lower and pct_id >= 99:
                            log.info(
                                "identify_variant: BLAST gene '%s' (%.1f%%) → %s",
                                rg.gene_name, pct_id, v.variant_id,
                            )
                            return {
                                **base, "type": "known", "variant": v,
                                "source": "blast_match", "ncbi_metadata": ncbi_meta,
                                "matched_on": (
                                    f"BLAST gene '{rg.gene_name}' in "
                                    f"'{hit_title}' ({pct_id:.1f}% identity)"
                                ),
                            }
        except Exception as exc:
            log.debug("identify_variant: BLAST failed (non-fatal): %s", exc)

    # ── Step 4: Genuinely new ────────────────────────────────────────────
    log.info("identify_variant: '%s' → new variant (no DB/NCBI/BLAST match)", input_str)
    return {**base, "ncbi_metadata": ncbi_meta}


# ── Shared pipeline helpers (mirrors main.py) ─────────────────────────────

def _generate_patient_id() -> str:
    return f"NP-{datetime.now().strftime('%Y%m%d-%H%M%S')}"


def _fetch_data_quality(gene_names):
    all_quality = {}
    all_warnings = []
    for gene in gene_names:
        try:
            combined, quality_scores, warnings = fetch_all_gene_data(gene)
            all_quality.update(quality_scores)
            all_warnings.extend(warnings)
        except Exception as e:
            log.warning("Data fetch failed for %s: %s", gene, e)
            fallback_qs = build_quality_score(DataSource.LOCAL_DB, review_status="fetch_failed")
            all_quality[f"{gene} (Local fallback)"] = fallback_qs
            all_warnings.append(f"All external sources failed for {gene} — using local DB only")
    for label, qs in all_quality.items():
        if qs.raw_score < QUALITY_GATE_THRESHOLD:
            all_warnings.append(
                f"Insufficient data quality for {label}: "
                f"{qs.raw_score:.0f}/100 ({qs.confidence}) — "
                f"recommend manual laboratory verification"
            )
    summary = generate_quality_summary(all_quality)
    return all_quality, all_warnings, summary


def _gate_treatment_on_quality(quality_scores, treatment_rec):
    gate_warnings = []
    gated = False
    low_sources = {
        label: qs for label, qs in quality_scores.items()
        if qs.raw_score < QUALITY_GATE_THRESHOLD
    }
    if low_sources:
        gated = True
        for label, qs in low_sources.items():
            gate_warnings.append(
                f"QUALITY GATE: {label} scored {qs.raw_score:.0f}/100 — "
                f"treatment recommendation requires manual verification"
            )
        treatment_rec.setdefault("warnings", []).append(
            "DATA QUALITY WARNING: One or more data sources scored below 60/100. "
            "Treatment recommendation must be verified by laboratory susceptibility testing "
            "before clinical use."
        )
        treatment_rec["quality_verified"] = False
        confidence = treatment_rec.get("confidence_score", 1.0)
        treatment_rec["confidence_score"] = min(confidence, 0.5)
        treatment_rec["confidence"] = "LOW"
    return treatment_rec, gated, gate_warnings


def _build_report(variant, risk, treatment_rec, comparison, patient_id="",
                  quality_scores=None, quality_warnings=None,
                  quality_text="", treatment_gated=False):
    treatment_text = format_treatment_summary(treatment_rec)
    doctor_notes = (
        f"Resistance class: {treatment_rec.get('resistance_class', 'Unknown')}. "
        f"Risk score: {risk.score:.1f}% ({risk.level.value}). "
        f"See treatment recommendation for regimen details."
    )
    if quality_text:
        doctor_notes += "\n\n" + quality_text
    return ClinicalReport(
        patient_id=patient_id or _generate_patient_id(),
        variant=variant,
        comparison_result=comparison,
        risk_score=risk,
        treatment=treatment_text,
        doctor_report=doctor_notes,
        patient_report=(
            f"Your TB sample has been analysed. "
            f"Risk level: {risk.level.value}. "
            f"Please follow the prescribed treatment plan carefully."
        ),
        data_quality=quality_scores or {},
        quality_warnings=quality_warnings or [],
        treatment_gated=treatment_gated,
    )


def _build_variant_from_sequence(sequence, gene_name, lineage="", source="",
                                 ncbi_metadata=None):
    db = load_resistance_db()
    gene_data = db["genes"].get(gene_name, {})
    ref_seq = gene_data.get("reference_sequence", "")
    mutations = []
    drug_resistance = {}
    if ref_seq and sequence:
        mutations = find_all_mutations(sequence, ref_seq, gene_name)
        for m in mutations:
            if m.is_resistance_conferring and m.drug_affected:
                for drug in m.drug_affected.split(","):
                    drug_resistance[drug.strip()] = DrugSensitivity.RESISTANT
    resistance_genes = []
    if mutations:
        resistance_genes.append(ResistanceGene(
            gene_name=gene_name, mutations=mutations,
            drug_target=gene_data.get("drug_target", ""),
        ))
    variant_id = f"TB_NEW_{datetime.now().strftime('%Y%m%d%H%M%S')}"
    is_protein = all(c.upper() in "ACDEFGHIKLMNPQRSTVWY*" for c in sequence if c.isalpha())
    protein_seqs = {gene_name: sequence} if is_protein else {}
    nuc_seqs = {} if is_protein else {gene_name: sequence}

    # Derive proper name from NCBI metadata
    meta = ncbi_metadata or {}
    strain = meta.get("strain", "")
    accession = meta.get("accession", "")
    ncbi_gene = meta.get("gene", "")
    if strain:
        variant_name = f"{strain}_{ncbi_gene or gene_name}"
    elif accession:
        variant_name = f"Unknown_{accession}"
    else:
        variant_name = f"Unknown_{gene_name}_{variant_id}"

    if not lineage and meta.get("organism"):
        lineage = meta["organism"]

    return TBVariant(
        variant_id=variant_id,
        name=variant_name,
        lineage=lineage, drug_resistance=drug_resistance,
        resistance_genes=resistance_genes,
        protein_sequences=protein_seqs, nucleotide_sequences=nuc_seqs,
        source=source or "Web analysis",
    )


# ── Analysis threads ──────────────────────────────────────────────────────

def _advance(job_id, step_label):
    job = jobs[job_id]
    job["progress"].append(step_label)


def _run_known_analysis(job_id, variant_id, patient_id):
    job = jobs[job_id]
    try:
        _advance(job_id, KNOWN_STEPS[0])
        variant = get_variant_by_id(variant_id)
        if not variant:
            resolved = _resolve_variant_input(variant_id)
            if not resolved:
                raise ValueError(f"Variant '{variant_id}' not found")
            variant = resolved

        _advance(job_id, KNOWN_STEPS[1])
        gene_names = [rg.gene_name for rg in variant.resistance_genes] or ["rpoB"]
        quality_scores, quality_warnings, quality_text = _fetch_data_quality(gene_names)

        _advance(job_id, KNOWN_STEPS[2])
        risk = calculate_risk_score(variant)

        _advance(job_id, KNOWN_STEPS[3])
        treatment_rec = get_treatment_for_known_variant(variant.variant_id)
        if "error" in treatment_rec:
            resistance_class = predict_resistance_level(risk.score, variant.all_mutations())
            treatment_rec = generate_treatment_recommendation(variant, resistance_class)
        treatment_rec, treatment_gated, gate_warnings = _gate_treatment_on_quality(
            quality_scores, treatment_rec,
        )
        quality_warnings.extend(gate_warnings)

        _advance(job_id, KNOWN_STEPS[4])
        known_variants = load_known_variants()
        comparisons = full_similarity_analysis(variant, known_variants)
        closest = comparisons[0] if comparisons else None

        _advance(job_id, KNOWN_STEPS[5])
        pid = patient_id or _generate_patient_id()
        report = _build_report(
            variant, risk, treatment_rec, closest, pid,
            quality_scores=quality_scores,
            quality_warnings=quality_warnings,
            quality_text=quality_text,
            treatment_gated=treatment_gated,
        )
        doctor_pdf = generate_doctor_report(report)
        patient_pdf = generate_patient_report(report)

        _advance(job_id, KNOWN_STEPS[6])
        chart_paths = {}
        try:
            chart_paths = generate_all_charts(report, comparisons, str(CHARTS_DIR))
            chart_paths["resistance_profile"] = resistance_profile_chart(variant)
            chart_paths["mutation_frequency"] = mutation_frequency_chart(variant)
            chart_paths["risk_gauge_png"] = risk_score_gauge(risk.score, variant.variant_id)
            if comparisons:
                chart_paths["comparison_heatmap"] = comparison_heatmap(comparisons)
            chart_paths["treatment_availability"] = treatment_availability_chart(treatment_rec)
        except Exception as e:
            log.warning("Chart generation failed: %s", e)

        job["status"] = "complete"
        job["result"] = _serialize_result(
            report, comparisons, treatment_rec,
            doctor_pdf, patient_pdf, chart_paths,
        )

    except Exception as e:
        log.exception("Known analysis failed for job %s", job_id)
        job["status"] = "error"
        job["error"] = str(e)


def _run_new_analysis(job_id, sequence, accession, gene_name, lineage, patient_id):
    job = jobs[job_id]
    try:
        # ── Variant identification gate ──────────────────────────────────
        # Before running the full new-variant pipeline, check whether the
        # input actually matches a known variant via local DB, NCBI metadata,
        # or BLAST similarity.
        lookup_input = accession or ""
        ncbi_meta: dict = {}
        if lookup_input:
            identification = identify_variant(lookup_input)
            ncbi_meta = identification.get("ncbi_metadata", {})

            if ncbi_meta.get("organism") or ncbi_meta.get("strain"):
                log.info(
                    "NCBI metadata for '%s': organism=%r strain=%r gene=%r title=%r seq_len=%d",
                    lookup_input, ncbi_meta.get("organism", ""),
                    ncbi_meta.get("strain", ""), ncbi_meta.get("gene", ""),
                    ncbi_meta.get("title", ""), len(ncbi_meta.get("sequence", "")),
                )

            if identification["type"] == "known":
                matched = identification["variant"]
                source_map = {
                    "local_db": "local database",
                    "ncbi_metadata": "NCBI metadata",
                    "blast_match": "BLAST sequence match",
                }
                source_label = source_map.get(
                    identification["source"], identification["source"],
                )
                log.info(
                    "New-variant request redirected → known %s (source: %s, match: %s)",
                    matched.variant_id, identification["source"],
                    identification.get("matched_on", ""),
                )
                job["steps"] = KNOWN_STEPS
                job["progress"] = []
                notification = (
                    f"This input was recognised as {matched.variant_id} "
                    f"({matched.name}) via {source_label}."
                )
                if identification.get("matched_on"):
                    notification += f" Matched on: {identification['matched_on']}."
                notification += " Running known-variant analysis instead."
                job["notification"] = notification
                _run_known_analysis(job_id, matched.variant_id, patient_id)
                return

        gene = gene_name or "rpoB"
        # If NCBI metadata detected a gene, use that as a hint
        if ncbi_meta.get("gene") and not gene_name:
            gene = ncbi_meta["gene"]
            log.info("Gene inferred from NCBI metadata: %s", gene)

        _advance(job_id, NEW_STEPS[0])
        # Re-use sequence from metadata if already fetched
        if ncbi_meta.get("sequence") and not sequence:
            sequence = ncbi_meta["sequence"]
            log.info("Using sequence from NCBI metadata: %d chars", len(sequence))
        elif accession and not sequence:
            sequence = fetch_ncbi_gene_sequence(accession)
            if not sequence:
                gene_id = search_ncbi_gene(gene, "Mycobacterium tuberculosis")
                if gene_id:
                    sequence = fetch_ncbi_gene_sequence(gene_id)
            if not sequence:
                raise ValueError(f"Could not fetch sequence for accession '{accession}'")

        _advance(job_id, NEW_STEPS[1])
        variant = _build_variant_from_sequence(
            sequence=sequence, gene_name=gene, lineage=lineage,
            source=f"NCBI:{accession}" if accession else "Direct sequence input",
            ncbi_metadata=ncbi_meta,
        )

        _advance(job_id, NEW_STEPS[2])
        gene_names = [rg.gene_name for rg in variant.resistance_genes] or [gene]
        quality_scores, quality_warnings, quality_text = _fetch_data_quality(gene_names)
        try:
            uniprot_id = TB_REVIEWED_UNIPROT.get(gene)
            if uniprot_id:
                protein_seq = fetch_uniprot_sequence(uniprot_id)
                if protein_seq:
                    variant.protein_sequences[gene] = protein_seq
                fetch_alphafold_prediction(uniprot_id)
        except Exception:
            pass

        _advance(job_id, NEW_STEPS[3])
        known_variants = load_known_variants()
        closest_protein, protein_score = find_closest_protein_match(variant, known_variants)

        _advance(job_id, NEW_STEPS[4])
        reference = closest_protein or (known_variants[0] if known_variants else None)
        if reference:
            try:
                full_gene_analysis(variant, reference)
            except Exception:
                pass

        _advance(job_id, NEW_STEPS[5])
        comparisons = full_similarity_analysis(variant, known_variants)
        closest = comparisons[0] if comparisons else None
        risk = calculate_risk_score(variant)
        resistance_class = predict_resistance_level(risk.score, variant.all_mutations())

        _advance(job_id, NEW_STEPS[6])
        if closest and closest.weighted_final_score >= 70:
            treatment_rec = infer_treatment_from_similarity(closest, variant)
        else:
            treatment_rec = generate_treatment_recommendation(variant, resistance_class)
        treatment_rec, treatment_gated, gate_warnings = _gate_treatment_on_quality(
            quality_scores, treatment_rec,
        )
        quality_warnings.extend(gate_warnings)
        interactions = check_drug_interactions(treatment_rec.get("recommended_regimen", []))
        if interactions:
            treatment_rec.setdefault("warnings", []).extend(interactions)

        _advance(job_id, NEW_STEPS[7])
        pid = patient_id or _generate_patient_id()
        report = _build_report(
            variant, risk, treatment_rec, closest, pid,
            quality_scores=quality_scores,
            quality_warnings=quality_warnings,
            quality_text=quality_text,
            treatment_gated=treatment_gated,
        )
        doctor_pdf = generate_doctor_report(report)
        patient_pdf = generate_patient_report(report)

        _advance(job_id, NEW_STEPS[8])
        chart_paths = {}
        try:
            chart_paths = generate_all_charts(report, comparisons, str(CHARTS_DIR))
            chart_paths["resistance_profile"] = resistance_profile_chart(variant)
            chart_paths["mutation_frequency"] = mutation_frequency_chart(variant)
            chart_paths["risk_gauge_png"] = risk_score_gauge(risk.score, variant.variant_id)
            if comparisons:
                chart_paths["comparison_heatmap"] = comparison_heatmap(comparisons)
            chart_paths["treatment_availability"] = treatment_availability_chart(treatment_rec)
        except Exception as e:
            log.warning("Chart generation failed: %s", e)

        _advance(job_id, NEW_STEPS[9])
        try:
            add_variant(variant)
        except Exception:
            pass

        job["status"] = "complete"
        job["result"] = _serialize_result(
            report, comparisons, treatment_rec,
            doctor_pdf, patient_pdf, chart_paths,
        )

    except Exception as e:
        log.exception("New analysis failed for job %s", job_id)
        job["status"] = "error"
        job["error"] = str(e)


def _run_compare(job_id, variant_id, second_input, gene_name):
    job = jobs[job_id]
    try:
        _advance(job_id, COMPARE_STEPS[0])
        known = get_variant_by_id(variant_id)
        if not known:
            resolved = _resolve_variant_input(variant_id)
            if not resolved:
                raise ValueError(f"Variant '{variant_id}' not found")
            known = resolved

        _advance(job_id, COMPARE_STEPS[1])
        second_variant = _resolve_variant_input(second_input)
        new_accession = None
        seq = None
        if second_variant:
            seqs = second_variant.protein_sequences or second_variant.nucleotide_sequences
            seq = list(seqs.values())[0] if seqs else "PLACEHOLDER"
        else:
            new_accession = second_input
            seq = fetch_ncbi_gene_sequence(new_accession)
            if not seq:
                raise ValueError(f"Could not fetch sequence for '{second_input}'")

        gene = gene_name or "rpoB"
        new_variant = _build_variant_from_sequence(
            seq, gene, source=f"Comparison:{new_accession or 'known'}",
        )

        _advance(job_id, COMPARE_STEPS[2])
        comparisons = full_similarity_analysis(new_variant, [known])
        comp = comparisons[0] if comparisons else None
        recommendation = interpret_and_recommend(comp) if comp else {}

        _advance(job_id, COMPARE_STEPS[3])
        known_risk = calculate_risk_score(known)
        new_risk = calculate_risk_score(new_variant)

        _advance(job_id, COMPARE_STEPS[4])
        gene_names = [rg.gene_name for rg in new_variant.resistance_genes] or [gene]
        quality_scores, quality_warnings, quality_text = _fetch_data_quality(gene_names)

        job["status"] = "complete"
        job["result"] = {
            "mode": "compare",
            "known_variant": {
                "variant_id": known.variant_id,
                "name": known.name,
                "lineage": known.lineage,
                "risk_score": known_risk.score,
                "risk_level": known_risk.level.value,
            },
            "new_variant": {
                "variant_id": new_variant.variant_id,
                "name": new_variant.name,
            },
            "comparison": {
                "similarity": comp.weighted_final_score if comp else 0,
                "protein_similarity": comp.protein_similarity_score if comp else 0,
                "gene_similarity": comp.gene_similarity_score if comp else 0,
                "confidence": comp.confidence_level.value if comp else "LOW",
            } if comp else None,
            "known_risk": {"score": known_risk.score, "level": known_risk.level.value},
            "new_risk": {"score": new_risk.score, "level": new_risk.level.value},
            "risk_delta": new_risk.score - known_risk.score,
            "recommendation": recommendation,
            "quality_scores": {
                label: {"score": qs.raw_score, "confidence": qs.confidence,
                        "source": qs.source.value}
                for label, qs in quality_scores.items()
            },
            "quality_warnings": quality_warnings,
        }

    except Exception as e:
        log.exception("Compare failed for job %s", job_id)
        job["status"] = "error"
        job["error"] = str(e)


def _run_demo(job_id):
    job = jobs[job_id]
    try:
        _advance(job_id, DEMO_STEPS[0])
        variants = load_known_variants()

        _advance(job_id, DEMO_STEPS[1])
        result = get_variant_by_id("TB_VAR_002")
        risk = calculate_risk_score(result) if result else None
        treatment_rec = get_treatment_for_known_variant("TB_VAR_002") if result else {}

        _advance(job_id, DEMO_STEPS[2])
        db = load_resistance_db()
        by_gene = {}
        drug_resistance = {}
        raw_mutations = [
            {"gene": "rpoB", "position": 531, "ref_amino_acid": "S",
             "alt_amino_acid": "L", "mutation_type": "missense"},
            {"gene": "katG", "position": 315, "ref_amino_acid": "S",
             "alt_amino_acid": "T", "mutation_type": "missense"},
            {"gene": "rpoB", "position": 514, "ref_nucleotide": "TTC",
             "alt_nucleotide": "TTT", "ref_amino_acid": "F",
             "alt_amino_acid": "F", "mutation_type": "silent"},
        ]
        for m in raw_mutations:
            gene = m["gene"]
            ref_aa = m.get("ref_amino_acid", "")
            alt_aa = m.get("alt_amino_acid", "")
            is_syn = (ref_aa == alt_aa) and ref_aa != ""
            is_res = m.get("mutation_type") not in ("silent", "promoter") and not is_syn
            gene_data = db["genes"].get(gene, {})
            mut_code = f"{ref_aa}{m['position']}{alt_aa}"
            mut_db_info = gene_data.get("mutations", {}).get(mut_code, {})
            drug_affected = ", ".join(mut_db_info.get("drugs", []))
            one_step = False
            one_step_drug = None
            if is_syn:
                for _, p in gene_data.get("silent_precursors", {}).items():
                    if p.get("position") == m["position"]:
                        one_step = True
                        one_step_drug = gene_data.get("drug_target", "")
            if is_res and mut_db_info.get("drugs"):
                for d in mut_db_info["drugs"]:
                    drug_resistance[d] = DrugSensitivity.RESISTANT
            mutation = Mutation(
                position=m["position"],
                reference_codon=m.get("ref_nucleotide", m.get("ref_codon", "")),
                mutant_codon=m.get("alt_nucleotide", m.get("alt_codon", "")),
                reference_amino_acid=ref_aa, mutant_amino_acid=alt_aa,
                is_synonymous=is_syn, is_resistance_conferring=is_res,
                drug_affected=drug_affected,
                resistance_level=(ResistanceLevel.HIGH
                                  if mut_db_info.get("resistance") == "high"
                                  else ResistanceLevel.MODERATE),
                one_step_away_risk=one_step, one_step_away_drug=one_step_drug,
            )
            by_gene.setdefault(gene, []).append(mutation)

        resistance_genes = []
        for gn, gm in by_gene.items():
            gd = db["genes"].get(gn, {})
            resistance_genes.append(ResistanceGene(
                gene_name=gn, mutations=gm, drug_target=gd.get("drug_target", ""),
            ))
        new_variant = TBVariant(
            variant_id=f"TB_NEW_{datetime.now().strftime('%Y%m%d%H%M%S')}",
            name="Nepal_Unknown_2024", lineage="Lineage 2 (Beijing)",
            drug_resistance=drug_resistance,
            resistance_genes=resistance_genes, source="Demo analysis",
        )

        _advance(job_id, DEMO_STEPS[3])
        nr = calculate_risk_score(new_variant)
        new_class = predict_resistance_level(nr.score, new_variant.all_mutations())
        new_comparisons = full_similarity_analysis(new_variant, load_known_variants())
        new_treatment = generate_treatment_recommendation(new_variant, new_class)

        _advance(job_id, DEMO_STEPS[4])
        demo_gene_names = list(by_gene.keys()) or ["rpoB"]
        demo_quality, demo_qwarnings, demo_qtxt = _fetch_data_quality(demo_gene_names)
        new_treatment, demo_gated, demo_gate_warnings = _gate_treatment_on_quality(
            demo_quality, new_treatment,
        )
        demo_qwarnings.extend(demo_gate_warnings)

        _advance(job_id, DEMO_STEPS[5])
        closest_comp = new_comparisons[0] if new_comparisons else None
        demo_report = _build_report(
            new_variant, nr, new_treatment, closest_comp, "NP-2024-001",
            quality_scores=demo_quality, quality_warnings=demo_qwarnings,
            quality_text=demo_qtxt, treatment_gated=demo_gated,
        )
        doctor_pdf = generate_doctor_report(demo_report)
        patient_pdf = generate_patient_report(demo_report)
        chart_paths = {}
        try:
            chart_paths["resistance_profile"] = resistance_profile_chart(new_variant)
            chart_paths["mutation_frequency"] = mutation_frequency_chart(new_variant)
            chart_paths["risk_gauge"] = risk_score_gauge(nr.score, new_variant.variant_id)
            if new_comparisons:
                chart_paths["comparison_heatmap"] = comparison_heatmap(new_comparisons)
            chart_paths["treatment_availability"] = treatment_availability_chart(new_treatment)
        except Exception:
            pass

        job["status"] = "complete"
        job["result"] = _serialize_result(
            demo_report, new_comparisons, new_treatment,
            doctor_pdf, patient_pdf, chart_paths,
        )

    except Exception as e:
        log.exception("Demo failed for job %s", job_id)
        job["status"] = "error"
        job["error"] = str(e)


# ── Result serializer ─────────────────────────────────────────────────────

def _serialize_result(report, comparisons, treatment_rec,
                      doctor_pdf, patient_pdf, chart_paths):
    variant = report.variant
    risk = report.risk_score

    if variant.has_xdr_profile():
        mdr_status = "XDR-TB CONFIRMED"
    elif variant.has_mdr_profile():
        mdr_status = "MDR-TB CONFIRMED"
    elif variant.resistant_drugs():
        mdr_status = f"Mono/Poly-resistant"
    else:
        mdr_status = "Drug-susceptible"

    mutations_list = []
    for rg in variant.resistance_genes:
        for m in rg.mutations:
            if m.is_resistance_conferring:
                mut_type = "resistance"
            elif m.is_synonymous and m.one_step_away_risk:
                mut_type = "precursor"
            elif m.is_synonymous:
                mut_type = "silent"
            else:
                mut_type = "missense"
            mutations_list.append({
                "gene": rg.gene_name,
                "position": m.position,
                "change": m.short_code(),
                "type": mut_type,
                "drug_affected": m.drug_affected or "-",
                "risk": (m.resistance_level.value.upper()
                         if m.is_resistance_conferring
                         else ("PRECURSOR" if m.one_step_away_risk else "low")),
            })

    comparisons_list = []
    for c in (comparisons or [])[:5]:
        comparisons_list.append({
            "variant_id": c.matched_variant_id,
            "weighted_score": c.weighted_final_score,
            "protein_score": c.protein_similarity_score,
            "gene_score": c.gene_similarity_score,
            "confidence": c.confidence_level.value,
        })

    quality_list = []
    avg_quality = 0
    for label, qs in report.data_quality.items():
        quality_list.append({
            "label": label,
            "source": qs.source.value,
            "score": qs.raw_score,
            "confidence": qs.confidence,
            "use_for_analysis": qs.use_for_analysis,
        })
        avg_quality += qs.raw_score
    if quality_list:
        avg_quality /= len(quality_list)

    chart_urls = {}
    for name, path in chart_paths.items():
        if path and os.path.exists(path):
            chart_urls[name] = os.path.relpath(path, str(PROJECT_ROOT)).replace("\\", "/")

    return {
        "mode": "analysis",
        "patient_id": report.patient_id,
        "variant": {
            "variant_id": variant.variant_id,
            "name": variant.name,
            "lineage": variant.lineage or "Unknown",
        },
        "risk": {
            "score": risk.score,
            "level": risk.level.value,
            "color": risk.color.value,
            "factors": risk.factors,
        },
        "mdr_status": mdr_status,
        "resistant_drugs": variant.resistant_drugs(),
        "treatment": {
            "class": treatment_rec.get("resistance_class", "N/A"),
            "regimen": treatment_rec.get("regimen_string", "N/A"),
            "duration_weeks": treatment_rec.get("duration_weeks", "?"),
            "ntp_compliant": treatment_rec.get("nepal_ntp_compliant", "N/A"),
            "warnings": treatment_rec.get("warnings", []),
        },
        "doctor_report": report.doctor_report,
        "patient_report": report.patient_report,
        "treatment_text": report.treatment,
        "mutations": mutations_list,
        "comparisons": comparisons_list,
        "quality": {
            "scores": quality_list,
            "avg_score": round(avg_quality, 1),
            "warnings": report.quality_warnings,
            "treatment_gated": report.treatment_gated,
        },
        "charts": chart_urls,
        "files": {
            "doctor_pdf": doctor_pdf,
            "patient_pdf": patient_pdf,
        },
    }


# ── Flask routes ──────────────────────────────────────────────────────────

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/analyze")
def analyze_form():
    mode = request.args.get("mode", "known")
    known_variants = []
    try:
        known_variants = [
            {"id": v.variant_id, "name": v.name, "lineage": v.lineage,
             "mdr": v.has_mdr_profile(), "xdr": v.has_xdr_profile()}
            for v in load_known_variants()
        ]
    except Exception:
        pass
    return render_template("index.html", mode=mode, known_variants=known_variants)


@app.route("/run-analysis", methods=["POST"])
def run_analysis():
    data = request.get_json() or request.form.to_dict()
    mode = data.get("mode", "known")
    job_id = str(uuid.uuid4())

    if mode == "known":
        step_list = KNOWN_STEPS
    elif mode == "new":
        step_list = NEW_STEPS
    elif mode == "compare":
        step_list = COMPARE_STEPS
    elif mode == "demo":
        step_list = DEMO_STEPS
    else:
        return jsonify({"error": f"Unknown mode: {mode}"}), 400

    jobs[job_id] = {
        "status": "running",
        "progress": [],
        "steps": step_list,
        "result": None,
        "error": None,
        "notification": None,
    }

    if mode == "known":
        variant_input = data.get("variant_id", "").strip()
        patient_id = data.get("patient_id", "").strip()
        if not variant_input:
            return jsonify({"error": "variant_id is required"}), 400
        t = threading.Thread(
            target=_run_known_analysis,
            args=(job_id, variant_input, patient_id),
            daemon=True,
        )
    elif mode == "new":
        sequence = data.get("sequence", "").strip()
        accession = data.get("accession", "").strip()
        gene_name = data.get("gene_name", "rpoB").strip()
        lineage = data.get("lineage", "").strip()
        patient_id = data.get("patient_id", "").strip()
        if not sequence and not accession:
            return jsonify({"error": "sequence or accession required"}), 400
        t = threading.Thread(
            target=_run_new_analysis,
            args=(job_id, sequence or None, accession or None,
                  gene_name, lineage, patient_id),
            daemon=True,
        )
    elif mode == "compare":
        variant_id = data.get("variant_id", "").strip()
        second_input = data.get("second_input", "").strip()
        gene_name = data.get("gene_name", "rpoB").strip()
        if not variant_id or not second_input:
            return jsonify({"error": "both variant IDs are required"}), 400
        t = threading.Thread(
            target=_run_compare,
            args=(job_id, variant_id, second_input, gene_name),
            daemon=True,
        )
    elif mode == "demo":
        t = threading.Thread(target=_run_demo, args=(job_id,), daemon=True)

    t.start()
    return jsonify({"job_id": job_id, "steps": step_list})


@app.route("/status/<job_id>")
def job_status(job_id):
    job = jobs.get(job_id)
    if not job:
        return jsonify({"error": "Job not found"}), 404
    resp = {
        "status": job["status"],
        "progress": job["progress"],
        "steps": job["steps"],
        "error": job["error"],
    }
    if job.get("notification"):
        resp["notification"] = job["notification"]
    return jsonify(resp)


@app.route("/results/<job_id>")
def results_page(job_id):
    job = jobs.get(job_id)
    if not job:
        abort(404)
    if job["status"] != "complete":
        return render_template("index.html", error="Analysis not yet complete.")
    return render_template(
        "results.html", job_id=job_id, data=job["result"],
        notification=job.get("notification"),
    )


@app.route("/api/results/<job_id>")
def results_json(job_id):
    job = jobs.get(job_id)
    if not job:
        return jsonify({"error": "Job not found"}), 404
    if job["status"] != "complete":
        return jsonify({"error": "Not complete"}), 202
    return jsonify(job["result"])


@app.route("/download/<job_id>/<file_type>")
def download_file(job_id, file_type):
    job = jobs.get(job_id)
    if not job or job["status"] != "complete":
        abort(404)

    result = job["result"]

    if file_type == "doctor_pdf":
        path = result["files"]["doctor_pdf"]
        if os.path.exists(path):
            return send_file(path, as_attachment=True,
                             download_name=f"clinical_report_{result['patient_id']}.pdf")
    elif file_type == "patient_pdf":
        path = result["files"]["patient_pdf"]
        if os.path.exists(path):
            return send_file(path, as_attachment=True,
                             download_name=f"patient_summary_{result['patient_id']}.pdf")
    elif file_type == "charts_zip":
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            for name, rel_path in result.get("charts", {}).items():
                full_path = str(PROJECT_ROOT / rel_path)
                if os.path.exists(full_path):
                    zf.write(full_path, f"charts/{os.path.basename(full_path)}")
        buf.seek(0)
        return send_file(buf, as_attachment=True,
                         download_name=f"charts_{result['patient_id']}.zip",
                         mimetype="application/zip")

    abort(404)


@app.route("/chart/<path:chart_path>")
def serve_chart(chart_path):
    full_path = PROJECT_ROOT / chart_path
    if full_path.exists() and full_path.suffix in (".png", ".html", ".svg"):
        return send_file(str(full_path))
    abort(404)


# ── Entry point ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 44)
    print("  TBAnalytica Web UI")
    print("  Open browser at: http://localhost:5000")
    print("=" * 44)
    app.run(debug=True, port=5000)
