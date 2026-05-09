"""
TBAnalytica — TB Variant Analysis and Treatment Recommendation System
Clinical Decision Support for Nepal

Interactive entry point that orchestrates all modules.
Provides both an interactive console interface and a function API for web integration.
"""

import logging
import os
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from schema import (
    TBVariant, Mutation, ResistanceGene, ResistanceLevel,
    DrugSensitivity, ComparisonResult, RiskScore, ClinicalReport,
    DataQualityScore, DataSource, QUALITY_GATE_THRESHOLD,
)

# ── Module imports ──────────────────────────────────────────────────────────

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
from modules.gene_analysis import (
    full_gene_analysis, find_all_mutations,
    assess_silent_mutation_risk,
)
from modules.similarity_score import (
    full_similarity_analysis, calculate_risk_score,
    interpret_and_recommend, predict_resistance_level,
)
from modules.treatment import (
    get_treatment_for_known_variant, generate_treatment_recommendation,
    infer_treatment_from_similarity, format_treatment_summary,
    get_treatment_centers, check_drug_interactions,
)
from modules.report_generator import (
    generate_doctor_report, generate_patient_report,
)
from modules.charts import (
    generate_all_charts, resistance_profile_chart,
    mutation_frequency_chart, comparison_heatmap,
    treatment_availability_chart, risk_score_gauge,
)

# ── Logging ─────────────────────────────────────────────────────────────────

LOG_DIR = Path(__file__).resolve().parent / "logs"
LOG_DIR.mkdir(exist_ok=True)

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

CHARTS_DIR = OUTPUT_DIR / "charts"
CHARTS_DIR.mkdir(exist_ok=True)


def _setup_logging(verbose: bool = False) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    log_file = LOG_DIR / f"tbanalytica_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    logger = logging.getLogger("TBAnalytica")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    fh = logging.FileHandler(str(log_file), encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S",
    ))
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(level)
    ch.setFormatter(logging.Formatter("%(levelname)-8s | %(message)s"))
    logger.addHandler(ch)

    return logger


log = _setup_logging()


# ── Interactive input helpers ──────────────────────────────────────────────

def _ask_choice(prompt: str, valid: range) -> int:
    while True:
        try:
            raw = input(prompt).strip()
            val = int(raw)
            if val in valid:
                return val
            print(f"  Please enter a number between {valid.start} and {valid.stop - 1}.")
        except ValueError:
            print(f"  Invalid input. Please enter a number between {valid.start} and {valid.stop - 1}.")
        except EOFError:
            raise SystemExit(0)


def _ask_text(prompt: str, allow_empty: bool = False) -> str:
    while True:
        try:
            raw = input(prompt).strip()
            if raw or allow_empty:
                return raw
            print("  This field is required. Please enter a value.")
        except EOFError:
            raise SystemExit(0)


def _ask_yes_no(prompt: str) -> bool:
    while True:
        try:
            raw = input(prompt).strip().lower()
            if raw in ("yes", "y"):
                return True
            if raw in ("no", "n"):
                return False
            print("  Please enter yes or no.")
        except EOFError:
            raise SystemExit(0)


# ── Progress display ──────────────────────────────────────────────────────

class StepProgress:
    """Simple step-by-step progress display with checkmarks."""

    def __init__(self, steps: list[str]):
        self._steps = steps
        self._total = len(steps)
        self._current = 0

    def step(self, label: str | None = None):
        self._current += 1
        desc = label or (self._steps[self._current - 1] if self._current <= self._total else "")
        print(f"  [{self._current}/{self._total}] {desc}...", end="", flush=True)

    def done(self):
        print("  done")

    def fail(self, msg: str = ""):
        print(f"  FAILED{': ' + msg if msg else ''}")


# ── Core helpers ───────────────────────────────────────────────────────────

def _generate_patient_id() -> str:
    return f"NP-{datetime.now().strftime('%Y%m%d-%H%M%S')}"


def _build_variant_from_sequence(
    sequence: str,
    gene_name: str,
    lineage: str = "",
    source: str = "",
    ncbi_metadata: dict | None = None,
) -> TBVariant:
    db = load_resistance_db()
    gene_data = db["genes"].get(gene_name, {})
    ref_seq = gene_data.get("reference_sequence", "")

    mutations: list[Mutation] = []
    drug_resistance: dict[str, DrugSensitivity] = {}

    if ref_seq and sequence:
        mutations = find_all_mutations(sequence, ref_seq, gene_name)
        for m in mutations:
            if m.is_resistance_conferring and m.drug_affected:
                for drug in m.drug_affected.split(","):
                    drug_resistance[drug.strip()] = DrugSensitivity.RESISTANT

    resistance_genes = []
    if mutations:
        resistance_genes.append(ResistanceGene(
            gene_name=gene_name,
            mutations=mutations,
            drug_target=gene_data.get("drug_target", ""),
        ))

    variant_id = f"TB_NEW_{datetime.now().strftime('%Y%m%d%H%M%S')}"
    is_protein = all(c.upper() in "ACDEFGHIKLMNPQRSTVWY*" for c in sequence if c.isalpha())
    if is_protein:
        protein_seqs = {gene_name: sequence}
        nuc_seqs: dict[str, str] = {}
    else:
        protein_seqs = {}
        nuc_seqs = {gene_name: sequence}

    # ── Derive a proper name from NCBI metadata ─────────────────────────
    # Priority: strain name → accession-based ID → fallback
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

    # Use lineage from NCBI organism field if not explicitly provided
    if not lineage and meta.get("organism"):
        lineage = meta["organism"]

    return TBVariant(
        variant_id=variant_id,
        name=variant_name,
        lineage=lineage,
        drug_resistance=drug_resistance,
        resistance_genes=resistance_genes,
        protein_sequences=protein_seqs,
        nucleotide_sequences=nuc_seqs,
        source=source or "Interactive analysis",
    )


def _fetch_data_quality(gene_names: list[str]) -> tuple[
    dict[str, DataQualityScore], list[str], str,
]:
    all_quality: dict[str, DataQualityScore] = {}
    all_warnings: list[str] = []

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
            log.warning("Quality gate FAILED for %s: %.0f/100", label, qs.raw_score)

    summary = generate_quality_summary(all_quality)
    return all_quality, all_warnings, summary


def _gate_treatment_on_quality(
    quality_scores: dict[str, DataQualityScore],
    treatment_rec: dict,
) -> tuple[dict, bool, list[str]]:
    gate_warnings: list[str] = []
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

        log.warning(
            "Treatment gated: %d source(s) below threshold — %s",
            len(low_sources), ", ".join(low_sources.keys()),
        )

    return treatment_rec, gated, gate_warnings


def _build_report(
    variant: TBVariant,
    risk: RiskScore,
    treatment_rec: dict,
    comparison: ComparisonResult | None,
    patient_id: str = "",
    quality_scores: dict[str, DataQualityScore] | None = None,
    quality_warnings: list[str] | None = None,
    quality_text: str = "",
    treatment_gated: bool = False,
) -> ClinicalReport:
    treatment_text = format_treatment_summary(treatment_rec)
    doctor_notes = (
        f"Resistance class: {treatment_rec.get('resistance_class', 'Unknown')}. "
        f"Risk score: {risk.score:.1f}% ({risk.level.value}). "
        f"See treatment recommendation for regimen details."
    )
    if quality_text:
        doctor_notes += "\n\n" + quality_text

    qw = quality_warnings or []
    qs = quality_scores or {}

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
        data_quality=qs,
        quality_warnings=qw,
        treatment_gated=treatment_gated,
    )


def _resolve_variant_input(user_input: str) -> TBVariant | None:
    """Resolve user input to a TBVariant — accepts variant ID or lineage name."""
    variant = get_variant_by_id(user_input)
    if variant:
        return variant

    lineage_aliases = {
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

    search_lineage = lineage_aliases.get(user_input.lower(), user_input)

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

    Performs a multi-step check before committing to the new-variant pipeline:
      1. Direct local DB lookup (variant ID, lineage alias, strain name)
      2. NCBI metadata fetch + cross-reference (strain, gene, organism)
      3. BLAST sequence similarity against known variants (optional, slow)
      4. Fall-through: genuinely new — return NCBI metadata for proper naming

    Returns dict with keys:
        type          — "known" or "new"
        variant       — TBVariant if known, else None
        source        — "local_db" | "ncbi_metadata" | "blast_match" | None
        ncbi_metadata — dict from fetch_ncbi_metadata (always populated for
                        accession inputs, empty dict otherwise)
        matched_on    — human-readable string describing what matched (or "")
    """
    base = {
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

        # Build keyword list from NCBI metadata to match against local DB
        keywords: list[str] = []
        strain = ncbi_meta.get("strain", "")
        organism = ncbi_meta.get("organism", "")
        title = ncbi_meta.get("title", "")
        gene = ncbi_meta.get("gene", "")

        if strain:
            keywords.append(strain)
        if organism:
            keywords.append(organism)
        if title:
            keywords.append(title)
        # Extract well-known lineage hints
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
            # Try resolving through the lineage-alias-aware function
            resolved = _resolve_variant_input(kw)
            if resolved:
                log.info(
                    "identify_variant: '%s' → NCBI metadata keyword '%s' matched %s",
                    input_str, kw, resolved.variant_id,
                )
                return {**base, "type": "known", "variant": resolved,
                        "source": "ncbi_metadata", "ncbi_metadata": ncbi_meta,
                        "matched_on": f"NCBI metadata keyword '{kw}'"}

            # Also try substring match against all variant names / lineages
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

                # Match gene name against variant's resistance genes
                if gene and any(gene.lower() == rg.gene_name.lower()
                                for rg in v.resistance_genes):
                    # Gene alone isn't enough — only match if strain also matches
                    pass

    # ── Step 3: BLAST sequence similarity (only for accession-based input) ──
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
                    continue  # only high-identity hits

                # Check if BLAST hit title matches any known variant name/strain
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
                    # Check strain keywords
                    for rg in v.resistance_genes:
                        if rg.gene_name.lower() in hit_lower:
                            if pct_id >= 99:
                                log.info(
                                    "identify_variant: BLAST hit gene '%s' (%.1f%%) → %s",
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


# ── 1. analyze_known_variant ───────────────────────────────────────────────

def analyze_known_variant(
    variant_id: str,
    patient_id: str = "",
    generate_charts_flag: bool = True,
) -> ClinicalReport:
    """Full analysis pipeline for a known variant."""
    log.info("Starting analysis for known variant: %s", variant_id)

    steps = [
        "Loading variant from database",
        "Checking data quality",
        "Calculating risk score",
        "Generating treatment recommendation",
        "Comparing to known variants",
        "Generating PDF reports",
        "Generating charts",
    ]
    progress = StepProgress(steps)

    # Step 1 — Load variant
    progress.step("Loading variant from database")
    variant = get_variant_by_id(variant_id)
    if not variant:
        progress.fail("not found")
        log.error("Variant %s not found in database", variant_id)
        raise ValueError(f"Variant '{variant_id}' not found in the known variants database")
    log.info("Loaded variant: %s (%s)", variant.variant_id, variant.name)
    progress.done()

    # Step 2 — Data quality check
    progress.step("Checking data quality")
    gene_names = [rg.gene_name for rg in variant.resistance_genes] or ["rpoB"]
    quality_scores, quality_warnings, quality_text = _fetch_data_quality(gene_names)
    log.info("Data quality: %d sources checked, %d warnings", len(quality_scores), len(quality_warnings))
    progress.done()

    # Step 3 — Risk score
    progress.step("Calculating risk score")
    risk = calculate_risk_score(variant)
    log.info("Risk score: %.1f%% [%s]", risk.score, risk.level.value)
    progress.done()

    # Step 4 — Treatment (gated on quality)
    progress.step("Generating treatment recommendation")
    treatment_rec = get_treatment_for_known_variant(variant_id)
    if "error" in treatment_rec:
        log.warning("Treatment lookup returned error: %s", treatment_rec["error"])
        resistance_class = predict_resistance_level(risk.score, variant.all_mutations())
        treatment_rec = generate_treatment_recommendation(variant, resistance_class)

    treatment_rec, treatment_gated, gate_warnings = _gate_treatment_on_quality(
        quality_scores, treatment_rec,
    )
    quality_warnings.extend(gate_warnings)
    log.info("Treatment regimen: %s (gated=%s)", treatment_rec.get("regimen_string", "N/A"), treatment_gated)
    progress.done()

    # Step 5 — Similarity comparison
    progress.step("Comparing to known variants")
    known_variants = load_known_variants()
    comparisons = full_similarity_analysis(variant, known_variants)
    closest = comparisons[0] if comparisons else None
    if closest:
        log.info("Closest match: %s (%.1f%%)", closest.matched_variant_id, closest.weighted_final_score)
    progress.done()

    # Step 6 — Reports
    progress.step("Generating PDF reports")
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
    log.info("Doctor report: %s", doctor_pdf)
    log.info("Patient report: %s", patient_pdf)
    progress.done()

    # Step 7 — Charts
    progress.step("Generating charts")
    chart_paths: dict[str, str] = {}
    if generate_charts_flag:
        try:
            chart_paths = generate_all_charts(report, comparisons, str(CHARTS_DIR))
            chart_paths["resistance_profile"] = resistance_profile_chart(variant)
            chart_paths["mutation_frequency"] = mutation_frequency_chart(variant)
            chart_paths["risk_gauge_png"] = risk_score_gauge(risk.score, variant.variant_id)
            if comparisons:
                chart_paths["comparison_heatmap"] = comparison_heatmap(comparisons)
            chart_paths["treatment_availability"] = treatment_availability_chart(treatment_rec)
            log.info("Generated %d charts", len(chart_paths))
        except Exception as e:
            log.warning("Chart generation failed: %s", e)
    progress.done()

    print()
    print("  Done! Reports saved to output/")
    print()

    _print_results(report, doctor_pdf, patient_pdf, chart_paths, comparisons, treatment_rec)
    _offer_open_files(doctor_pdf, patient_pdf, chart_paths)

    return report


# ── 2. analyze_new_variant ─────────────────────────────────────────────────

def analyze_new_variant(
    sequence: str | None = None,
    accession: str | None = None,
    gene_name: str | None = None,
    lineage: str = "",
    patient_id: str = "",
    generate_charts_flag: bool = True,
) -> ClinicalReport:
    """Full analysis pipeline for a new/unknown variant."""
    if not sequence and not accession:
        raise ValueError("Provide either 'sequence' or 'accession'")

    # ── Variant identification gate ──────────────────────────────────────
    # Before running the full new-variant pipeline, check whether the input
    # actually matches a known variant in our database or via NCBI metadata.
    lookup_input = accession or ""
    ncbi_meta: dict = {}
    if lookup_input:
        identification = identify_variant(lookup_input)
        ncbi_meta = identification.get("ncbi_metadata", {})

        # Debug: show what NCBI returned for this accession
        if ncbi_meta.get("organism") or ncbi_meta.get("strain"):
            print()
            print(f"  [DEBUG] NCBI metadata for '{lookup_input}':")
            print(f"          organism : {ncbi_meta.get('organism', '(none)')}")
            print(f"          strain   : {ncbi_meta.get('strain', '(none)')}")
            print(f"          gene     : {ncbi_meta.get('gene', '(none)')}")
            print(f"          title    : {ncbi_meta.get('title', '(none)')}")
            print(f"          db       : {ncbi_meta.get('db', '(none)')}")
            print(f"          seq len  : {len(ncbi_meta.get('sequence', ''))}")

        if identification["type"] == "known":
            matched = identification["variant"]
            source_map = {
                "local_db": "local database",
                "ncbi_metadata": "NCBI metadata",
                "blast_match": "BLAST sequence match",
            }
            source_label = source_map.get(identification["source"], identification["source"])
            print()
            print(f"  ** This input was recognised as {matched.variant_id} "
                  f"({matched.name}) via {source_label}.")
            if identification.get("matched_on"):
                print(f"     Matched on: {identification['matched_on']}")
            print(f"  ** Redirecting to known-variant analysis pipeline.")
            print()
            log.info(
                "New-variant request redirected → known variant %s (source: %s, match: %s)",
                matched.variant_id, identification["source"],
                identification.get("matched_on", ""),
            )
            return analyze_known_variant(
                variant_id=matched.variant_id,
                patient_id=patient_id,
                generate_charts_flag=generate_charts_flag,
            )

    gene = gene_name or "rpoB"
    # If NCBI metadata detected a gene, use that as a hint
    if ncbi_meta.get("gene") and not gene_name:
        gene = ncbi_meta["gene"]
        log.info("Gene inferred from NCBI metadata: %s", gene)

    log.info("Starting analysis for new variant (gene=%s, accession=%s)", gene, accession or "N/A")

    steps = [
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
    progress = StepProgress(steps)

    # Step 1 — Fetch sequence
    progress.step("Fetching sequence data from NCBI")
    # Re-use the sequence from metadata if already fetched
    if ncbi_meta.get("sequence") and not sequence:
        sequence = ncbi_meta["sequence"]
        log.info("Using sequence from NCBI metadata: %d characters", len(sequence))
    elif accession and not sequence:
        log.info("Fetching sequence from NCBI: %s", accession)
        sequence = fetch_ncbi_gene_sequence(accession)
        if not sequence:
            gene_id = search_ncbi_gene(gene, "Mycobacterium tuberculosis")
            if gene_id:
                sequence = fetch_ncbi_gene_sequence(gene_id)
        if not sequence:
            progress.fail("could not fetch sequence")
            log.error("Could not fetch sequence for accession %s", accession)
            raise ValueError(f"Could not fetch sequence for accession '{accession}'")
        log.info("Fetched sequence: %d characters", len(sequence))
    progress.done()

    # Step 2 — Build variant
    progress.step("Building variant profile")
    variant = _build_variant_from_sequence(
        sequence=sequence,
        gene_name=gene,
        lineage=lineage,
        source=f"NCBI:{accession}" if accession else "Direct sequence input",
        ncbi_metadata=ncbi_meta,
    )
    log.info("Built variant: %s with %d mutations", variant.variant_id, len(variant.all_mutations()))
    progress.done()

    # Step 3 — External data + quality check
    progress.step("Fetching protein data from UniProt")
    gene_names = [rg.gene_name for rg in variant.resistance_genes] or [gene]
    quality_scores, quality_warnings, quality_text = _fetch_data_quality(gene_names)
    log.info("Data quality: %d sources checked, %d warnings", len(quality_scores), len(quality_warnings))

    try:
        uniprot_id = TB_REVIEWED_UNIPROT.get(gene)
        if uniprot_id:
            protein_seq = fetch_uniprot_sequence(uniprot_id)
            if protein_seq:
                variant.protein_sequences[gene] = protein_seq
                log.info("UniProt sequence: %d aa", len(protein_seq))

            af_data = fetch_alphafold_prediction(uniprot_id)
            if af_data:
                log.info("AlphaFold data retrieved")
    except Exception as e:
        log.warning("External data fetch failed (non-fatal): %s", e)
    progress.done()

    # Step 4 — Protein comparison
    progress.step("Running protein comparison")
    known_variants = load_known_variants()
    closest_protein, protein_score = find_closest_protein_match(variant, known_variants)
    if closest_protein:
        log.info("Closest protein match: %s (%.1f%%)", closest_protein.variant_id, protein_score * 100)
    progress.done()

    # Step 5 — Gene analysis
    progress.step("Analyzing gene mutations")
    reference = closest_protein or (known_variants[0] if known_variants else None)
    gene_report: dict = {}
    if reference:
        try:
            gene_report = full_gene_analysis(variant, reference)
            log.info("Gene analysis complete: %d gene changes", len(gene_report.get("all_mutations", [])))
        except Exception as e:
            log.warning("Gene analysis failed (non-fatal): %s", e)
    progress.done()

    # Step 6 — Similarity analysis
    progress.step("Generating similarity scores")
    comparisons = full_similarity_analysis(variant, known_variants)
    closest = comparisons[0] if comparisons else None
    recommendation: dict = {}
    if closest:
        recommendation = interpret_and_recommend(closest)
        log.info(
            "Top match: %s (%.1f%%, confidence=%s)",
            closest.matched_variant_id, closest.weighted_final_score,
            recommendation.get("confidence", "N/A"),
        )
    risk = calculate_risk_score(variant)
    resistance_class = predict_resistance_level(risk.score, variant.all_mutations())
    log.info("Risk: %.1f%% [%s], class: %s", risk.score, risk.level.value, resistance_class)
    progress.done()

    # Step 7 — Treatment (gated on quality)
    progress.step("Generating treatment recommendation")
    if closest and closest.weighted_final_score >= 70:
        treatment_rec = infer_treatment_from_similarity(closest, variant)
        log.info("Treatment inferred from similarity to %s", closest.matched_variant_id)
    else:
        treatment_rec = generate_treatment_recommendation(variant, resistance_class)
        log.info("Treatment generated from resistance class: %s", resistance_class)

    treatment_rec, treatment_gated, gate_warnings = _gate_treatment_on_quality(
        quality_scores, treatment_rec,
    )
    quality_warnings.extend(gate_warnings)

    interactions = check_drug_interactions(treatment_rec.get("recommended_regimen", []))
    if interactions:
        treatment_rec.setdefault("warnings", []).extend(interactions)
        log.info("Drug interactions detected: %d", len(interactions))
    log.info("Treatment regimen: %s (gated=%s)", treatment_rec.get("regimen_string", "N/A"), treatment_gated)
    progress.done()

    # Step 8 — Reports
    progress.step("Generating PDF reports")
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
    log.info("Doctor report: %s", doctor_pdf)
    log.info("Patient report: %s", patient_pdf)
    progress.done()

    # Step 9 — Charts
    progress.step("Generating charts")
    chart_paths: dict[str, str] = {}
    if generate_charts_flag:
        try:
            chart_paths = generate_all_charts(report, comparisons, str(CHARTS_DIR))
            chart_paths["resistance_profile"] = resistance_profile_chart(variant)
            chart_paths["mutation_frequency"] = mutation_frequency_chart(variant)
            chart_paths["risk_gauge_png"] = risk_score_gauge(risk.score, variant.variant_id)
            if comparisons:
                chart_paths["comparison_heatmap"] = comparison_heatmap(comparisons)
            chart_paths["treatment_availability"] = treatment_availability_chart(treatment_rec)
            log.info("Generated %d charts", len(chart_paths))
        except Exception as e:
            log.warning("Chart generation failed: %s", e)
    progress.done()

    # Step 10 — Save to database
    progress.step("Saving variant to database")
    try:
        add_variant(variant)
        log.info("Variant %s saved to database", variant.variant_id)
    except Exception as e:
        log.warning("Could not save variant to database: %s", e)
    progress.done()

    print()
    print("  Done! Reports saved to output/")
    print()

    _print_results(report, doctor_pdf, patient_pdf, chart_paths, comparisons, treatment_rec)
    _offer_open_files(doctor_pdf, patient_pdf, chart_paths)

    return report


# ── 3. compare_variants ────────────────────────────────────────────────────

def compare_variants(
    variant_id: str,
    new_accession: str | None = None,
    new_sequence: str | None = None,
    gene_name: str = "rpoB",
) -> dict:
    """Compare a known variant to a new sequence/accession."""
    log.info("Comparing variant %s to new data", variant_id)

    known = get_variant_by_id(variant_id)
    if not known:
        raise ValueError(f"Known variant '{variant_id}' not found")

    seq = new_sequence
    if new_accession and not seq:
        seq = fetch_ncbi_gene_sequence(new_accession)
        if not seq:
            raise ValueError(f"Could not fetch sequence for accession '{new_accession}'")

    if not seq:
        raise ValueError("Provide either 'new_sequence' or 'new_accession'")

    new_variant = _build_variant_from_sequence(seq, gene_name, source=f"Comparison:{new_accession or 'direct'}")

    comparisons = full_similarity_analysis(new_variant, [known])
    comp = comparisons[0] if comparisons else None
    recommendation = interpret_and_recommend(comp) if comp else {}

    gene_report: dict = {}
    try:
        gene_report = full_gene_analysis(new_variant, known)
    except Exception as e:
        log.warning("Gene analysis in comparison failed: %s", e)

    known_risk = calculate_risk_score(known)
    new_risk = calculate_risk_score(new_variant)

    gene_names = [rg.gene_name for rg in new_variant.resistance_genes] or [gene_name]
    quality_scores, quality_warnings, quality_text = _fetch_data_quality(gene_names)

    return {
        "known_variant": known,
        "new_variant": new_variant,
        "comparison": comp,
        "recommendation": recommendation,
        "gene_report": gene_report,
        "known_risk": known_risk,
        "new_risk": new_risk,
        "risk_delta": new_risk.score - known_risk.score,
        "quality_scores": quality_scores,
        "quality_warnings": quality_warnings,
        "quality_text": quality_text,
    }


# ── 4. run_demo ────────────────────────────────────────────────────────────

def run_demo():
    """Demo mode showing the complete pipeline end-to-end."""

    print()
    print("  Running full demo pipeline...")
    print("  " + "-" * 40)
    print()

    # ── 1. Load database
    print("  [1] Loading known variants database...")
    variants = load_known_variants()
    print(f"      Loaded {len(variants)} known variants")
    for v in variants[:3]:
        print(f"      - {v.variant_id}: {v.name} (lineage={v.lineage}, MDR={v.has_mdr_profile()})")
    print()

    # ── 2. Analyse known MDR variant
    print("  [2] Analysing known MDR variant (TB_VAR_002)...")
    result = get_variant_by_id("TB_VAR_002")
    if result:
        risk = calculate_risk_score(result)
        resistance_class = predict_resistance_level(risk.score, result.all_mutations())
        treatment_rec = get_treatment_for_known_variant("TB_VAR_002")
        print(f"      Name: {result.name}")
        print(f"      Risk Score: {risk.score:.1f}% [{risk.level.value}] ({risk.color.value})")
        print(f"      Resistance Class: {resistance_class}")
        print(f"      Resistant drugs: {', '.join(result.resistant_drugs())}")
        print()
        print(format_treatment_summary(treatment_rec))
    print()

    # ── 3. JSON round-trip
    print("  [3] JSON serialization demo...")
    if result:
        print(f"      TBVariant.to_json() length: {len(result.to_json())} chars")
        roundtrip = TBVariant.from_json(result.to_json())
        print(f"      Round-trip OK: variant_id={roundtrip.variant_id}, mutations={len(roundtrip.all_mutations())}")
        print(f"      RiskScore.to_json(): {risk.to_json()[:120]}...")
    print()

    # ── 4. Analyse NEW unknown variant
    print("  [4] Analysing a NEW unknown variant...")
    db = load_resistance_db()
    by_gene: dict[str, list[Mutation]] = {}
    drug_resistance: dict[str, DrugSensitivity] = {}

    raw_mutations = [
        {"gene": "rpoB", "position": 531, "ref_amino_acid": "S", "alt_amino_acid": "L", "mutation_type": "missense"},
        {"gene": "katG", "position": 315, "ref_amino_acid": "S", "alt_amino_acid": "T", "mutation_type": "missense"},
        {"gene": "rpoB", "position": 514, "ref_nucleotide": "TTC", "alt_nucleotide": "TTT",
         "ref_amino_acid": "F", "alt_amino_acid": "F", "mutation_type": "silent"},
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
            precursors = gene_data.get("silent_precursors", {})
            for _, p in precursors.items():
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
            reference_amino_acid=ref_aa,
            mutant_amino_acid=alt_aa,
            is_synonymous=is_syn,
            is_resistance_conferring=is_res,
            drug_affected=drug_affected,
            resistance_level=(
                ResistanceLevel.HIGH if mut_db_info.get("resistance") == "high"
                else ResistanceLevel.MODERATE
            ),
            one_step_away_risk=one_step,
            one_step_away_drug=one_step_drug,
        )
        by_gene.setdefault(gene, []).append(mutation)

    resistance_genes = []
    for gene_name_key, gene_mutations in by_gene.items():
        gd = db["genes"].get(gene_name_key, {})
        resistance_genes.append(ResistanceGene(
            gene_name=gene_name_key,
            mutations=gene_mutations,
            drug_target=gd.get("drug_target", ""),
        ))

    new_variant = TBVariant(
        variant_id=f"TB_NEW_{datetime.now().strftime('%Y%m%d%H%M%S')}",
        name="Nepal_Unknown_2024",
        lineage="Lineage 2 (Beijing)",
        drug_resistance=drug_resistance,
        resistance_genes=resistance_genes,
        source="Demo analysis",
    )

    nr = calculate_risk_score(new_variant)
    new_class = predict_resistance_level(nr.score, new_variant.all_mutations())
    new_comparisons = full_similarity_analysis(new_variant, load_known_variants())
    new_treatment = generate_treatment_recommendation(new_variant, new_class)

    demo_gene_names = list(by_gene.keys()) or ["rpoB"]
    demo_quality, demo_qwarnings, demo_qtxt = _fetch_data_quality(demo_gene_names)
    new_treatment, demo_gated, demo_gate_warnings = _gate_treatment_on_quality(demo_quality, new_treatment)
    demo_qwarnings.extend(demo_gate_warnings)

    silent_assessment = []
    for rg in new_variant.resistance_genes:
        for mut in rg.mutations:
            if mut.is_synonymous:
                assessed = assess_silent_mutation_risk(rg.gene_name, [{
                    "codon_position": mut.position,
                    "ref_codon": mut.reference_codon,
                    "alt_codon": mut.mutant_codon,
                    "amino_acid": mut.reference_amino_acid,
                    "nucleotide_position": mut.position,
                }])
                silent_assessment.extend(assessed)

    print(f"      New Variant ID: {new_variant.variant_id}")
    print(f"      Risk: {nr.score:.1f}% [{nr.level.value}] ({nr.color.value})")
    print(f"      Resistance Class: {new_class}")
    print(f"      Resistant drugs: {', '.join(new_variant.resistant_drugs())}")
    if silent_assessment:
        print(f"      Silent mutation warnings: {len(silent_assessment)}")
        for sa in silent_assessment:
            print(f"        - Position {sa.get('codon_position')}: risk={sa.get('risk_level')} | {sa.get('note', '')}")
    if nr.factors:
        print(f"      Risk factors:")
        for f in nr.factors:
            print(f"        - {f}")
    print()

    # ── 4b. Data quality summary
    print("  [4b] Data quality summary...")
    for label, qs in demo_quality.items():
        icon = {"HIGH": "[OK]", "MODERATE": "[!!]", "LOW": "[??]", "REJECT": "[XX]"}.get(qs.confidence, "[??]")
        print(f"      {icon} {label}: {qs.raw_score:.0f}/100 ({qs.confidence})")
    if demo_gated:
        print("      *** TREATMENT GATED — manual laboratory verification required ***")
    if demo_qwarnings:
        for w in demo_qwarnings[:3]:
            print(f"      ! {w}")
    print()

    # ── 5. Reports and charts
    print("  [5] Generating reports and charts...")
    closest_comp = new_comparisons[0] if new_comparisons else None
    demo_report = _build_report(
        new_variant, nr, new_treatment, closest_comp, "NP-2024-001",
        quality_scores=demo_quality,
        quality_warnings=demo_qwarnings,
        quality_text=demo_qtxt,
        treatment_gated=demo_gated,
    )

    clinical_pdf = generate_doctor_report(demo_report)
    patient_pdf = generate_patient_report(demo_report)

    chart_paths: dict[str, str] = {}
    chart_paths["resistance_profile"] = resistance_profile_chart(new_variant)
    chart_paths["mutation_frequency"] = mutation_frequency_chart(new_variant)
    chart_paths["risk_gauge"] = risk_score_gauge(nr.score, new_variant.variant_id)
    if new_comparisons:
        chart_paths["comparison_heatmap"] = comparison_heatmap(new_comparisons)
    chart_paths["treatment_availability"] = treatment_availability_chart(new_treatment)

    print(f"      Clinical Report: {clinical_pdf}")
    print(f"      Patient Summary: {patient_pdf}")
    for name, path in chart_paths.items():
        if path:
            print(f"      Chart [{name}]: {path}")
    print()

    # ── 6. Treatment centres
    print("  [6] Treatment centres for this resistance class...")
    centres = get_treatment_centers(new_class)
    for c in centres:
        print(f"      - {c['name']} ({c['location']})")
    print()

    print("  " + "=" * 40)
    print("  Demo complete.")
    print("  " + "=" * 40)


# ── Console output ──────────────────────────────────────────────────────────

def _print_welcome():
    print()
    print("  ================================")
    print("       TBAnalytica v1.0")
    print("    TB Variant Analysis System")
    print("  ================================")
    print()


def _print_results(
    report: ClinicalReport,
    doctor_pdf: str,
    patient_pdf: str,
    chart_paths: dict,
    comparisons: list[ComparisonResult],
    treatment_rec: dict,
):
    variant = report.variant
    risk = report.risk_score

    print("  ================================")
    print("  Analysis Complete")
    print("  ================================")

    # Data Quality
    if report.data_quality:
        print()
        print("  Data Quality:")
        for label, qs in report.data_quality.items():
            icon = {"HIGH": "[OK]", "MODERATE": "[!!]", "LOW": "[??]", "REJECT": "[XX]"}.get(qs.confidence, "[??]")
            print(f"    {icon} {label}: {qs.raw_score:.0f}/100 ({qs.confidence})")
        if report.treatment_gated:
            print("    *** TREATMENT GATED — manual verification required ***")
        if report.quality_warnings:
            for w in report.quality_warnings[:3]:
                print(f"    ! {w}")
            if len(report.quality_warnings) > 3:
                print(f"    ... and {len(report.quality_warnings) - 3} more")

    print()
    print(f"  Patient ID:     {report.patient_id}")
    print(f"  Variant:        {variant.name or variant.variant_id}")
    print(f"  Lineage:        {variant.lineage or 'Unknown'}")
    print(f"  Risk Level:     {risk.level.value} (score: {risk.score:.0f}/100)")

    if variant.has_xdr_profile():
        print(f"  MDR Status:     XDR-TB CONFIRMED")
    elif variant.has_mdr_profile():
        print(f"  MDR Status:     MDR-TB CONFIRMED")
    elif variant.resistant_drugs():
        print(f"  MDR Status:     Mono/Poly-resistant ({', '.join(variant.resistant_drugs())})")
    else:
        print(f"  MDR Status:     Drug-susceptible")

    if comparisons:
        top = comparisons[0]
        print(f"  Closest Match:  {top.matched_variant_id} ({top.weighted_final_score:.1f}% similarity)")

    print(f"  Treatment:      {treatment_rec.get('resistance_class', 'N/A')} — {treatment_rec.get('regimen_string', 'N/A')}")

    warnings = treatment_rec.get("warnings", [])
    if warnings:
        print(f"  Warnings:       {len(warnings)}")
        for w in warnings[:2]:
            print(f"                  ! {w}")

    print()
    print("  Reports saved:")
    print(f"    -> {doctor_pdf}")
    print(f"    -> {patient_pdf}")
    if chart_paths:
        for name, path in chart_paths.items():
            if path:
                print(f"    -> {path}")
    print()


def _offer_open_files(doctor_pdf: str, patient_pdf: str, chart_paths: dict):
    try:
        open_now = _ask_yes_no("  Open reports now? (yes/no): ")
    except (KeyboardInterrupt, SystemExit):
        return

    if open_now:
        try:
            if os.path.exists(doctor_pdf):
                os.startfile(doctor_pdf)
            if os.path.exists(patient_pdf):
                os.startfile(patient_pdf)
            for path in chart_paths.values():
                if path and os.path.exists(path) and path.endswith(".png"):
                    os.startfile(path)
        except Exception as e:
            print(f"  Could not open files: {e}")
            print(f"  Please open them manually from the output/ directory.")


# ── Interactive menu handlers ──────────────────────────────────────────────

def _handle_known_variant():
    print()
    known = load_known_variants()
    if known:
        print("  Available variants in database:")
        for v in known:
            mdr_tag = " [MDR]" if v.has_mdr_profile() else ""
            xdr_tag = " [XDR]" if v.has_xdr_profile() else ""
            print(f"    - {v.variant_id}: {v.name} ({v.lineage}){mdr_tag}{xdr_tag}")
        print()

    variant_input = _ask_text(
        "  Enter variant ID or lineage name\n"
        "  (e.g. TB_VAR_002, Beijing, CAS, Indo-Oceanic): "
    )

    variant = _resolve_variant_input(variant_input)
    if not variant:
        print(f"\n  Could not find a variant matching '{variant_input}'.")
        print("  Please check the variant ID or lineage name and try again.")
        return

    print(f"  -> Found: {variant.variant_id} ({variant.name})")
    print()

    patient_id = _ask_text("  Enter patient ID (or press Enter to skip): ", allow_empty=True)

    print()
    analyze_known_variant(
        variant_id=variant.variant_id,
        patient_id=patient_id,
    )


def _handle_new_variant():
    print()
    has_accession = _ask_yes_no("  Do you have an NCBI accession number? (yes/no): ")

    accession = None
    sequence = None

    if has_accession:
        accession = _ask_text("  Enter accession number: ")

        # ── Early identification check ───────────────────────────────────
        identification = identify_variant(accession)
        ncbi_meta = identification.get("ncbi_metadata", {})

        # Debug: show NCBI metadata
        if ncbi_meta.get("organism") or ncbi_meta.get("strain"):
            print()
            print(f"  [DEBUG] NCBI metadata for '{accession}':")
            print(f"          organism : {ncbi_meta.get('organism', '(none)')}")
            print(f"          strain   : {ncbi_meta.get('strain', '(none)')}")
            print(f"          gene     : {ncbi_meta.get('gene', '(none)')}")
            print(f"          title    : {ncbi_meta.get('title', '(none)')}")
            print(f"          seq len  : {len(ncbi_meta.get('sequence', ''))}")

        if identification["type"] == "known":
            matched = identification["variant"]
            source_map = {
                "local_db": "local database",
                "ncbi_metadata": "NCBI metadata",
                "blast_match": "BLAST sequence match",
            }
            source_label = source_map.get(identification["source"], identification["source"])
            print()
            print(f"  ** We recognised this as {matched.variant_id} "
                  f"({matched.name}) via {source_label}.")
            if identification.get("matched_on"):
                print(f"     Matched on: {identification['matched_on']}")
            print(f"  ** Running known-variant analysis instead.")
            print()
            patient_id = _ask_text(
                "  Enter patient ID (or press Enter to skip): ", allow_empty=True,
            )
            print()
            analyze_known_variant(
                variant_id=matched.variant_id,
                patient_id=patient_id,
            )
            return
    else:
        print()
        print("  Enter protein or nucleotide sequence directly.")
        print("  (Paste the full sequence, then press Enter)")
        sequence = _ask_text("  Sequence: ")

    print()
    print("  Which gene is this?")
    print("    1. rpoB (rifampicin resistance)")
    print("    2. katG (isoniazid resistance)")
    print("    3. inhA (isoniazid resistance)")
    print("    4. gyrA (fluoroquinolone resistance)")
    print("    5. Other")
    gene_choice = _ask_choice("  Enter choice (1-5): ", range(1, 6))

    gene_map = {1: "rpoB", 2: "katG", 3: "inhA", 4: "gyrA"}
    if gene_choice == 5:
        gene_name = _ask_text("  Enter gene name: ")
    else:
        gene_name = gene_map[gene_choice]

    print()
    patient_id = _ask_text("  Enter patient ID (or press Enter to skip): ", allow_empty=True)

    print()
    analyze_new_variant(
        sequence=sequence,
        accession=accession,
        gene_name=gene_name,
        patient_id=patient_id,
    )


def _handle_compare():
    print()
    known = load_known_variants()
    if known:
        print("  Available variants in database:")
        for v in known:
            print(f"    - {v.variant_id}: {v.name}")
        print()

    first_id = _ask_text("  Enter first variant ID: ")

    variant = _resolve_variant_input(first_id)
    if not variant:
        print(f"\n  Could not find variant '{first_id}'.")
        return

    print(f"  -> Found: {variant.variant_id} ({variant.name})")
    print()

    second_input = _ask_text("  Enter second variant ID or NCBI accession: ")

    second_variant = _resolve_variant_input(second_input)
    new_accession = None
    new_sequence = None

    if second_variant:
        print(f"  -> Found known variant: {second_variant.variant_id}")
        seqs = second_variant.protein_sequences or second_variant.nucleotide_sequences
        if seqs:
            new_sequence = list(seqs.values())[0]
        else:
            print("  Warning: second variant has no stored sequences. Comparison may be limited.")
            new_sequence = "PLACEHOLDER"
    else:
        new_accession = second_input
        print(f"  -> Will fetch from NCBI: {new_accession}")

    print()
    steps = StepProgress([
        "Loading variants",
        "Fetching sequences",
        "Running comparison",
        "Analyzing differences",
        "Checking data quality",
    ])

    try:
        result = compare_variants(
            variant_id=variant.variant_id,
            new_accession=new_accession,
            new_sequence=new_sequence,
        )
    except ValueError as e:
        print(f"\n  Error: {e}")
        return

    comp = result["comparison"]

    print()
    print("  ================================")
    print("  Comparison Complete")
    print("  ================================")
    print()
    print(f"  First Variant:  {result['known_variant'].variant_id} ({result['known_variant'].name})")
    print(f"  Second Variant: {result['new_variant'].variant_id}")
    if comp:
        print(f"  Similarity:     {comp.weighted_final_score:.1f}%")
        print(f"  Confidence:     {comp.confidence_level.value}")
    print(f"  Risk (first):   {result['known_risk'].score:.1f}% [{result['known_risk'].level.value}]")
    print(f"  Risk (second):  {result['new_risk'].score:.1f}% [{result['new_risk'].level.value}]")
    print(f"  Risk delta:     {result['risk_delta']:+.1f}%")

    rec = result.get("recommendation", {})
    if rec:
        print(f"  Action:         {rec.get('confidence', 'N/A')} — flags: {rec.get('flags', [])}")

    # Data quality
    qs = result.get("quality_scores", {})
    if qs:
        print()
        print("  Data Quality:")
        for label, q in qs.items():
            icon = {"HIGH": "[OK]", "MODERATE": "[!!]", "LOW": "[??]", "REJECT": "[XX]"}.get(q.confidence, "[??]")
            print(f"    {icon} {label}: {q.raw_score:.0f}/100 ({q.confidence})")

    print()


# ── main (interactive) ────────────────────────────────────────────────────

def main():
    """Interactive entry point. Shows menu and routes to analysis functions."""

    _print_welcome()

    while True:
        print("  What would you like to do?")
        print()
        print("    1. Analyze a known TB variant")
        print("    2. Analyze a new/unknown variant")
        print("    3. Compare two variants")
        print("    4. Run demo")
        print("    5. Exit")
        print()

        choice = _ask_choice("  Enter choice (1-5): ", range(1, 6))

        try:
            if choice == 1:
                _handle_known_variant()
            elif choice == 2:
                _handle_new_variant()
            elif choice == 3:
                _handle_compare()
            elif choice == 4:
                run_demo()
            elif choice == 5:
                print()
                print("  Thank you for using TBAnalytica.")
                print()
                break

        except ValueError as e:
            log.error("Input error: %s", e)
            print(f"\n  ERROR: {e}")
            print("  Please try again.\n")
        except KeyboardInterrupt:
            print("\n\n  Analysis cancelled.")
            print()
        except Exception as e:
            log.exception("Unexpected error during analysis")
            print(f"\n  ERROR: {e}")
            print("  Check logs/ directory for details.\n")

        if choice != 5:
            try:
                again = _ask_yes_no("\n  Analyze another variant? (yes/no): ")
            except (KeyboardInterrupt, SystemExit):
                print("\n\n  Thank you for using TBAnalytica.\n")
                break
            if again:
                print()
                continue
            else:
                print()
                print("  Thank you for using TBAnalytica.")
                print()
                break


if __name__ == "__main__":
    main()
