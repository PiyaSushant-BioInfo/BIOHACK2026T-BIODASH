"""
TBAnalytica Similarity Score Module
Calculates final weighted similarity scores combining protein-level,
gene-level, binding-site, and resistance-profile comparisons.
Produces RiskScore and ComparisonResult schema objects.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from schema import (
    TBVariant, Mutation, ComparisonResult, ResistanceLevel,
    ConfidenceLevel, RiskScore, RiskLevel, DrugSensitivity,
)

from modules.variant_db import load_known_variants, load_resistance_db
from modules.protein_compare import (
    calculate_sequence_identity,
    compare_drug_binding_sites,
    calculate_protein_similarity_score,
    align_proteins,
)
from modules.gene_analysis import (
    find_silent_mutations, assess_silent_mutation_risk,
    find_all_mutations, generate_gene_report,
)


# ---------------------------------------------------------------------------
# Resistance classification thresholds
# ---------------------------------------------------------------------------

_BEIJING_LINEAGE_KEYWORDS = {"beijing", "lineage 2", "east asian"}

_RISK_BASE_SCORES = {
    "SUSCEPTIBLE": (20, 30),
    "MONO": (31, 50),
    "MDR": (51, 75),
    "PRE_XDR": (76, 89),
    "XDR": (90, 100),
}


# ---------------------------------------------------------------------------
# 1. calculate_weighted_score
# ---------------------------------------------------------------------------

def calculate_weighted_score(
    protein_score: float,
    gene_score: float,
    binding_site_score: float,
    resistance_mutation_match: float,
) -> float:
    """Compute the final weighted similarity score.

    Weights:
      - Protein sequence identity:       30%
      - Gene sequence identity:           20%
      - Drug binding site comparison:     30%
      - Resistance mutation profile:      20%

    All inputs should be 0.0–1.0. Returns 0.0–1.0.
    """
    score = (
        0.30 * _clamp(protein_score)
        + 0.20 * _clamp(gene_score)
        + 0.30 * _clamp(binding_site_score)
        + 0.20 * _clamp(resistance_mutation_match)
    )
    return round(score, 4)


def _clamp(v: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, v))


# ---------------------------------------------------------------------------
# 2. calculate_risk_score
# ---------------------------------------------------------------------------

def calculate_risk_score(variant: TBVariant) -> RiskScore:
    """Score the resistance risk for a variant.

    Base score brackets:
      Drug-sensitive  20–30  LOW    green
      Mono-resistant  31–50  MODERATE yellow
      MDR-TB          51–75  HIGH   red
      Pre-XDR         76–89  HIGH   red
      XDR-TB          90–100 CRITICAL black

    Bonus factors (+5 each, capped at 100):
      - Beijing lineage
      - Precursor mutations present
      - Multiple resistance genes affected
      - High-level resistance mutations
    """
    factors: list[str] = []
    db = load_resistance_db()

    resistant_drugs = set()
    for rg in variant.resistance_genes:
        gene_data = db["genes"].get(rg.gene_name, {})
        all_mutations_db = gene_data.get("mutations", {})

        for m in rg.mutations:
            mut_code = f"{m.reference_amino_acid}{m.position}{m.mutant_amino_acid}"
            mut_info = all_mutations_db.get(mut_code, {})

            for drug in mut_info.get("drugs", []):
                resistant_drugs.add(drug.lower())

            if m.is_resistance_conferring and m.drug_affected:
                factors.append(
                    f"{rg.gene_name} {m.short_code()}: resistance to {m.drug_affected}"
                )

            if m.is_synonymous and m.one_step_away_risk:
                factors.append(
                    f"{rg.gene_name} {m.short_code()}: silent mutation, "
                    f"one-step risk for {m.one_step_away_drug}"
                )

    for drug, status in variant.drug_resistance.items():
        if status == DrugSensitivity.RESISTANT:
            resistant_drugs.add(drug.lower())

    has_rif = "rifampicin" in resistant_drugs
    has_inh = "isoniazid" in resistant_drugs
    has_fq = any(d in resistant_drugs for d in ["levofloxacin", "moxifloxacin"])
    has_inj = any(d in resistant_drugs for d in ["amikacin", "kanamycin", "capreomycin"])

    if has_rif and has_inh and has_fq and has_inj:
        category = "XDR"
    elif has_rif and has_inh and has_fq:
        category = "PRE_XDR"
    elif has_rif and has_inh:
        category = "MDR"
    elif resistant_drugs:
        category = "MONO"
    else:
        category = "SUSCEPTIBLE"

    lo, hi = _RISK_BASE_SCORES[category]
    base = _scale_within_range(variant, db, lo, hi)

    bonus = 0
    lineage_lower = variant.lineage.lower()
    if any(kw in lineage_lower for kw in _BEIJING_LINEAGE_KEYWORDS):
        bonus += 5
        factors.append("Beijing/Lineage 2: higher transmissibility in Nepal")

    precursors = [
        m for m in variant.all_mutations()
        if m.is_synonymous and m.one_step_away_risk
    ]
    if precursors:
        bonus += 5
        factors.append(f"{len(precursors)} precursor mutation(s) detected")

    affected_genes = {rg.gene_name for rg in variant.resistance_genes if rg.mutations}
    if len(affected_genes) >= 3:
        bonus += 5
        factors.append(f"Mutations across {len(affected_genes)} resistance genes")

    high_level = [
        m for m in variant.all_mutations()
        if m.is_resistance_conferring and m.resistance_level == ResistanceLevel.HIGH
    ]
    if high_level:
        bonus += 5
        factors.append(f"{len(high_level)} high-level resistance mutation(s)")

    final = min(base + bonus, 100.0)

    return RiskScore(
        variant_id=variant.variant_id,
        score=round(final, 2),
        factors=factors,
    )


def _scale_within_range(
    variant: TBVariant, db: dict, lo: float, hi: float,
) -> float:
    """Position the score within its category range based on mutation severity."""
    total_weight = 0.0
    max_weight = 0.0

    for rg in variant.resistance_genes:
        gene_data = db["genes"].get(rg.gene_name, {})
        mutations_db = gene_data.get("mutations", {})

        for m in rg.mutations:
            mut_code = f"{m.reference_amino_acid}{m.position}{m.mutant_amino_acid}"
            mut_info = mutations_db.get(mut_code, {})

            severity = {"high": 1.0, "moderate": 0.6, "low": 0.3}
            w = severity.get(mut_info.get("resistance", "low"), 0.1)
            freq = mut_info.get("frequency", 0.1)
            total_weight += w * (1.0 + freq)

            if m.is_synonymous:
                for _, prec in gene_data.get("silent_precursors", {}).items():
                    if prec.get("position") == m.position:
                        prec_w = {"high": 0.4, "moderate": 0.2, "low": 0.1}
                        total_weight += prec_w.get(prec.get("risk", "low"), 0.05)

            max_weight += 2.0

    ratio = (total_weight / max_weight) if max_weight > 0 else 0.5
    return lo + ratio * (hi - lo)


# ---------------------------------------------------------------------------
# 3. compare_resistance_profiles
# ---------------------------------------------------------------------------

def compare_resistance_profiles(
    new_variant: TBVariant,
    known_variant: TBVariant,
) -> float:
    """Compare drug resistance pattern similarity between two variants.

    Checks every drug present in either variant's resistance profile and
    scores how many have the same classification (resistant/sensitive).
    Returns 0.0–1.0.
    """
    new_dr = {d.lower(): s for d, s in new_variant.drug_resistance.items()}
    known_dr = {d.lower(): s for d, s in known_variant.drug_resistance.items()}

    new_mut_drugs = _drugs_from_mutations(new_variant)
    known_mut_drugs = _drugs_from_mutations(known_variant)

    for d in new_mut_drugs:
        new_dr.setdefault(d, DrugSensitivity.RESISTANT)
    for d in known_mut_drugs:
        known_dr.setdefault(d, DrugSensitivity.RESISTANT)

    all_drugs = set(new_dr.keys()) | set(known_dr.keys())
    if not all_drugs:
        return 1.0

    matches = 0
    for drug in all_drugs:
        new_status = new_dr.get(drug, DrugSensitivity.SENSITIVE)
        known_status = known_dr.get(drug, DrugSensitivity.SENSITIVE)
        if new_status == known_status:
            matches += 1

    return round(matches / len(all_drugs), 4)


def _drugs_from_mutations(variant: TBVariant) -> set[str]:
    """Extract drug names from resistance-conferring mutations."""
    db = load_resistance_db()
    drugs = set()
    for rg in variant.resistance_genes:
        gene_data = db["genes"].get(rg.gene_name, {})
        mutations_db = gene_data.get("mutations", {})
        for m in rg.mutations:
            if m.drug_affected:
                for d in m.drug_affected.split(","):
                    drugs.add(d.strip().lower())
            mut_code = f"{m.reference_amino_acid}{m.position}{m.mutant_amino_acid}"
            for d in mutations_db.get(mut_code, {}).get("drugs", []):
                drugs.add(d.lower())
    return drugs


# ---------------------------------------------------------------------------
# 4. full_similarity_analysis
# ---------------------------------------------------------------------------

def full_similarity_analysis(
    new_variant: TBVariant,
    known_variants: list[TBVariant],
) -> list[ComparisonResult]:
    """Run complete similarity analysis against all known variants.

    For each known variant computes:
      - Protein similarity (sequence identity + binding sites)
      - Gene similarity (nucleotide identity)
      - Resistance profile match
      - Weighted final score
    Returns top 5 matches sorted by weighted_final_score descending.
    """
    results: list[ComparisonResult] = []

    new_mutations = new_variant.all_mutations()
    silent = [m for m in new_mutations if m.is_synonymous]
    resistance = [m for m in new_mutations if m.is_resistance_conferring]

    for known in known_variants:
        if known.variant_id == new_variant.variant_id:
            continue

        protein_sim = _compute_protein_similarity(new_variant, known)
        gene_sim = _compute_gene_similarity(new_variant, known)
        binding_sim = _compute_binding_site_similarity(new_variant, known)
        profile_match = compare_resistance_profiles(new_variant, known)

        weighted = calculate_weighted_score(
            protein_score=protein_sim,
            gene_score=gene_sim,
            binding_site_score=binding_sim,
            resistance_mutation_match=profile_match,
        )

        known_mutations = known.all_mutations()
        overlap = _compute_mutation_overlap(new_mutations, known_mutations)
        novel = overlap["unique_to_query"]

        res_class = _classify_resistance(new_variant)

        results.append(ComparisonResult(
            new_variant_id=new_variant.variant_id,
            matched_variant_id=known.variant_id,
            protein_similarity_score=round(protein_sim * 100, 2),
            gene_similarity_score=round(gene_sim * 100, 2),
            weighted_final_score=round(weighted * 100, 2),
            gene_changes=new_mutations,
            silent_mutations=silent,
            resistance_mutations=resistance,
            novel_mutations=novel,
            treatment_recommendation=res_class,
        ))

    results.sort(key=lambda r: r.weighted_final_score, reverse=True)
    return results[:5]


def _compute_protein_similarity(new: TBVariant, known: TBVariant) -> float:
    """Protein similarity using sequence identity across shared proteins."""
    shared_genes = set()
    for rg in new.resistance_genes:
        shared_genes.add(rg.gene_name)
    for rg in known.resistance_genes:
        shared_genes.add(rg.gene_name)

    scores = []
    for gene in shared_genes:
        new_seq = new.protein_sequences.get(gene, "") or new.protein_sequences.get(gene.lower(), "")
        known_seq = known.protein_sequences.get(gene, "") or known.protein_sequences.get(gene.lower(), "")
        if new_seq and known_seq:
            scores.append(calculate_sequence_identity(new_seq, known_seq))

    if scores:
        return sum(scores) / len(scores)
    return _mutation_jaccard(new, known)


def _compute_gene_similarity(new: TBVariant, known: TBVariant) -> float:
    """Gene similarity using nucleotide sequence identity."""
    shared_genes = set()
    for rg in new.resistance_genes:
        shared_genes.add(rg.gene_name)
    for rg in known.resistance_genes:
        shared_genes.add(rg.gene_name)

    scores = []
    for gene in shared_genes:
        new_seq = new.nucleotide_sequences.get(gene, "") or new.nucleotide_sequences.get(gene.lower(), "")
        known_seq = known.nucleotide_sequences.get(gene, "") or known.nucleotide_sequences.get(gene.lower(), "")
        if new_seq and known_seq:
            min_len = min(len(new_seq), len(known_seq))
            if min_len > 0:
                matches = sum(1 for a, b in zip(new_seq.upper(), known_seq.upper()) if a == b)
                scores.append(matches / max(len(new_seq), len(known_seq)))

    if scores:
        return sum(scores) / len(scores)
    return _mutation_jaccard(new, known)


def _compute_binding_site_similarity(new: TBVariant, known: TBVariant) -> float:
    """Binding site similarity across shared resistance proteins."""
    shared_genes = set()
    for rg in new.resistance_genes:
        shared_genes.add(rg.gene_name)
    for rg in known.resistance_genes:
        shared_genes.add(rg.gene_name)

    scores = []
    for gene in shared_genes:
        new_seq = new.protein_sequences.get(gene, "") or new.protein_sequences.get(gene.lower(), "")
        known_seq = known.protein_sequences.get(gene, "") or known.protein_sequences.get(gene.lower(), "")
        if new_seq and known_seq:
            result = compare_drug_binding_sites(new_seq, known_seq, gene)
            scores.append(result["identity_at_binding_sites"])

    if scores:
        return sum(scores) / len(scores)

    return compare_resistance_profiles(new, known)


def _mutation_jaccard(new: TBVariant, known: TBVariant) -> float:
    """Mutation-level Jaccard similarity as fallback."""
    new_set = {
        (rg.gene_name, m.position, m.mutant_amino_acid)
        for rg in new.resistance_genes for m in rg.mutations
    }
    known_set = {
        (rg.gene_name, m.position, m.mutant_amino_acid)
        for rg in known.resistance_genes for m in rg.mutations
    }
    if not new_set and not known_set:
        return 1.0
    union = new_set | known_set
    intersection = new_set & known_set
    return len(intersection) / len(union) if union else 1.0


def _classify_resistance(variant: TBVariant) -> str:
    """Return a resistance classification string."""
    db = load_resistance_db()
    resistant_drugs: set[str] = set()

    for rg in variant.resistance_genes:
        gene_data = db["genes"].get(rg.gene_name, {})
        mutations_db = gene_data.get("mutations", {})
        for m in rg.mutations:
            mut_code = f"{m.reference_amino_acid}{m.position}{m.mutant_amino_acid}"
            for d in mutations_db.get(mut_code, {}).get("drugs", []):
                resistant_drugs.add(d.lower())

    for drug, status in variant.drug_resistance.items():
        if status == DrugSensitivity.RESISTANT:
            resistant_drugs.add(drug.lower())

    has_rif = "rifampicin" in resistant_drugs
    has_inh = "isoniazid" in resistant_drugs
    has_fq = any(d in resistant_drugs for d in ["levofloxacin", "moxifloxacin"])
    has_inj = any(d in resistant_drugs for d in ["amikacin", "kanamycin", "capreomycin"])

    if has_rif and has_inh and has_fq and has_inj:
        return "XDR"
    if has_rif and has_inh and has_fq:
        return "Pre-XDR"
    if has_rif and has_inh:
        return "MDR"
    if resistant_drugs:
        return "MONO-RESISTANT"
    return "SUSCEPTIBLE"


# ---------------------------------------------------------------------------
# 5. interpret_and_recommend
# ---------------------------------------------------------------------------

def interpret_and_recommend(comparison: ComparisonResult) -> dict:
    """Generate clinical recommendation from a ComparisonResult.

    Score thresholds (weighted_final_score is 0–100):
      >90  — same treatment as matched variant, HIGH confidence
      75–90 — same treatment with monitoring, MODERATE confidence
      60–75 — modified treatment, LOW confidence, verify with AST
      <60  — novel variant, do not assume, full susceptibility testing
    """
    score = comparison.weighted_final_score
    resistance_class = comparison.treatment_recommendation
    novel_count = len(comparison.novel_mutations)

    warnings: list[str] = []

    if score > 90:
        confidence = "HIGH"
        recommendation = (
            f"Follow treatment protocol for matched variant "
            f"{comparison.matched_variant_id} ({resistance_class}). "
            f"High similarity supports direct protocol adoption."
        )
        flags = ["protocol_match"]
    elif score >= 75:
        confidence = "MODERATE"
        recommendation = (
            f"Treatment protocol for {comparison.matched_variant_id} "
            f"({resistance_class}) may be applied with close monitoring. "
            f"Schedule follow-up culture at 2 months."
        )
        flags = ["monitor_closely", "follow_up_culture"]
        if novel_count > 0:
            warnings.append(
                f"{novel_count} novel mutation(s) detected — monitor for "
                f"treatment failure"
            )
    elif score >= 60:
        confidence = "LOW"
        recommendation = (
            f"Modified treatment based on {comparison.matched_variant_id} "
            f"({resistance_class}). Significant differences detected — "
            f"confirm with phenotypic antimicrobial susceptibility testing (AST)."
        )
        flags = ["modified_treatment", "ast_required", "specialist_review"]
        warnings.append("Partial match only — phenotypic AST strongly recommended")
        if novel_count > 0:
            warnings.append(
                f"{novel_count} novel mutation(s) with unknown resistance impact"
            )
    else:
        confidence = "LOW"
        recommendation = (
            f"Novel variant — insufficient similarity to any known variant. "
            f"Do NOT assume susceptibility pattern. Full phenotypic susceptibility "
            f"testing required before initiating targeted therapy. "
            f"Start empiric regimen per WHO guidelines until AST results available."
        )
        flags = ["novel_variant", "full_ast_required", "empiric_therapy", "specialist_review"]
        warnings.append("Novel variant: no reliable match in database")
        if comparison.resistance_mutations:
            warnings.append(
                f"{len(comparison.resistance_mutations)} known resistance "
                f"mutation(s) present despite low overall similarity"
            )

    if comparison.silent_mutations:
        precursors = [m for m in comparison.silent_mutations if m.one_step_away_risk]
        if precursors:
            warnings.append(
                f"{len(precursors)} silent precursor mutation(s) — "
                f"monitor for resistance emergence"
            )
            if "monitor_precursors" not in flags:
                flags.append("monitor_precursors")

    return {
        "matched_variant_id": comparison.matched_variant_id,
        "weighted_score": score,
        "resistance_class": resistance_class,
        "confidence": confidence,
        "recommendation": recommendation,
        "flags": flags,
        "warnings": warnings,
        "novel_mutations_count": novel_count,
        "resistance_mutations_count": len(comparison.resistance_mutations),
        "silent_mutations_count": len(comparison.silent_mutations),
    }


# ---------------------------------------------------------------------------
# Mutation overlap helper
# ---------------------------------------------------------------------------

def _compute_mutation_overlap(
    query_mutations: list[Mutation],
    ref_mutations: list[Mutation],
) -> dict:
    query_set = {(m.position, m.mutant_amino_acid) for m in query_mutations}
    ref_set = {(m.position, m.mutant_amino_acid) for m in ref_mutations}

    shared_keys = query_set & ref_set
    unique_query_keys = query_set - ref_set
    unique_ref_keys = ref_set - query_set

    shared = [m for m in query_mutations if (m.position, m.mutant_amino_acid) in shared_keys]
    unique_to_query = [m for m in query_mutations if (m.position, m.mutant_amino_acid) in unique_query_keys]
    unique_to_ref = [m for m in ref_mutations if (m.position, m.mutant_amino_acid) in unique_ref_keys]

    total = len(query_set | ref_set)
    jaccard = len(shared_keys) / total if total > 0 else 0.0

    return {
        "shared": shared,
        "unique_to_query": unique_to_query,
        "unique_to_reference": unique_to_ref,
        "jaccard_similarity": round(jaccard, 4),
        "overlap_count": len(shared_keys),
    }


# ---------------------------------------------------------------------------
# Backward-compatible public API
# ---------------------------------------------------------------------------

def compute_mutation_overlap(
    query_mutations: list[Mutation],
    ref_mutations: list[Mutation],
) -> dict:
    """Legacy public alias."""
    return _compute_mutation_overlap(query_mutations, ref_mutations)


def compute_resistance_risk_score(variant: TBVariant) -> float:
    """Legacy: return numeric risk score only."""
    risk = calculate_risk_score(variant)
    return risk.score


def compute_risk_score(variant: TBVariant) -> RiskScore:
    """Legacy alias for calculate_risk_score."""
    return calculate_risk_score(variant)


def predict_resistance_level(risk_score: float, all_mutations: list[Mutation]) -> str:
    """Legacy: classify resistance from score and mutations."""
    db = load_resistance_db()
    resistant_drugs: set[str] = set()
    for m in all_mutations:
        for gene_name, gene_data in db["genes"].items():
            mut_code = f"{m.reference_amino_acid}{m.position}{m.mutant_amino_acid}"
            for d in gene_data.get("mutations", {}).get(mut_code, {}).get("drugs", []):
                resistant_drugs.add(d.lower())

    has_rif = "rifampicin" in resistant_drugs
    has_inh = "isoniazid" in resistant_drugs
    has_fq = any(d in resistant_drugs for d in ["levofloxacin", "moxifloxacin"])
    has_inj = any(d in resistant_drugs for d in ["amikacin", "kanamycin", "capreomycin"])

    if has_rif and has_inh and has_fq and has_inj:
        return "XDR"
    if has_rif and has_inh and has_fq:
        return "Pre-XDR"
    if has_rif and has_inh:
        return "MDR"
    if risk_score > 60:
        return "HIGH"
    if risk_score > 35:
        return "MODERATE"
    if risk_score > 15:
        return "LOW"
    return "SUSCEPTIBLE"


def compare_variant_to_known(query: TBVariant) -> list[ComparisonResult]:
    """Legacy: compare against all known variants in the database."""
    known = load_known_variants()
    return full_similarity_analysis(query, known)


def find_closest_known_variant(query: TBVariant) -> ComparisonResult:
    """Legacy: find the single best match."""
    known = load_known_variants()
    results = full_similarity_analysis(query, known)
    if not results:
        risk = compute_resistance_risk_score(query)
        all_muts = query.all_mutations()
        return ComparisonResult(
            new_variant_id=query.variant_id,
            matched_variant_id="none",
            treatment_recommendation=predict_resistance_level(risk, all_muts),
        )
    return results[0]
