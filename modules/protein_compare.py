"""
TBAnalytica Protein Comparison Module
Handles protein-level sequence alignment, drug binding site comparison,
and similarity scoring between TB variants using BioPython.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqUtils.ProtParam import ProteinAnalysis as BioProtAnalysis
from Bio.Data.IUPACData import protein_weights

from schema import TBVariant

from modules.api_calls import (
    fetch_uniprot_sequence,
    fetch_uniprot_features,
    fetch_alphafold_prediction,
)
from modules.variant_db import load_resistance_db


# ---------------------------------------------------------------------------
# Drug binding sites — positions critical for drug interaction
# ---------------------------------------------------------------------------

DRUG_BINDING_SITES: dict[str, dict] = {
    "rpoB": {
        "drug": "rifampicin",
        "positions": [516, 526, 531],
        "labels": {
            516: "Asp516 — rifampicin contact residue",
            526: "His526 — rifampicin binding pocket",
            531: "Ser531 — primary rifampicin anchor",
        },
        "rrdr": (507, 533),
    },
    "katG": {
        "drug": "isoniazid",
        "positions": [315],
        "labels": {
            315: "Ser315 — isoniazid activation site",
        },
        "active_site": (278, 320),
    },
    "inhA": {
        "drug": "isoniazid",
        "positions": [94, 194],
        "labels": {
            94: "Ser94 — NADH/INH binding site",
            194: "Ile194 — substrate binding loop",
        },
    },
    "gyrA": {
        "drug": "fluoroquinolones",
        "positions": [90, 91, 94],
        "labels": {
            90: "Ala90 — quinolone resistance-determining region",
            91: "Ser91 — QRDR contact residue",
            94: "Asp94 — primary fluoroquinolone anchor",
        },
        "qrdr": (88, 94),
    },
    "gyrB": {
        "drug": "fluoroquinolones",
        "positions": [538, 540],
        "labels": {
            538: "Asn538 — GyrB QRDR",
            540: "Glu540 — GyrB QRDR",
        },
    },
    "pncA": {
        "drug": "pyrazinamide",
        "positions": [49, 57, 135, 139],
        "labels": {
            49: "Asp49 — active site metal coordination",
            57: "His57 — catalytic triad",
            135: "Thr135 — substrate channel",
            139: "Val139 — core packing",
        },
    },
    "embB": {
        "drug": "ethambutol",
        "positions": [306, 406, 497],
        "labels": {
            306: "Met306 — primary ethambutol binding",
            406: "Gly406 — secondary contact",
            497: "Gln497 — transmembrane domain",
        },
    },
    "rrs": {
        "drug": "aminoglycosides",
        "positions": [1401, 1402, 1484],
        "labels": {
            1401: "A1401 — decoding site (16S rRNA)",
            1402: "C1402 — decoding site (16S rRNA)",
            1484: "G1484 — helix 44",
        },
    },
}

AMINO_ACID_CLASSES = {
    "hydrophobic": set("AILMFWVP"),
    "positive_charge": set("RKH"),
    "negative_charge": set("DE"),
    "polar": set("STNQYC"),
    "special": set("G"),
}


def _classify_amino_acid(aa: str) -> str:
    aa = aa.upper()
    for cls, members in AMINO_ACID_CLASSES.items():
        if aa in members:
            return cls
    return "special"


# ---------------------------------------------------------------------------
# 1. align_protein_sequences — global + local alignment
# ---------------------------------------------------------------------------

def align_protein_sequences(seq1: str, seq2: str) -> dict:
    """Align two protein sequences using BioPython PairwiseAligner.

    Performs both global (Needleman-Wunsch) and local (Smith-Waterman)
    alignments. Returns scores, identity, gaps, and aligned sequences.
    """
    seq1 = seq1.upper().replace("*", "")
    seq2 = seq2.upper().replace("*", "")

    if not seq1 or not seq2:
        return {
            "global": _empty_alignment_result(),
            "local": _empty_alignment_result(),
            "identity": 0.0,
            "query_length": len(seq1),
            "reference_length": len(seq2),
        }

    try:
        blosum62 = substitution_matrices.load("BLOSUM62")
    except Exception:
        blosum62 = None

    global_result = _run_alignment(seq1, seq2, "global", blosum62)
    local_result = _run_alignment(seq1, seq2, "local", blosum62)

    return {
        "global": global_result,
        "local": local_result,
        "identity": global_result["identity_percent"],
        "query_length": len(seq1),
        "reference_length": len(seq2),
    }


def _run_alignment(seq1: str, seq2: str, mode: str, matrix) -> dict:
    aligner = PairwiseAligner()
    aligner.mode = mode

    if matrix is not None:
        aligner.substitution_matrix = matrix
        aligner.open_gap_score = -11
        aligner.extend_gap_score = -1
    else:
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1, seq2)
    if not alignments:
        return _empty_alignment_result()

    best = alignments[0]

    indices = best.indices
    target_indices = indices[0]
    query_indices = indices[1]
    alignment_length = len(target_indices)

    aligned_seq1_chars = []
    aligned_seq2_chars = []
    matches = 0
    mismatches = 0
    gaps = 0

    for i in range(alignment_length):
        ti = target_indices[i]
        qi = query_indices[i]

        if ti < 0 or ti >= len(seq1):
            a = "-"
        else:
            a = seq1[ti]

        if qi < 0 or qi >= len(seq2):
            b = "-"
        else:
            b = seq2[qi]

        aligned_seq1_chars.append(a)
        aligned_seq2_chars.append(b)

        if a == "-" or b == "-":
            gaps += 1
        elif a == b:
            matches += 1
        else:
            mismatches += 1

    identity = (matches / alignment_length * 100) if alignment_length > 0 else 0.0
    aligned_seq1 = "".join(aligned_seq1_chars)
    aligned_seq2 = "".join(aligned_seq2_chars)

    return {
        "score": float(best.score),
        "identity_percent": round(identity, 2),
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "alignment_length": alignment_length,
        "aligned_query": aligned_seq1[:200] + ("..." if len(aligned_seq1) > 200 else ""),
        "aligned_reference": aligned_seq2[:200] + ("..." if len(aligned_seq2) > 200 else ""),
    }


def _empty_alignment_result() -> dict:
    return {
        "score": 0.0,
        "identity_percent": 0.0,
        "matches": 0,
        "mismatches": 0,
        "gaps": 0,
        "alignment_length": 0,
        "aligned_query": "",
        "aligned_reference": "",
    }


# ---------------------------------------------------------------------------
# 2. calculate_sequence_identity
# ---------------------------------------------------------------------------

def calculate_sequence_identity(seq1: str, seq2: str) -> float:
    """Return fraction of identical amino acid positions (0.0 to 1.0).

    For sequences of different length, aligns from the start and counts
    identity over the longer sequence length.
    """
    seq1 = seq1.upper().replace("*", "")
    seq2 = seq2.upper().replace("*", "")

    if not seq1 or not seq2:
        return 0.0

    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    total = max(len(seq1), len(seq2))
    return round(matches / total, 4)


_BINDING_KEY_MAP = {k.lower(): k for k in DRUG_BINDING_SITES}


def _get_binding_info(protein_name: str) -> dict | None:
    canonical = _BINDING_KEY_MAP.get(protein_name.lower())
    if canonical:
        return DRUG_BINDING_SITES[canonical]
    return None


# ---------------------------------------------------------------------------
# 3. compare_drug_binding_sites
# ---------------------------------------------------------------------------

def compare_drug_binding_sites(
    new_seq: str,
    known_seq: str,
    protein_name: str,
) -> dict:
    """Compare amino acids at known drug binding positions.

    Returns position-by-position comparison at binding sites for the
    given protein. These positions are clinically critical — mutations
    here directly affect drug efficacy.
    """
    new_seq = new_seq.upper().replace("*", "")
    known_seq = known_seq.upper().replace("*", "")

    binding_info = _get_binding_info(protein_name)
    if not binding_info:
        return {
            "protein": protein_name,
            "drug": "unknown",
            "binding_sites_checked": 0,
            "identical": 0,
            "mutated": 0,
            "identity_at_binding_sites": 1.0,
            "site_comparisons": [],
            "warnings": [f"No binding site data for {protein_name}"],
        }

    positions = binding_info["positions"]
    labels = binding_info.get("labels", {})
    drug = binding_info["drug"]

    site_comparisons = []
    identical = 0
    mutated = 0

    for pos in positions:
        idx = pos - 1
        new_aa = new_seq[idx] if idx < len(new_seq) else "?"
        known_aa = known_seq[idx] if idx < len(known_seq) else "?"

        is_identical = (new_aa == known_aa) and new_aa != "?"
        if is_identical:
            identical += 1
        else:
            mutated += 1

        new_class = _classify_amino_acid(new_aa) if new_aa != "?" else "unknown"
        known_class = _classify_amino_acid(known_aa) if known_aa != "?" else "unknown"
        class_change = new_class != known_class

        severity = "none"
        if not is_identical:
            if class_change:
                severity = "high"
            else:
                severity = "moderate"

        site_comparisons.append({
            "position": pos,
            "label": labels.get(pos, f"Position {pos}"),
            "new_aa": new_aa,
            "known_aa": known_aa,
            "identical": is_identical,
            "new_class": new_class,
            "known_class": known_class,
            "class_change": class_change,
            "severity": severity,
        })

    checked = identical + mutated
    identity = identical / checked if checked > 0 else 1.0

    warnings = []
    for sc in site_comparisons:
        if sc["severity"] == "high":
            warnings.append(
                f"{protein_name} position {sc['position']}: "
                f"{sc['known_aa']}->{sc['new_aa']} "
                f"({sc['known_class']}->{sc['new_class']}) — "
                f"may affect {drug} binding"
            )

    return {
        "protein": protein_name,
        "drug": drug,
        "binding_sites_checked": checked,
        "identical": identical,
        "mutated": mutated,
        "identity_at_binding_sites": round(identity, 4),
        "site_comparisons": site_comparisons,
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# 4. calculate_protein_similarity_score
# ---------------------------------------------------------------------------

def calculate_protein_similarity_score(
    new_variant: TBVariant,
    known_variant: TBVariant,
) -> float:
    """Compute weighted protein similarity between two variants.

    For each shared resistance protein:
      - Sequence identity:        40% weight
      - Drug binding site match:  40% weight
      - Overall alignment score:  20% weight

    Returns average across all compared proteins (0.0 to 1.0).
    If no protein sequences are available, falls back to mutation-based
    comparison using the resistance_genes data.
    """
    new_proteins = new_variant.protein_sequences
    known_proteins = known_variant.protein_sequences

    scores: list[float] = []

    shared_genes = set()
    for rg in new_variant.resistance_genes:
        shared_genes.add(rg.gene_name.lower())
    for rg in known_variant.resistance_genes:
        shared_genes.add(rg.gene_name.lower())

    for gene in shared_genes:
        new_seq = new_proteins.get(gene, "") or new_proteins.get(gene.upper(), "")
        known_seq = known_proteins.get(gene, "") or known_proteins.get(gene.upper(), "")

        if new_seq and known_seq:
            identity = calculate_sequence_identity(new_seq, known_seq)

            binding = compare_drug_binding_sites(new_seq, known_seq, gene)
            binding_identity = binding["identity_at_binding_sites"]

            alignment = align_protein_sequences(new_seq, known_seq)
            max_possible_score = max(len(new_seq), len(known_seq)) * 2
            norm_score = min(alignment["global"]["score"] / max_possible_score, 1.0) if max_possible_score > 0 else 0.0
            norm_score = max(norm_score, 0.0)

            weighted = 0.4 * identity + 0.4 * binding_identity + 0.2 * norm_score
            scores.append(weighted)
        else:
            score = _mutation_based_similarity(new_variant, known_variant, gene)
            scores.append(score)

    if not scores:
        return _mutation_based_similarity(new_variant, known_variant)

    return round(sum(scores) / len(scores), 4)


def _mutation_based_similarity(
    new_variant: TBVariant,
    known_variant: TBVariant,
    gene_filter: str | None = None,
) -> float:
    """Fallback similarity when protein sequences are unavailable.
    Compares mutations by position and amino acid change."""
    new_muts = set()
    known_muts = set()

    for rg in new_variant.resistance_genes:
        if gene_filter and rg.gene_name.lower() != gene_filter.lower():
            continue
        for m in rg.mutations:
            new_muts.add((rg.gene_name.lower(), m.position, m.mutant_amino_acid))

    for rg in known_variant.resistance_genes:
        if gene_filter and rg.gene_name.lower() != gene_filter.lower():
            continue
        for m in rg.mutations:
            known_muts.add((rg.gene_name.lower(), m.position, m.mutant_amino_acid))

    if not new_muts and not known_muts:
        return 1.0

    union = new_muts | known_muts
    intersection = new_muts & known_muts
    return len(intersection) / len(union) if union else 1.0


# ---------------------------------------------------------------------------
# 5. find_closest_known_variant
# ---------------------------------------------------------------------------

def find_closest_protein_match(
    new_variant: TBVariant,
    known_variants: list[TBVariant],
) -> tuple[TBVariant | None, float]:
    """Compare new variant against all known variants by protein similarity.

    Returns the closest match and its similarity score.
    If known_variants is empty, returns (None, 0.0).
    """
    if not known_variants:
        return None, 0.0

    best_variant = None
    best_score = -1.0

    for kv in known_variants:
        score = calculate_protein_similarity_score(new_variant, kv)
        if score > best_score:
            best_score = score
            best_variant = kv

    return best_variant, round(best_score, 4)


# ---------------------------------------------------------------------------
# 6. interpret_similarity
# ---------------------------------------------------------------------------

def interpret_similarity(score: float) -> dict:
    """Interpret a protein similarity score into clinical guidance.

    Thresholds:
      >0.90 — HIGH: same treatment protocol likely effective
      0.60–0.90 — MODERATE: use with caution, verify resistance profile
      <0.60 — LOW: novel variant, full workup required
    """
    if score > 0.90:
        return {
            "level": "HIGH",
            "interpretation": "HIGH similarity — same treatment protocol likely effective",
            "confidence": "high",
            "recommendation": "follow_existing_protocol",
            "requires_full_workup": False,
        }
    elif score >= 0.60:
        return {
            "level": "MODERATE",
            "interpretation": "MODERATE similarity — use with caution, verify resistance profile",
            "confidence": "moderate",
            "recommendation": "verify_before_treating",
            "requires_full_workup": False,
        }
    else:
        return {
            "level": "LOW",
            "interpretation": "LOW similarity — novel variant, full workup required",
            "confidence": "low",
            "recommendation": "full_workup_required",
            "requires_full_workup": True,
        }


# ---------------------------------------------------------------------------
# Mutation impact analysis
# ---------------------------------------------------------------------------

def analyze_mutation_impact(gene: str, position: int, ref_aa: str, alt_aa: str) -> dict:
    """Assess the structural and functional impact of a single amino acid change."""
    ref_class = _classify_amino_acid(ref_aa)
    alt_class = _classify_amino_acid(alt_aa)
    class_change = ref_class != alt_class

    ref_weight = protein_weights.get(ref_aa.upper(), 0)
    alt_weight = protein_weights.get(alt_aa.upper(), 0)
    size_change = abs(ref_weight - alt_weight)

    db = load_resistance_db()
    gene_data = db["genes"].get(gene, {})

    critical_positions = set()
    if "rrdr_region" in gene_data:
        start, end = gene_data["rrdr_region"]
        critical_positions = set(range(start, end + 1))
    elif "qrdr_region" in gene_data:
        start, end = gene_data["qrdr_region"]
        critical_positions = set(range(start, end + 1))

    binding_info = _get_binding_info(gene) or {}
    binding_positions = set(binding_info.get("positions", []))
    critical_positions |= binding_positions

    in_critical = position in critical_positions
    at_binding_site = position in binding_positions

    severity = "low"
    if at_binding_site and class_change:
        severity = "critical"
    elif at_binding_site or (in_critical and class_change):
        severity = "high"
    elif class_change or in_critical:
        severity = "moderate"

    drug = binding_info.get("drug", gene_data.get("drug_target", ""))

    return {
        "gene": gene,
        "position": position,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "short_code": f"{ref_aa}{position}{alt_aa}",
        "ref_class": ref_class,
        "alt_class": alt_class,
        "class_change": class_change,
        "size_change_daltons": round(size_change, 1),
        "in_critical_region": in_critical,
        "at_binding_site": at_binding_site,
        "drug_affected": drug,
        "predicted_severity": severity,
    }


# ---------------------------------------------------------------------------
# Protein analysis helpers (fetch + physicochemical)
# ---------------------------------------------------------------------------

def fetch_protein_analysis(gene: str) -> dict:
    """Fetch protein info from UniProt/AlphaFold and compute properties."""
    db = load_resistance_db()
    gene_data = db["genes"].get(gene, {})

    uniprot_id = gene_data.get("uniprot_id", "")
    pdb_id = gene_data.get("pdb_id")

    sequence = ""
    if uniprot_id:
        sequence = fetch_uniprot_sequence(uniprot_id) or ""

    features = []
    if uniprot_id:
        features = fetch_uniprot_features(uniprot_id)

    domains = [f for f in features if f["type"] in ("Domain", "Region")]
    active_sites_raw = [f for f in features if f["type"] == "Active site"]
    active_site_positions = []
    for site in active_sites_raw:
        loc = site.get("location", {})
        start = loc.get("start", {}).get("value")
        if start:
            active_site_positions.append(start)

    alphafold_id = None
    if uniprot_id:
        af_data = fetch_alphafold_prediction(uniprot_id)
        if af_data:
            alphafold_id = af_data.get("entryId")

    physico = compute_physicochemical_properties(sequence)

    binding_info = _get_binding_info(gene) or {}

    return {
        "uniprot_id": uniprot_id,
        "pdb_id": pdb_id,
        "alphafold_id": alphafold_id,
        "protein_name": gene_data.get("full_name", ""),
        "drug_target": gene_data.get("drug_target", ""),
        "sequence": sequence,
        "length": len(sequence),
        "domains": domains,
        "active_sites": active_site_positions,
        "binding_sites": binding_info.get("positions", []),
        "physicochemical": physico,
    }


def compute_physicochemical_properties(sequence: str) -> dict:
    """Compute physicochemical properties using BioPython ProteinAnalysis."""
    if not sequence:
        return {}

    clean = "".join(c for c in sequence.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
    if not clean:
        return {}

    analysis = BioProtAnalysis(clean)
    return {
        "molecular_weight": round(analysis.molecular_weight(), 2),
        "isoelectric_point": round(analysis.isoelectric_point(), 2),
        "instability_index": round(analysis.instability_index(), 2),
        "gravy": round(analysis.gravy(), 4),
        "aromaticity": round(analysis.aromaticity(), 4),
        "amino_acid_percent": {
            k: round(v, 4) for k, v in (
                analysis.amino_acids_percent if hasattr(analysis, "amino_acids_percent")
                else analysis.get_amino_acids_percent()
            ).items()
        },
    }


# ---------------------------------------------------------------------------
# Multi-variant protein comparison report
# ---------------------------------------------------------------------------

def compare_variants_protein_profile(
    new_variant: TBVariant,
    known_variants: list[TBVariant],
) -> dict:
    """Full protein comparison report: per-variant scores, binding site
    analysis, closest match, and clinical interpretation."""
    comparisons = []
    for kv in known_variants:
        score = calculate_protein_similarity_score(new_variant, kv)
        comparisons.append({
            "variant_id": kv.variant_id,
            "variant_name": kv.name,
            "similarity_score": score,
            "interpretation": interpret_similarity(score),
        })

    comparisons.sort(key=lambda c: c["similarity_score"], reverse=True)

    closest, closest_score = find_closest_protein_match(new_variant, known_variants)

    binding_report = {}
    for rg in new_variant.resistance_genes:
        gene = rg.gene_name.lower()
        if not _get_binding_info(gene):
            continue
        new_seq = new_variant.protein_sequences.get(gene, "")
        if closest and closest.protein_sequences.get(gene, ""):
            ref_seq = closest.protein_sequences[gene]
            binding_report[gene] = compare_drug_binding_sites(new_seq, ref_seq, gene)

    mutation_impacts = []
    for rg in new_variant.resistance_genes:
        for m in rg.mutations:
            if m.reference_amino_acid and m.mutant_amino_acid and not m.is_synonymous:
                impact = analyze_mutation_impact(
                    rg.gene_name, m.position,
                    m.reference_amino_acid, m.mutant_amino_acid,
                )
                mutation_impacts.append(impact)

    return {
        "new_variant_id": new_variant.variant_id,
        "comparisons": comparisons,
        "closest_match": {
            "variant_id": closest.variant_id if closest else None,
            "score": closest_score,
            "interpretation": interpret_similarity(closest_score),
        },
        "binding_site_analysis": binding_report,
        "mutation_impacts": mutation_impacts,
    }


# ---------------------------------------------------------------------------
# Backward-compatible aliases
# ---------------------------------------------------------------------------

def align_proteins(seq1: str, seq2: str) -> dict:
    """Legacy wrapper — returns the global alignment result dict."""
    result = align_protein_sequences(seq1, seq2)
    return result["global"]


def compare_protein_sequences(query_seq: str, reference_seq: str) -> dict:
    """Legacy wrapper — returns alignment plus per-position differences."""
    alignment = align_proteins(query_seq, reference_seq)

    differences = []
    for i, (q, r) in enumerate(zip(query_seq.upper(), reference_seq.upper())):
        if q != r:
            differences.append({
                "position": i + 1,
                "query_aa": q,
                "reference_aa": r,
                "query_class": _classify_amino_acid(q),
                "reference_class": _classify_amino_acid(r),
                "class_change": _classify_amino_acid(q) != _classify_amino_acid(r),
            })

    return {
        "alignment": alignment,
        "differences": differences,
        "total_differences": len(differences),
        "query_length": len(query_seq),
        "reference_length": len(reference_seq),
    }
