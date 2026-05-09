"""
TBAnalytica Gene Analysis Module
Handles nucleotide-level gene comparison, codon translation,
silent mutation detection, one-step-away resistance precursor analysis,
and comprehensive gene-level mutation reporting.
"""

from pathlib import Path
from collections import Counter
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from schema import Mutation, ResistanceLevel, TBVariant

from Bio.Align import PairwiseAligner

from modules.api_calls import fetch_ncbi_gene_sequence, search_ncbi_gene
from modules.variant_db import load_resistance_db, get_silent_precursors


# ---------------------------------------------------------------------------
# Standard genetic code — all 64 codons
# ---------------------------------------------------------------------------

CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}

NUCLEOTIDES = "ACGT"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _translate_codon(codon: str) -> str:
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3:
        return "X"
    return CODON_TABLE.get(codon, "X")


def _single_nt_neighbours(codon: str) -> list[str]:
    """Return all codons reachable by a single nucleotide substitution."""
    codon = codon.upper()
    neighbours = []
    for i in range(3):
        for nt in NUCLEOTIDES:
            if nt != codon[i]:
                neighbours.append(codon[:i] + nt + codon[i + 1:])
    return neighbours


def _get_critical_region(gene_data: dict) -> set[int]:
    positions = set()
    if "rrdr_region" in gene_data:
        s, e = gene_data["rrdr_region"]
        positions = set(range(s, e + 1))
    elif "qrdr_region" in gene_data:
        s, e = gene_data["qrdr_region"]
        positions = set(range(s, e + 1))
    return positions


# ---------------------------------------------------------------------------
# 1. translate_sequence
# ---------------------------------------------------------------------------

def translate_sequence(nucleotide_seq: str) -> str:
    """Translate a nucleotide sequence to amino acids using the standard
    genetic code. Handles all 64 codons. Stop codons are rendered as '*'.
    Incomplete trailing codons are ignored.
    """
    seq = nucleotide_seq.upper().replace("U", "T")
    amino_acids = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if len(codon) < 3:
            break
        aa = CODON_TABLE.get(codon, "X")
        amino_acids.append(aa)
    return "".join(amino_acids)


# ---------------------------------------------------------------------------
# 2. find_all_mutations
# ---------------------------------------------------------------------------

def find_all_mutations(
    new_seq: str,
    reference_seq: str,
    gene_name: str,
) -> list[Mutation]:
    """Align nucleotide sequences and find every codon-level difference.

    For each differing codon position:
      - Translates both codons
      - Determines synonymous vs non-synonymous
      - Checks the resistance database
      - Flags one-step-away precursors
    Returns a list of Mutation schema objects.
    """
    new_seq = new_seq.upper().replace("U", "T")
    reference_seq = reference_seq.upper().replace("U", "T")

    db = load_resistance_db()
    gene_data = db["genes"].get(gene_name, {})
    known_mutations = gene_data.get("mutations", {})

    min_len = min(len(new_seq), len(reference_seq))
    mutations: list[Mutation] = []

    for i in range(0, min_len - 2, 3):
        ref_codon = reference_seq[i:i + 3]
        new_codon = new_seq[i:i + 3]

        if ref_codon == new_codon:
            continue

        ref_aa = _translate_codon(ref_codon)
        new_aa = _translate_codon(new_codon)
        codon_pos = (i // 3) + 1

        is_syn = (ref_aa == new_aa) and ref_aa != "X"

        mut_code = f"{ref_aa}{codon_pos}{new_aa}"
        db_info = known_mutations.get(mut_code, {})
        is_resistance = bool(db_info.get("drugs")) and not is_syn
        drug_affected = ", ".join(db_info.get("drugs", []))
        res_level = ResistanceLevel.HIGH if db_info.get("resistance") == "high" else ResistanceLevel.MODERATE

        one_step = False
        one_step_drug = None
        if is_syn:
            step_info = check_one_step_away(gene_name, codon_pos, new_codon)
            one_step = step_info["is_precursor"]
            one_step_drug = step_info.get("target_drug")

        mutations.append(Mutation(
            position=codon_pos,
            reference_codon=ref_codon,
            mutant_codon=new_codon,
            reference_amino_acid=ref_aa,
            mutant_amino_acid=new_aa,
            is_synonymous=is_syn,
            is_resistance_conferring=is_resistance,
            drug_affected=drug_affected,
            resistance_level=res_level if is_resistance else ResistanceLevel.LOW,
            one_step_away_risk=one_step,
            one_step_away_drug=one_step_drug,
        ))

    return mutations


# ---------------------------------------------------------------------------
# 3. classify_mutation
# ---------------------------------------------------------------------------

def classify_mutation(mutation: Mutation) -> Mutation:
    """Classify a mutation and return an updated copy.

    Categories:
      - synonymous (silent): same amino acid, different codon
      - nonsense: premature stop codon
      - missense: different amino acid
      - frameshift: flagged externally (indel detection)
    """
    ref_aa = mutation.reference_amino_acid
    mut_aa = mutation.mutant_amino_acid

    if mut_aa == "*" and ref_aa != "*":
        mutation.is_synonymous = False
        mutation.is_resistance_conferring = True
        if not mutation.drug_affected:
            db = load_resistance_db()
            gene_data = db["genes"].get("", {})
            mutation.drug_affected = gene_data.get("drug_target", "")
        return mutation

    if ref_aa and mut_aa and ref_aa == mut_aa:
        mutation.is_synonymous = True
        mutation.is_resistance_conferring = False
        return mutation

    if ref_aa and mut_aa and ref_aa != mut_aa:
        mutation.is_synonymous = False
        return mutation

    return mutation


def get_mutation_type_label(mutation: Mutation) -> str:
    """Return a human-readable classification label."""
    if mutation.mutant_amino_acid == "*" and mutation.reference_amino_acid != "*":
        return "nonsense"
    if mutation.is_synonymous:
        return "silent"
    if mutation.reference_amino_acid and mutation.mutant_amino_acid:
        return "missense"
    return "unknown"


# ---------------------------------------------------------------------------
# 4. check_resistance_database
# ---------------------------------------------------------------------------

def check_resistance_database(
    gene_name: str,
    position: int,
    amino_acid_change: str,
) -> dict:
    """Check the resistance database for a specific mutation.

    Looks up the WHO TB mutation catalogue data (via resistance_mutations.json)
    and returns whether this is a known resistance mutation, which drug is
    affected, and the resistance level.

    Args:
        gene_name: Gene symbol (e.g. "rpoB")
        position: Codon position
        amino_acid_change: Mutation code like "S531L" or just "L" (alt AA)
    """
    db = load_resistance_db()
    gene_data = db["genes"].get(gene_name, {})
    known = gene_data.get("mutations", {})

    if len(amino_acid_change) > 2 and amino_acid_change[0].isalpha() and amino_acid_change[-1].isalpha():
        mut_code = amino_acid_change
    else:
        for code, info in known.items():
            if code[1:-1].lstrip("-").isdigit() and int(code[1:-1]) == position:
                if len(amino_acid_change) == 1 and code[-1] == amino_acid_change:
                    mut_code = code
                    break
        else:
            mut_code = amino_acid_change

    db_entry = known.get(mut_code, {})

    critical_region = _get_critical_region(gene_data)

    return {
        "gene": gene_name,
        "position": position,
        "mutation_code": mut_code,
        "is_known_resistance": bool(db_entry.get("drugs")),
        "drug_affected": db_entry.get("drugs", []),
        "resistance_level": db_entry.get("resistance", "unknown"),
        "frequency": db_entry.get("frequency", 0.0),
        "mic_fold_change": db_entry.get("mic_fold_change", ""),
        "in_critical_region": position in critical_region,
        "drug_target": gene_data.get("drug_target", ""),
    }


# ---------------------------------------------------------------------------
# 5. check_one_step_away
# ---------------------------------------------------------------------------

def check_one_step_away(
    gene_name: str,
    position: int,
    current_codon: str,
) -> dict:
    """Check if a single nucleotide change from the current codon would
    produce a known resistance mutation at this position.

    Generates all 9 possible single-nt neighbours of the codon, translates
    each, and checks if the resulting amino acid change is in the resistance
    database for this gene.
    """
    current_codon = current_codon.upper().replace("U", "T")
    current_aa = _translate_codon(current_codon)

    db = load_resistance_db()
    gene_data = db["genes"].get(gene_name, {})
    known = gene_data.get("mutations", {})
    drug_target = gene_data.get("drug_target", "")

    resistance_at_pos: dict[str, dict] = {}
    for code, info in known.items():
        if len(code) > 2 and code[0].isalpha() and code[-1].isalpha():
            pos_str = code[1:-1]
            if pos_str.lstrip("-").isdigit() and int(pos_str) == position:
                target_aa = code[-1]
                resistance_at_pos[target_aa] = {
                    "mutation_code": code,
                    "drugs": info.get("drugs", []),
                    "resistance": info.get("resistance", "unknown"),
                }

    if not resistance_at_pos:
        precursors = gene_data.get("silent_precursors", {})
        for _, pinfo in precursors.items():
            if pinfo.get("position") == position:
                return {
                    "is_precursor": True,
                    "position": position,
                    "current_codon": current_codon,
                    "current_aa": current_aa,
                    "risk_level": pinfo.get("risk", "moderate"),
                    "target_mutations": [],
                    "target_drug": drug_target,
                    "note": pinfo.get("note", "Known silent precursor position"),
                }

        return {
            "is_precursor": False,
            "position": position,
            "current_codon": current_codon,
            "current_aa": current_aa,
            "risk_level": "none",
            "target_mutations": [],
            "target_drug": None,
            "note": "",
        }

    neighbours = _single_nt_neighbours(current_codon)
    reachable_targets = []

    for neighbour_codon in neighbours:
        neighbour_aa = _translate_codon(neighbour_codon)
        if neighbour_aa in resistance_at_pos and neighbour_aa != current_aa:
            target = resistance_at_pos[neighbour_aa]
            reachable_targets.append({
                "target_codon": neighbour_codon,
                "target_aa": neighbour_aa,
                "mutation_code": target["mutation_code"],
                "drugs": target["drugs"],
                "resistance": target["resistance"],
                "nucleotide_change": _describe_nt_change(current_codon, neighbour_codon),
            })

    if reachable_targets:
        max_risk = max(
            ("high" if t["resistance"] == "high" else "moderate" for t in reachable_targets),
            key=lambda r: {"high": 2, "moderate": 1}.get(r, 0),
        )
        drugs = list({d for t in reachable_targets for d in t["drugs"]})
        return {
            "is_precursor": True,
            "position": position,
            "current_codon": current_codon,
            "current_aa": current_aa,
            "risk_level": max_risk,
            "target_mutations": reachable_targets,
            "target_drug": ", ".join(drugs) if drugs else drug_target,
            "note": f"One SNP away from {len(reachable_targets)} resistance mutation(s)",
        }

    return {
        "is_precursor": False,
        "position": position,
        "current_codon": current_codon,
        "current_aa": current_aa,
        "risk_level": "none",
        "target_mutations": [],
        "target_drug": None,
        "note": "",
    }


def _describe_nt_change(from_codon: str, to_codon: str) -> str:
    changes = []
    for i, (f, t) in enumerate(zip(from_codon, to_codon)):
        if f != t:
            changes.append(f"{f}{i + 1}{t}")
    return "; ".join(changes)


# ---------------------------------------------------------------------------
# 6. generate_gene_report
# ---------------------------------------------------------------------------

def generate_gene_report(
    mutations: list[Mutation],
    gene_name: str,
) -> dict:
    """Categorize and summarize all mutations in a gene.

    Categories:
      - resistance: known drug-resistance mutations
      - silent_precursor: synonymous mutations that are one step away
      - silent_benign: synonymous mutations with no precursor risk
      - novel: non-synonymous mutations not in the resistance database
      - nonsense: premature stop codons
    """
    db = load_resistance_db()
    gene_data = db["genes"].get(gene_name, {})
    drug_target = gene_data.get("drug_target", "")
    critical_region = _get_critical_region(gene_data)

    resistance = []
    silent_precursor = []
    silent_benign = []
    novel = []
    nonsense = []

    for m in mutations:
        label = get_mutation_type_label(m)

        if label == "nonsense":
            nonsense.append(m)
        elif m.is_resistance_conferring:
            resistance.append(m)
        elif m.is_synonymous and m.one_step_away_risk:
            silent_precursor.append(m)
        elif m.is_synonymous:
            silent_benign.append(m)
        else:
            novel.append(m)

    risk_level = "none"
    if resistance or nonsense:
        risk_level = "high"
    elif silent_precursor:
        risk_level = "moderate"
    elif novel:
        in_critical = any(m.position in critical_region for m in novel)
        risk_level = "moderate" if in_critical else "low"
    elif silent_benign:
        risk_level = "low"

    return {
        "gene_name": gene_name,
        "drug_target": drug_target,
        "total_mutations": len(mutations),
        "resistance": [_mutation_summary(m) for m in resistance],
        "silent_precursor": [_mutation_summary(m) for m in silent_precursor],
        "silent_benign": [_mutation_summary(m) for m in silent_benign],
        "novel": [_mutation_summary(m) for m in novel],
        "nonsense": [_mutation_summary(m) for m in nonsense],
        "resistance_count": len(resistance),
        "silent_precursor_count": len(silent_precursor),
        "silent_benign_count": len(silent_benign),
        "novel_count": len(novel),
        "nonsense_count": len(nonsense),
        "overall_risk": risk_level,
    }


def _mutation_summary(m: Mutation) -> dict:
    return {
        "position": m.position,
        "short_code": m.short_code(),
        "ref_codon": m.reference_codon,
        "mut_codon": m.mutant_codon,
        "ref_aa": m.reference_amino_acid,
        "mut_aa": m.mutant_amino_acid,
        "type": get_mutation_type_label(m),
        "drug_affected": m.drug_affected,
        "resistance_level": m.resistance_level.value,
        "one_step_risk": m.one_step_away_risk,
        "one_step_drug": m.one_step_away_drug,
    }


# ---------------------------------------------------------------------------
# 7. full_gene_analysis
# ---------------------------------------------------------------------------

def full_gene_analysis(
    new_variant: TBVariant,
    reference_variant: TBVariant,
) -> dict:
    """Run complete gene-level analysis across ALL resistance genes.

    Compares nucleotide sequences for every gene present in both variants,
    finds all mutations, classifies them, and produces a comprehensive
    cross-gene report.
    """
    gene_reports: list[dict] = []
    all_mutations: list[Mutation] = []
    total_resistance = 0
    total_silent = 0
    total_precursor = 0
    total_novel = 0
    total_nonsense = 0
    warnings: list[str] = []

    genes_in_new = {rg.gene_name for rg in new_variant.resistance_genes}
    genes_in_ref = {rg.gene_name for rg in reference_variant.resistance_genes}
    all_gene_names = genes_in_new | genes_in_ref

    new_nuc = new_variant.nucleotide_sequences
    ref_nuc = reference_variant.nucleotide_sequences

    for gene_name in sorted(all_gene_names):
        new_seq = new_nuc.get(gene_name, "") or new_nuc.get(gene_name.lower(), "")
        ref_seq = ref_nuc.get(gene_name, "") or ref_nuc.get(gene_name.lower(), "")

        if new_seq and ref_seq:
            mutations = find_all_mutations(new_seq, ref_seq, gene_name)
        else:
            mutations = _mutations_from_resistance_genes(new_variant, reference_variant, gene_name)
            if not new_seq and not ref_seq:
                warnings.append(f"{gene_name}: no nucleotide sequences available, using mutation-level comparison")

        all_mutations.extend(mutations)
        report = generate_gene_report(mutations, gene_name)
        gene_reports.append(report)

        total_resistance += report["resistance_count"]
        total_silent += report["silent_benign_count"]
        total_precursor += report["silent_precursor_count"]
        total_novel += report["novel_count"]
        total_nonsense += report["nonsense_count"]

    overall_risk = "none"
    if total_resistance > 0 or total_nonsense > 0:
        overall_risk = "high"
    elif total_precursor > 0:
        overall_risk = "moderate"
    elif total_novel > 0:
        overall_risk = "low-moderate"
    elif total_silent > 0:
        overall_risk = "low"

    drugs_affected = list({
        m.drug_affected for m in all_mutations
        if m.is_resistance_conferring and m.drug_affected
    })

    precursor_drugs = list({
        m.one_step_away_drug for m in all_mutations
        if m.one_step_away_risk and m.one_step_away_drug
    })

    return {
        "new_variant_id": new_variant.variant_id,
        "reference_variant_id": reference_variant.variant_id,
        "genes_analyzed": len(gene_reports),
        "gene_reports": gene_reports,
        "summary": {
            "total_mutations": len(all_mutations),
            "resistance_mutations": total_resistance,
            "silent_mutations": total_silent,
            "precursor_mutations": total_precursor,
            "novel_mutations": total_novel,
            "nonsense_mutations": total_nonsense,
            "drugs_affected": drugs_affected,
            "precursor_drugs": precursor_drugs,
            "overall_risk": overall_risk,
        },
        "all_mutations": [_mutation_summary(m) for m in all_mutations],
        "warnings": warnings,
    }


def _mutations_from_resistance_genes(
    new_variant: TBVariant,
    ref_variant: TBVariant,
    gene_name: str,
) -> list[Mutation]:
    """Fallback: build mutation list from ResistanceGene objects when
    nucleotide sequences are not available."""
    new_muts: dict[int, Mutation] = {}
    ref_muts: dict[int, Mutation] = {}

    for rg in new_variant.resistance_genes:
        if rg.gene_name == gene_name:
            for m in rg.mutations:
                new_muts[m.position] = m

    for rg in ref_variant.resistance_genes:
        if rg.gene_name == gene_name:
            for m in rg.mutations:
                ref_muts[m.position] = m

    novel_positions = set(new_muts.keys()) - set(ref_muts.keys())
    return [new_muts[p] for p in sorted(novel_positions)]


# ---------------------------------------------------------------------------
# Legacy / backward-compatible functions
# ---------------------------------------------------------------------------

def find_silent_mutations(query_seq: str, reference_seq: str) -> list[dict]:
    """Find synonymous (silent) mutations between two nucleotide sequences."""
    query_seq = query_seq.upper().replace("U", "T")
    reference_seq = reference_seq.upper().replace("U", "T")
    min_len = min(len(query_seq), len(reference_seq))
    silent = []

    for i in range(0, min_len - 2, 3):
        ref_codon = reference_seq[i:i + 3]
        query_codon = query_seq[i:i + 3]

        if ref_codon != query_codon:
            ref_aa = _translate_codon(ref_codon)
            query_aa = _translate_codon(query_codon)

            if ref_aa == query_aa and ref_aa != "X":
                codon_pos = (i // 3) + 1
                silent.append({
                    "codon_position": codon_pos,
                    "ref_codon": ref_codon,
                    "alt_codon": query_codon,
                    "amino_acid": ref_aa,
                    "nucleotide_position": i + 1,
                })

    return silent


def assess_silent_mutation_risk(gene: str, silent_mutations: list[dict]) -> list[dict]:
    """Assess risk of silent mutations using precursor database and proximity."""
    precursors = get_silent_precursors(gene)
    precursor_positions = {p["position"] for p in precursors}

    assessed = []
    for mut in silent_mutations:
        pos = mut["codon_position"]
        risk = "low"
        note = ""

        if pos in precursor_positions:
            matching = next((p for p in precursors if p["position"] == pos), None)
            if matching:
                risk = matching.get("risk", "moderate")
                note = matching.get("note", "Known silent precursor")
        else:
            codon = mut.get("alt_codon", "")
            if codon and len(codon) == 3:
                step = check_one_step_away(gene, pos, codon)
                if step["is_precursor"]:
                    risk = step["risk_level"]
                    note = step["note"]

            if risk == "low":
                db = load_resistance_db()
                gene_data = db["genes"].get(gene, {})
                known_positions = set()
                for k in gene_data.get("mutations", {}).keys():
                    mid = k[1:-1] if len(k) > 2 else ""
                    if mid.lstrip("-").isdigit():
                        known_positions.add(int(mid))
                distance = min((abs(pos - kp) for kp in known_positions), default=999)
                if distance <= 3:
                    risk = "moderate"
                    note = "Within 3 codons of known resistance position"
                elif distance <= 10:
                    risk = "low-moderate"
                    note = "Within 10 codons of known resistance position"

        assessed.append({
            **mut,
            "risk_level": risk,
            "note": note,
            "is_known_precursor": pos in precursor_positions,
        })

    return assessed


def perform_gene_analysis(gene_name: str, sequence: str = "") -> dict:
    """Analyze a gene: GC content, codon usage, known SNPs, precursors."""
    if not sequence:
        gene_id = search_ncbi_gene(gene_name)
        if gene_id:
            sequence = fetch_ncbi_gene_sequence(gene_id) or ""

    gc = _compute_gc_content(sequence)
    codon_usage = _compute_codon_usage(sequence) if sequence else {}

    db = load_resistance_db()
    gene_data = db["genes"].get(gene_name, {})
    ncbi_id = gene_data.get("ncbi_gene_id", "")

    snps = []
    for mut_code, info in gene_data.get("mutations", {}).items():
        if len(mut_code) > 2 and mut_code[0].isalpha() and mut_code[-1].isalpha():
            pos_str = mut_code[1:-1]
            if pos_str.lstrip("-").isdigit():
                snps.append(Mutation(
                    position=int(pos_str),
                    reference_amino_acid=mut_code[0],
                    mutant_amino_acid=mut_code[-1],
                    is_resistance_conferring=True,
                    drug_affected=", ".join(info.get("drugs", [])),
                    resistance_level=ResistanceLevel.HIGH if info.get("resistance") == "high" else ResistanceLevel.MODERATE,
                ))

    precursors = get_silent_precursors(gene_name)
    hotspots = [p["position"] for p in precursors]

    return {
        "gene_name": gene_name,
        "ncbi_id": ncbi_id,
        "sequence": sequence,
        "length": len(sequence),
        "gc_content": gc,
        "codon_usage": codon_usage,
        "snps": [s.model_dump(mode="json") for s in snps],
        "silent_mutation_hotspots": hotspots,
    }


def compare_gene_sequences(query_seq: str, reference_seq: str) -> dict:
    """Align two gene sequences and report SNPs and silent mutations."""
    query_seq = query_seq.upper().replace("U", "T")
    reference_seq = reference_seq.upper().replace("U", "T")

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(query_seq, reference_seq)
    best = alignments[0]

    snps = []
    min_len = min(len(query_seq), len(reference_seq))
    for i in range(min_len):
        if query_seq[i] != reference_seq[i]:
            snps.append({"position": i + 1, "ref": reference_seq[i], "alt": query_seq[i]})

    silent = find_silent_mutations(query_seq, reference_seq)

    return {
        "alignment_score": float(best.score),
        "total_snps": len(snps),
        "snps": snps[:50],
        "silent_mutations": silent,
        "query_length": len(query_seq),
        "reference_length": len(reference_seq),
        "identity_percent": round((min_len - len(snps)) / min_len * 100, 2) if min_len > 0 else 0,
    }


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def _compute_gc_content(sequence: str) -> float:
    if not sequence:
        return 0.0
    seq = sequence.upper()
    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq) * 100, 2)


def _compute_codon_usage(sequence: str) -> dict:
    seq = sequence.upper()
    codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
    counts = Counter(codons)
    total = sum(counts.values())
    if total == 0:
        return {}
    return {codon: round(count / total, 4) for codon, count in sorted(counts.items())}


def translate_codon(codon: str) -> str:
    """Legacy public alias."""
    return _translate_codon(codon)


def identify_mutation_type(ref_codon: str, alt_codon: str) -> str:
    """Return 'silent', 'nonsense', or 'missense'."""
    ref_aa = _translate_codon(ref_codon)
    alt_aa = _translate_codon(alt_codon)
    if ref_aa == alt_aa:
        return "silent"
    if alt_aa == "*":
        return "nonsense"
    return "missense"


def compute_gc_content(sequence: str) -> float:
    """Legacy public alias."""
    return _compute_gc_content(sequence)


def compute_codon_usage(sequence: str) -> dict:
    """Legacy public alias."""
    return _compute_codon_usage(sequence)
