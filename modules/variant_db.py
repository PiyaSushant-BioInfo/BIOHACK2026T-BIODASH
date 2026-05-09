import json
from pathlib import Path
from typing import Optional
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from schema import (
    TBVariant, Mutation, ResistanceGene, ResistanceLevel,
    DrugSensitivity, RiskScore,
)

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def _map_resistance_level(raw: str) -> ResistanceLevel:
    mapping = {
        "high": ResistanceLevel.HIGH,
        "moderate": ResistanceLevel.MODERATE,
        "low": ResistanceLevel.LOW,
        "critical": ResistanceLevel.CRITICAL,
    }
    return mapping.get(raw, ResistanceLevel.LOW)


def _build_drug_resistance(resistant: list[str], susceptible: list[str]) -> dict[str, DrugSensitivity]:
    result = {}
    for d in resistant:
        result[d] = DrugSensitivity.RESISTANT
    for d in susceptible:
        result[d] = DrugSensitivity.SENSITIVE
    return result


def _parse_variant_mutations(raw_mutations: list[dict]) -> dict[str, list[Mutation]]:
    """Group raw JSON mutations by gene, returning {gene_name: [Mutation, ...]}."""
    by_gene: dict[str, list[Mutation]] = {}
    for m in raw_mutations:
        gene = m["gene"]
        ref_aa = m.get("ref_amino_acid") or ""
        alt_aa = m.get("alt_amino_acid") or ""
        is_syn = (ref_aa == alt_aa) and ref_aa != ""
        mut_type_raw = m.get("mutation_type", "missense")
        is_res = mut_type_raw not in ("silent", "promoter") and not is_syn

        mutation = Mutation(
            position=m["position"],
            reference_codon=m.get("ref_nucleotide", ""),
            mutant_codon=m.get("alt_nucleotide", ""),
            reference_amino_acid=ref_aa,
            mutant_amino_acid=alt_aa,
            is_synonymous=is_syn,
            is_resistance_conferring=is_res,
            drug_affected=m.get("clinical_significance", "").split(" resistance")[0] if "resistance" in m.get("clinical_significance", "") else "",
            resistance_level=_map_resistance_level(m.get("resistance_level_str", "moderate")),
        )
        by_gene.setdefault(gene, []).append(mutation)
    return by_gene


def load_known_variants() -> list[TBVariant]:
    path = DATA_DIR / "known_variants.json"
    with open(path, "r") as f:
        data = json.load(f)

    db = load_resistance_db()
    variants = []
    for v in data["variants"]:
        grouped = _parse_variant_mutations(v["mutations"])
        resistance_genes = []
        for gene_name, mutations in grouped.items():
            gene_info = db["genes"].get(gene_name, {})
            resistance_genes.append(ResistanceGene(
                gene_name=gene_name,
                mutations=mutations,
                drug_target=gene_info.get("drug_target", ""),
            ))

        drug_resistance = _build_drug_resistance(
            v.get("resistant_drugs", []),
            v.get("susceptible_drugs", []),
        )

        variants.append(TBVariant(
            variant_id=v["variant_id"],
            name=v["strain_name"],
            lineage=v["lineage"],
            drug_resistance=drug_resistance,
            resistance_genes=resistance_genes,
            source=v.get("source", ""),
        ))
    return variants


def load_resistance_db() -> dict:
    path = DATA_DIR / "resistance_mutations.json"
    with open(path, "r") as f:
        return json.load(f)


def get_variant_by_id(variant_id: str) -> Optional[TBVariant]:
    for v in load_known_variants():
        if v.variant_id == variant_id:
            return v
    return None


def search_variants_by_gene(gene: str) -> list[TBVariant]:
    return [
        v for v in load_known_variants()
        if any(rg.gene_name == gene for rg in v.resistance_genes)
    ]


def search_variants_by_drug(drug: str) -> list[TBVariant]:
    return [
        v for v in load_known_variants()
        if drug.lower() in [d.lower() for d in v.resistant_drugs()]
    ]


def search_variants_by_resistance_drug(drug: str, status: DrugSensitivity = DrugSensitivity.RESISTANT) -> list[TBVariant]:
    return [
        v for v in load_known_variants()
        if v.drug_resistance.get(drug) == status
    ]


def add_variant(variant: TBVariant) -> None:
    path = DATA_DIR / "known_variants.json"
    with open(path, "r") as f:
        data = json.load(f)

    mutations_list = []
    for rg in variant.resistance_genes:
        for m in rg.mutations:
            mutations_list.append({
                "gene": rg.gene_name,
                "position": m.position,
                "ref_nucleotide": m.reference_codon,
                "alt_nucleotide": m.mutant_codon,
                "ref_amino_acid": m.reference_amino_acid,
                "alt_amino_acid": m.mutant_amino_acid,
                "codon_position": m.position,
                "mutation_type": "silent" if m.is_synonymous else "missense",
                "frequency": 0.0,
                "clinical_significance": f"{m.drug_affected} resistance" if m.is_resistance_conferring else "",
            })

    variant_dict = {
        "variant_id": variant.variant_id,
        "strain_name": variant.name,
        "lineage": variant.lineage,
        "mutations": mutations_list,
        "resistance_level": "multi_drug_resistant" if variant.has_mdr_profile() else "susceptible",
        "resistant_drugs": variant.resistant_drugs(),
        "susceptible_drugs": variant.susceptible_drugs(),
        "source": variant.source,
        "date_identified": variant.last_updated.isoformat(),
        "geographic_origin": "",
    }

    data["variants"].append(variant_dict)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def get_mutation_info(gene: str, mutation_code: str) -> Optional[dict]:
    db = load_resistance_db()
    gene_data = db["genes"].get(gene)
    if not gene_data:
        return None
    mut_info = gene_data.get("mutations", {}).get(mutation_code)
    if mut_info:
        return {**mut_info, "gene": gene, "mutation": mutation_code}
    return None


def get_silent_precursors(gene: str) -> list[dict]:
    db = load_resistance_db()
    gene_data = db["genes"].get(gene, {})
    return [
        {"mutation": k, **v}
        for k, v in gene_data.get("silent_precursors", {}).items()
    ]
