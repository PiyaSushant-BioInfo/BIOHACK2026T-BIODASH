"""Integration test for protein_compare module."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from modules.protein_compare import (
    align_protein_sequences, calculate_sequence_identity,
    compare_drug_binding_sites, calculate_protein_similarity_score,
    find_closest_protein_match, interpret_similarity,
    analyze_mutation_impact, compare_variants_protein_profile,
    align_proteins, compare_protein_sequences,
    fetch_protein_analysis, compute_physicochemical_properties,
)
from schema import TBVariant, ResistanceGene, Mutation, DrugSensitivity, ResistanceLevel
from modules.variant_db import load_known_variants


def test_align_protein_sequences():
    print("--- 1. align_protein_sequences ---")
    ref = "MSDELPHINSIRLGKSFHGRS"
    mut = "MSDELPHINLIRLGKSFHGRS"
    result = align_protein_sequences(ref, mut)
    assert result["global"]["identity_percent"] > 90, f"Expected >90%, got {result['global']['identity_percent']}"
    assert result["global"]["score"] > 0
    assert result["local"]["identity_percent"] > 0
    print(f"  Global: {result['global']['identity_percent']}%, score={result['global']['score']}")
    print(f"  Local:  {result['local']['identity_percent']}%, score={result['local']['score']}")

    empty = align_protein_sequences("", "ABCD")
    assert empty["identity"] == 0.0
    print("  Empty input handled OK")


def test_calculate_sequence_identity():
    print("--- 2. calculate_sequence_identity ---")
    ref = "MSDELPHINSIRLGKSFHGRS"
    mut = "MSDELPHINLIRLGKSFHGRS"
    ident = calculate_sequence_identity(ref, mut)
    assert 0.9 < ident <= 1.0, f"Expected ~0.95, got {ident}"
    print(f"  1-diff identity: {ident}")

    assert calculate_sequence_identity(ref, ref) == 1.0
    assert calculate_sequence_identity("", ref) == 0.0
    print("  Edge cases OK")


def test_compare_drug_binding_sites():
    print("--- 3. compare_drug_binding_sites ---")
    rpob_ref = "X" * 530 + "S" + "X" * 100
    rpob_mut = "X" * 530 + "L" + "X" * 100
    result = compare_drug_binding_sites(rpob_mut, rpob_ref, "rpoB")
    assert result["protein"] == "rpoB"
    assert result["drug"] == "rifampicin"
    assert result["binding_sites_checked"] == 3
    assert result["mutated"] == 1
    print(f"  Checked {result['binding_sites_checked']} sites, {result['mutated']} mutated")
    for sc in result["site_comparisons"]:
        status = "MUTATED" if not sc["identical"] else "ok"
        print(f"    pos {sc['position']}: {sc['known_aa']}->{sc['new_aa']} [{status}] sev={sc['severity']}")
    assert len(result["warnings"]) > 0
    print(f"  Warnings: {result['warnings'][0]}")

    unknown = compare_drug_binding_sites("ABC", "ABC", "unknownGene")
    assert unknown["binding_sites_checked"] == 0
    print("  Unknown gene handled OK")


def test_calculate_protein_similarity_score():
    print("--- 4. calculate_protein_similarity_score ---")
    v1 = TBVariant(
        variant_id="TEST_1",
        resistance_genes=[ResistanceGene(
            gene_name="rpoB",
            mutations=[Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                                is_resistance_conferring=True)],
            drug_target="rifampicin",
        )],
    )
    v2 = TBVariant(
        variant_id="TEST_2",
        resistance_genes=[ResistanceGene(
            gene_name="rpoB",
            mutations=[Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                                is_resistance_conferring=True)],
            drug_target="rifampicin",
        )],
    )
    score_same = calculate_protein_similarity_score(v1, v2)
    assert score_same == 1.0, f"Same mutations should yield 1.0, got {score_same}"
    print(f"  Same mutations (no seqs): {score_same}")

    v3 = TBVariant(
        variant_id="TEST_3",
        resistance_genes=[ResistanceGene(
            gene_name="rpoB",
            mutations=[Mutation(position=526, reference_amino_acid="H", mutant_amino_acid="Y",
                                is_resistance_conferring=True)],
            drug_target="rifampicin",
        )],
    )
    score_diff = calculate_protein_similarity_score(v1, v3)
    assert score_diff < 1.0, f"Different mutations should be < 1.0, got {score_diff}"
    print(f"  Different mutations: {score_diff}")


def test_find_closest_protein_match():
    print("--- 5. find_closest_protein_match ---")
    known = load_known_variants()
    test_var = TBVariant(
        variant_id="TEST_NEW",
        resistance_genes=[
            ResistanceGene(gene_name="rpoB", mutations=[
                Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                         is_resistance_conferring=True),
            ], drug_target="rifampicin"),
            ResistanceGene(gene_name="katG", mutations=[
                Mutation(position=315, reference_amino_acid="S", mutant_amino_acid="T",
                         is_resistance_conferring=True),
            ], drug_target="isoniazid"),
        ],
    )
    closest, cscore = find_closest_protein_match(test_var, known)
    assert closest is not None
    assert cscore > 0
    print(f"  Closest: {closest.variant_id} ({closest.name}), score={cscore}")

    none_var, none_score = find_closest_protein_match(test_var, [])
    assert none_var is None and none_score == 0.0
    print("  Empty list handled OK")


def test_interpret_similarity():
    print("--- 6. interpret_similarity ---")
    high = interpret_similarity(0.95)
    assert high["level"] == "HIGH"
    assert high["requires_full_workup"] is False

    mod = interpret_similarity(0.75)
    assert mod["level"] == "MODERATE"

    low = interpret_similarity(0.40)
    assert low["level"] == "LOW"
    assert low["requires_full_workup"] is True

    for s, expected in [(0.95, "HIGH"), (0.75, "MODERATE"), (0.40, "LOW")]:
        interp = interpret_similarity(s)
        print(f"  {s:.2f} -> {interp['level']}: {interp['recommendation']}")


def test_analyze_mutation_impact():
    print("--- 7. analyze_mutation_impact ---")
    impact = analyze_mutation_impact("rpoB", 531, "S", "L")
    assert impact["predicted_severity"] in ("critical", "high")
    assert impact["at_binding_site"] is True
    assert impact["class_change"] is True
    print(f"  {impact['short_code']}: severity={impact['predicted_severity']}")
    print(f"  class: {impact['ref_class']}->{impact['alt_class']}, binding={impact['at_binding_site']}")
    print(f"  drug={impact['drug_affected']}, size_delta={impact['size_change_daltons']}Da")

    impact2 = analyze_mutation_impact("katG", 315, "S", "T")
    assert impact2["at_binding_site"] is True
    print(f"  katG {impact2['short_code']}: severity={impact2['predicted_severity']}")


def test_legacy_wrappers():
    print("--- 8. Legacy wrappers ---")
    ref = "MSDELPHINSIRLGKSFHGRS"
    mut = "MSDELPHINLIRLGKSFHGRS"
    legacy = align_proteins(ref, mut)
    assert "identity_percent" in legacy
    print(f"  align_proteins: {legacy['identity_percent']}%")

    comp = compare_protein_sequences(ref, mut)
    assert comp["total_differences"] == 1
    for d in comp["differences"]:
        assert "class_change" in d
        print(f"    pos {d['position']}: {d['query_aa']}->{d['reference_aa']} class_change={d['class_change']}")


def test_physicochemical():
    print("--- 9. physicochemical ---")
    props = compute_physicochemical_properties("MSDELPHINSIRLGKSFHGRS")
    assert "molecular_weight" in props
    assert props["molecular_weight"] > 0
    print(f"  MW={props['molecular_weight']}, pI={props['isoelectric_point']}, GRAVY={props['gravy']}")

    empty = compute_physicochemical_properties("")
    assert empty == {}
    print("  Empty handled OK")


def test_compare_variants_protein_profile():
    print("--- 10. compare_variants_protein_profile ---")
    known = load_known_variants()
    test_var = TBVariant(
        variant_id="TEST_PROFILE",
        resistance_genes=[
            ResistanceGene(gene_name="rpoB", mutations=[
                Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                         is_resistance_conferring=True),
            ], drug_target="rifampicin"),
        ],
    )
    profile = compare_variants_protein_profile(test_var, known)
    assert profile["new_variant_id"] == "TEST_PROFILE"
    assert len(profile["comparisons"]) == len(known)
    assert profile["closest_match"]["variant_id"] is not None
    print(f"  Compared against {len(profile['comparisons'])} known variants")
    print(f"  Closest: {profile['closest_match']['variant_id']} (score={profile['closest_match']['score']})")
    print(f"  Mutation impacts: {len(profile['mutation_impacts'])}")
    for mi in profile["mutation_impacts"]:
        print(f"    {mi['short_code']}: {mi['predicted_severity']}")


if __name__ == "__main__":
    test_align_protein_sequences()
    test_calculate_sequence_identity()
    test_compare_drug_binding_sites()
    test_calculate_protein_similarity_score()
    test_find_closest_protein_match()
    test_interpret_similarity()
    test_analyze_mutation_impact()
    test_legacy_wrappers()
    test_physicochemical()
    test_compare_variants_protein_profile()
    print()
    print("ALL TESTS PASSED")
