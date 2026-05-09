"""Integration test for similarity_score module."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from modules.similarity_score import (
    calculate_weighted_score, calculate_risk_score,
    compare_resistance_profiles, full_similarity_analysis,
    interpret_and_recommend,
    compute_mutation_overlap, compute_resistance_risk_score,
    compute_risk_score, predict_resistance_level,
    compare_variant_to_known, find_closest_known_variant,
)
from schema import (
    TBVariant, ResistanceGene, Mutation, DrugSensitivity,
    ResistanceLevel, ComparisonResult, RiskScore,
)
from modules.variant_db import load_known_variants


def test_calculate_weighted_score():
    print("--- 1. calculate_weighted_score ---")
    perfect = calculate_weighted_score(1.0, 1.0, 1.0, 1.0)
    assert perfect == 1.0, f"All 1.0 should give 1.0, got {perfect}"
    print(f"  All 1.0 -> {perfect}")

    zero = calculate_weighted_score(0.0, 0.0, 0.0, 0.0)
    assert zero == 0.0
    print(f"  All 0.0 -> {zero}")

    mixed = calculate_weighted_score(0.9, 0.8, 0.7, 0.6)
    expected = 0.30 * 0.9 + 0.20 * 0.8 + 0.30 * 0.7 + 0.20 * 0.6
    assert abs(mixed - round(expected, 4)) < 0.001
    print(f"  (0.9, 0.8, 0.7, 0.6) -> {mixed} (expected {round(expected, 4)})")

    # Verify weights: protein 30%, gene 20%, binding 30%, resistance 20%
    protein_only = calculate_weighted_score(1.0, 0.0, 0.0, 0.0)
    assert abs(protein_only - 0.30) < 0.001
    gene_only = calculate_weighted_score(0.0, 1.0, 0.0, 0.0)
    assert abs(gene_only - 0.20) < 0.001
    binding_only = calculate_weighted_score(0.0, 0.0, 1.0, 0.0)
    assert abs(binding_only - 0.30) < 0.001
    resist_only = calculate_weighted_score(0.0, 0.0, 0.0, 1.0)
    assert abs(resist_only - 0.20) < 0.001
    print("  Individual weights verified: 30/20/30/20")

    clamped = calculate_weighted_score(1.5, -0.5, 1.0, 1.0)
    assert 0.0 <= clamped <= 1.0
    print(f"  Out-of-range inputs clamped: {clamped}")


def test_calculate_risk_score():
    print("--- 2. calculate_risk_score ---")

    # Susceptible variant
    susceptible = TBVariant(variant_id="SUSC_001", lineage="Lineage 4")
    risk_s = calculate_risk_score(susceptible)
    assert 20 <= risk_s.score <= 35
    print(f"  Susceptible: {risk_s.score} [{risk_s.level.value}] ({risk_s.color.value})")

    # MDR variant
    mdr = TBVariant(
        variant_id="MDR_001",
        lineage="Lineage 4 (Euro-American)",
        resistance_genes=[
            ResistanceGene(gene_name="rpoB", mutations=[
                Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                         is_resistance_conferring=True, drug_affected="rifampicin",
                         resistance_level=ResistanceLevel.HIGH),
            ], drug_target="rifampicin"),
            ResistanceGene(gene_name="katG", mutations=[
                Mutation(position=315, reference_amino_acid="S", mutant_amino_acid="T",
                         is_resistance_conferring=True, drug_affected="isoniazid",
                         resistance_level=ResistanceLevel.HIGH),
            ], drug_target="isoniazid"),
        ],
        drug_resistance={"rifampicin": DrugSensitivity.RESISTANT, "isoniazid": DrugSensitivity.RESISTANT},
    )
    risk_mdr = calculate_risk_score(mdr)
    assert 51 <= risk_mdr.score <= 85
    assert len(risk_mdr.factors) > 0
    print(f"  MDR: {risk_mdr.score} [{risk_mdr.level.value}] ({risk_mdr.color.value})")
    for f in risk_mdr.factors:
        print(f"    - {f}")

    # Beijing lineage MDR (should get bonus)
    beijing_mdr = TBVariant(
        variant_id="BEIJING_MDR",
        lineage="Lineage 2 (Beijing)",
        resistance_genes=mdr.resistance_genes,
        drug_resistance=mdr.drug_resistance,
    )
    risk_beijing = calculate_risk_score(beijing_mdr)
    assert risk_beijing.score > risk_mdr.score
    assert any("Beijing" in f for f in risk_beijing.factors)
    print(f"  Beijing MDR: {risk_beijing.score} [{risk_beijing.level.value}] (bonus applied)")

    # With precursor
    precursor_mdr = TBVariant(
        variant_id="PREC_MDR",
        lineage="Lineage 2 (Beijing)",
        resistance_genes=[
            ResistanceGene(gene_name="rpoB", mutations=[
                Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                         is_resistance_conferring=True, drug_affected="rifampicin",
                         resistance_level=ResistanceLevel.HIGH),
                Mutation(position=514, reference_amino_acid="F", mutant_amino_acid="F",
                         is_synonymous=True, one_step_away_risk=True,
                         one_step_away_drug="rifampicin"),
            ], drug_target="rifampicin"),
            ResistanceGene(gene_name="katG", mutations=[
                Mutation(position=315, reference_amino_acid="S", mutant_amino_acid="T",
                         is_resistance_conferring=True, drug_affected="isoniazid",
                         resistance_level=ResistanceLevel.HIGH),
            ], drug_target="isoniazid"),
        ],
        drug_resistance={"rifampicin": DrugSensitivity.RESISTANT, "isoniazid": DrugSensitivity.RESISTANT},
    )
    risk_prec = calculate_risk_score(precursor_mdr)
    assert risk_prec.score >= risk_beijing.score
    assert any("precursor" in f for f in risk_prec.factors)
    print(f"  Beijing MDR + precursor: {risk_prec.score} [{risk_prec.level.value}]")


def test_compare_resistance_profiles():
    print("--- 3. compare_resistance_profiles ---")
    v1 = TBVariant(
        variant_id="V1",
        drug_resistance={
            "rifampicin": DrugSensitivity.RESISTANT,
            "isoniazid": DrugSensitivity.RESISTANT,
            "ethambutol": DrugSensitivity.SENSITIVE,
        },
    )
    v2 = TBVariant(
        variant_id="V2",
        drug_resistance={
            "rifampicin": DrugSensitivity.RESISTANT,
            "isoniazid": DrugSensitivity.RESISTANT,
            "ethambutol": DrugSensitivity.SENSITIVE,
        },
    )
    score_same = compare_resistance_profiles(v1, v2)
    assert score_same == 1.0
    print(f"  Identical profiles: {score_same}")

    v3 = TBVariant(
        variant_id="V3",
        drug_resistance={
            "rifampicin": DrugSensitivity.SENSITIVE,
            "isoniazid": DrugSensitivity.SENSITIVE,
            "ethambutol": DrugSensitivity.SENSITIVE,
        },
    )
    score_diff = compare_resistance_profiles(v1, v3)
    assert score_diff < 1.0
    print(f"  Different profiles: {score_diff}")

    v_empty = TBVariant(variant_id="EMPTY")
    score_empty = compare_resistance_profiles(v_empty, v_empty)
    assert score_empty == 1.0
    print(f"  Both empty: {score_empty}")


def test_full_similarity_analysis():
    print("--- 4. full_similarity_analysis ---")
    known = load_known_variants()
    new_var = TBVariant(
        variant_id="TEST_NEW",
        lineage="Lineage 2 (Beijing)",
        resistance_genes=[
            ResistanceGene(gene_name="rpoB", mutations=[
                Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                         is_resistance_conferring=True, drug_affected="rifampicin",
                         resistance_level=ResistanceLevel.HIGH),
            ], drug_target="rifampicin"),
            ResistanceGene(gene_name="katG", mutations=[
                Mutation(position=315, reference_amino_acid="S", mutant_amino_acid="T",
                         is_resistance_conferring=True, drug_affected="isoniazid",
                         resistance_level=ResistanceLevel.HIGH),
            ], drug_target="isoniazid"),
        ],
        drug_resistance={"rifampicin": DrugSensitivity.RESISTANT, "isoniazid": DrugSensitivity.RESISTANT},
    )

    results = full_similarity_analysis(new_var, known)
    assert len(results) <= 5
    assert len(results) > 0

    # Verify sorted by weighted_final_score descending
    for i in range(len(results) - 1):
        assert results[i].weighted_final_score >= results[i + 1].weighted_final_score

    print(f"  Got {len(results)} results (top 5)")
    for r in results:
        print(f"    {r.matched_variant_id}: weighted={r.weighted_final_score}%, "
              f"protein={r.protein_similarity_score}%, gene={r.gene_similarity_score}%, "
              f"confidence={r.confidence_level.value}")

    # Best match should be TB_VAR_002 (same MDR profile)
    assert results[0].matched_variant_id == "TB_VAR_002"
    print(f"  Best match: {results[0].matched_variant_id} (expected TB_VAR_002)")


def test_interpret_and_recommend():
    print("--- 5. interpret_and_recommend ---")

    # High similarity
    high = ComparisonResult(
        new_variant_id="NEW",
        matched_variant_id="KNOWN_1",
        weighted_final_score=95.0,
        resistance_mutations=[
            Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                     is_resistance_conferring=True),
        ],
        treatment_recommendation="MDR",
    )
    rec_high = interpret_and_recommend(high)
    assert rec_high["confidence"] == "HIGH"
    assert "protocol_match" in rec_high["flags"]
    print(f"  95% -> {rec_high['confidence']}: {rec_high['flags']}")

    # Moderate similarity
    mod = ComparisonResult(
        new_variant_id="NEW",
        matched_variant_id="KNOWN_2",
        weighted_final_score=80.0,
        resistance_mutations=[
            Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                     is_resistance_conferring=True),
        ],
        novel_mutations=[
            Mutation(position=999, reference_amino_acid="A", mutant_amino_acid="V"),
        ],
        treatment_recommendation="MDR",
    )
    rec_mod = interpret_and_recommend(mod)
    assert rec_mod["confidence"] == "MODERATE"
    assert "monitor_closely" in rec_mod["flags"]
    assert rec_mod["novel_mutations_count"] == 1
    assert len(rec_mod["warnings"]) > 0
    print(f"  80% + novel -> {rec_mod['confidence']}: warnings={rec_mod['warnings']}")

    # Low similarity
    low = ComparisonResult(
        new_variant_id="NEW",
        matched_variant_id="KNOWN_3",
        weighted_final_score=65.0,
        treatment_recommendation="MONO-RESISTANT",
    )
    rec_low = interpret_and_recommend(low)
    assert rec_low["confidence"] == "LOW"
    assert "ast_required" in rec_low["flags"]
    print(f"  65% -> {rec_low['confidence']}: {rec_low['flags']}")

    # Novel variant
    novel = ComparisonResult(
        new_variant_id="NEW",
        matched_variant_id="KNOWN_4",
        weighted_final_score=40.0,
        resistance_mutations=[
            Mutation(position=315, reference_amino_acid="S", mutant_amino_acid="T",
                     is_resistance_conferring=True),
        ],
        treatment_recommendation="SUSCEPTIBLE",
    )
    rec_novel = interpret_and_recommend(novel)
    assert rec_novel["confidence"] == "LOW"
    assert "novel_variant" in rec_novel["flags"]
    assert "full_ast_required" in rec_novel["flags"]
    print(f"  40% -> {rec_novel['confidence']}: {rec_novel['flags']}")

    # With precursor mutations
    with_prec = ComparisonResult(
        new_variant_id="NEW",
        matched_variant_id="KNOWN_5",
        weighted_final_score=85.0,
        silent_mutations=[
            Mutation(position=514, reference_amino_acid="F", mutant_amino_acid="F",
                     is_synonymous=True, one_step_away_risk=True, one_step_away_drug="rifampicin"),
        ],
        treatment_recommendation="MDR",
    )
    rec_prec = interpret_and_recommend(with_prec)
    assert "monitor_precursors" in rec_prec["flags"]
    assert any("precursor" in w for w in rec_prec["warnings"])
    print(f"  85% + precursor -> flags={rec_prec['flags']}")


def test_legacy_compatibility():
    print("--- 6. Legacy compatibility ---")
    known = load_known_variants()
    mdr = None
    for v in known:
        if v.variant_id == "TB_VAR_002":
            mdr = v
            break
    assert mdr is not None

    # compute_resistance_risk_score returns float
    raw = compute_resistance_risk_score(mdr)
    assert isinstance(raw, float)
    assert raw > 0
    print(f"  compute_resistance_risk_score: {raw}")

    # compute_risk_score returns RiskScore
    risk = compute_risk_score(mdr)
    assert isinstance(risk, RiskScore)
    assert risk.score == raw
    print(f"  compute_risk_score: {risk.score} [{risk.level.value}]")

    # predict_resistance_level
    level = predict_resistance_level(risk.score, mdr.all_mutations())
    assert level in ("MDR", "HIGH", "Pre-XDR")
    print(f"  predict_resistance_level: {level}")

    # compare_variant_to_known
    comparisons = compare_variant_to_known(mdr)
    assert isinstance(comparisons, list)
    assert len(comparisons) > 0
    assert all(isinstance(c, ComparisonResult) for c in comparisons)
    print(f"  compare_variant_to_known: {len(comparisons)} results")

    # find_closest_known_variant
    closest = find_closest_known_variant(mdr)
    assert isinstance(closest, ComparisonResult)
    assert closest.matched_variant_id != ""
    print(f"  find_closest_known_variant: {closest.matched_variant_id}")

    # compute_mutation_overlap
    m1 = [Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L")]
    m2 = [Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L"),
          Mutation(position=315, reference_amino_acid="S", mutant_amino_acid="T")]
    overlap = compute_mutation_overlap(m1, m2)
    assert overlap["overlap_count"] == 1
    assert overlap["jaccard_similarity"] == 0.5
    print(f"  compute_mutation_overlap: jaccard={overlap['jaccard_similarity']}")


if __name__ == "__main__":
    test_calculate_weighted_score()
    test_calculate_risk_score()
    test_compare_resistance_profiles()
    test_full_similarity_analysis()
    test_interpret_and_recommend()
    test_legacy_compatibility()
    print()
    print("ALL TESTS PASSED")
