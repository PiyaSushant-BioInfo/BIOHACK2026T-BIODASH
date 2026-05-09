"""Integration test for treatment module."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from modules.treatment import (
    get_treatment_for_known_variant, infer_treatment_from_similarity,
    generate_treatment_recommendation, check_drug_interactions,
    format_regimen_string, format_treatment_summary,
    get_treatment_centers, determine_regimen_key,
)
from schema import (
    TBVariant, ResistanceGene, Mutation, DrugSensitivity,
    ResistanceLevel, ComparisonResult,
)
from modules.variant_db import load_known_variants


def test_get_treatment_for_known_variant():
    print("--- 1. get_treatment_for_known_variant ---")

    # MDR variant
    result = get_treatment_for_known_variant("TB_VAR_002")
    assert "error" not in result
    assert result["resistance_class"] == "MDR"
    assert "isoniazid" in result["contraindicated_drugs"]
    assert "rifampicin" in result["contraindicated_drugs"]
    assert len(result["recommended_regimen"]) > 0
    assert result["duration_weeks"] > 0
    assert len(result["monitoring_notes"]) > 0
    assert result["regimen_string"] != ""

    # No isoniazid or rifampicin in recommended
    for drug in result["recommended_regimen"]:
        assert drug.lower() not in ("isoniazid", "rifampicin"), f"Resistant drug {drug} in regimen!"
    print(f"  TB_VAR_002 (MDR): {result['resistance_class']}")
    print(f"    Regimen: {result['regimen_string']}")
    print(f"    Drugs: {result['recommended_regimen']}")
    print(f"    Contraindicated: {result['contraindicated_drugs']}")
    print(f"    Duration: {result['duration_weeks']} weeks")
    print(f"    NTP compliant: {result['nepal_ntp_compliant']}")
    if result["drug_interactions"]:
        print(f"    Interactions: {result['drug_interactions']}")

    # Susceptible variant
    result_s = get_treatment_for_known_variant("TB_VAR_001")
    assert "error" not in result_s
    # TB_VAR_001 has isoniazid resistance
    print(f"  TB_VAR_001: class={result_s['resistance_class']}, regimen={result_s['regimen_string']}")

    # Unknown variant
    result_u = get_treatment_for_known_variant("NONEXISTENT")
    assert "error" in result_u
    print("  Unknown variant: error handled OK")


def test_infer_treatment_from_similarity():
    print("--- 2. infer_treatment_from_similarity ---")

    # High similarity to MDR variant
    comparison = ComparisonResult(
        new_variant_id="NEW_001",
        matched_variant_id="TB_VAR_002",
        weighted_final_score=85.0,
        resistance_mutations=[
            Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                     is_resistance_conferring=True, drug_affected="rifampicin"),
        ],
        treatment_recommendation="MDR",
    )
    new_var = TBVariant(
        variant_id="NEW_001",
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
    result = infer_treatment_from_similarity(comparison, new_var)
    assert result["confidence"] in ("HIGH", "MODERATE")
    assert result["based_on_variant"] == "TB_VAR_002"
    for drug in result["recommended_regimen"]:
        assert drug.lower() not in ("isoniazid", "rifampicin")
    print(f"  High similarity inference:")
    print(f"    Based on: {result['based_on_variant']}")
    print(f"    Confidence: {result['confidence']}")
    print(f"    Regimen: {result['regimen_string']}")
    print(f"    Contraindicated: {result['contraindicated_drugs']}")

    # With additional resistance
    new_var_extra = TBVariant(
        variant_id="NEW_002",
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
            ResistanceGene(gene_name="gyrA", mutations=[
                Mutation(position=94, reference_amino_acid="D", mutant_amino_acid="G",
                         is_resistance_conferring=True, drug_affected="levofloxacin, moxifloxacin",
                         resistance_level=ResistanceLevel.HIGH),
            ], drug_target="fluoroquinolones"),
        ],
        drug_resistance={
            "rifampicin": DrugSensitivity.RESISTANT,
            "isoniazid": DrugSensitivity.RESISTANT,
            "levofloxacin": DrugSensitivity.RESISTANT,
            "moxifloxacin": DrugSensitivity.RESISTANT,
        },
    )
    result_extra = infer_treatment_from_similarity(comparison, new_var_extra)
    assert len(result_extra["warnings"]) > 0
    for drug in result_extra["recommended_regimen"]:
        assert drug.lower() not in ("isoniazid", "rifampicin", "levofloxacin", "moxifloxacin")
    print(f"  With additional FQ resistance:")
    print(f"    Contraindicated: {result_extra['contraindicated_drugs']}")
    print(f"    Regimen: {result_extra['regimen_string']}")
    if result_extra.get("substitutions"):
        for s in result_extra["substitutions"]:
            print(f"    Substitution: {s['removed']} -> {s['substituted']}")
    for w in result_extra["warnings"]:
        print(f"    Warning: {w}")


def test_generate_treatment_recommendation():
    print("--- 3. generate_treatment_recommendation ---")

    # Legacy call with string resistance_class
    mdr_var = TBVariant(
        variant_id="TEST_MDR",
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

    rec = generate_treatment_recommendation(mdr_var, "MDR")
    assert rec["resistance_class"] == "MDR"
    assert len(rec["recommended_regimen"]) > 0
    assert rec["duration_weeks"] > 0
    assert rec["duration_months"] > 0
    assert isinstance(rec["warnings"], list)
    assert isinstance(rec["next_steps"], list)
    assert rec["nepal_ntp_compliant"] is not None
    assert rec["regimen_string"] != ""
    assert "drugs_to_avoid" in rec
    assert "drugs_to_monitor" in rec
    assert "confidence" in rec

    print(f"  MDR recommendation:")
    print(f"    Regimen: {rec['regimen_string']}")
    print(f"    Duration: {rec['duration_months']} months ({rec['duration_weeks']} weeks)")
    print(f"    Confidence: {rec['confidence']} ({rec['confidence_score']*100:.0f}%)")
    print(f"    NTP compliant: {rec['nepal_ntp_compliant']}")
    print(f"    Drugs to avoid: {rec['drugs_to_avoid']}")
    print(f"    Drugs to monitor: {rec['drugs_to_monitor']}")
    print(f"    Next steps: {len(rec['next_steps'])}")

    # With precursor mutations — precursor for a drug NOT already resistant
    prec_var = TBVariant(
        variant_id="TEST_PREC",
        resistance_genes=[
            ResistanceGene(gene_name="katG", mutations=[
                Mutation(position=315, reference_amino_acid="S", mutant_amino_acid="T",
                         is_resistance_conferring=True, drug_affected="isoniazid",
                         resistance_level=ResistanceLevel.HIGH),
            ], drug_target="isoniazid"),
            ResistanceGene(gene_name="rpoB", mutations=[
                Mutation(position=514, reference_amino_acid="F", mutant_amino_acid="F",
                         is_synonymous=True, one_step_away_risk=True,
                         one_step_away_drug="rifampicin"),
            ], drug_target="rifampicin"),
        ],
        drug_resistance={"isoniazid": DrugSensitivity.RESISTANT},
    )
    rec_prec = generate_treatment_recommendation(prec_var, "MONO-RESISTANT")
    assert len(rec_prec["drugs_to_monitor"]) > 0
    assert "rifampicin" in rec_prec["drugs_to_monitor"]
    assert any("precursor" in w.lower() or "monitor" in w.lower() for w in rec_prec["warnings"])
    print(f"  With precursor: monitor={rec_prec['drugs_to_monitor']}")


def test_check_drug_interactions():
    print("--- 4. check_drug_interactions ---")

    # QT prolongation
    warnings = check_drug_interactions(["bedaquiline", "moxifloxacin", "linezolid"])
    assert len(warnings) > 0
    assert any("QT" in w for w in warnings)
    print(f"  Bdq + Mfx + Lzd: {len(warnings)} interaction(s)")
    for w in warnings:
        print(f"    {w}")

    # Major interaction
    warnings_rif_bdq = check_drug_interactions(["rifampicin", "bedaquiline"])
    assert len(warnings_rif_bdq) > 0
    assert any("HIGH" in w for w in warnings_rif_bdq)
    print(f"  Rif + Bdq: {warnings_rif_bdq[0]}")

    # No interactions
    warnings_safe = check_drug_interactions(["isoniazid", "rifampicin", "pyrazinamide", "ethambutol"])
    assert len(warnings_safe) == 0
    print("  HRZE: no interactions")

    # Empty
    assert check_drug_interactions([]) == []
    print("  Empty: OK")


def test_format_regimen_string():
    print("--- 5. format_regimen_string ---")

    # Standard first-line
    s1 = format_regimen_string(["isoniazid", "rifampicin", "pyrazinamide", "ethambutol"])
    assert s1 == "2HRZE/4HR"
    print(f"  First-line: {s1}")

    # Second-line
    s2 = format_regimen_string(["bedaquiline", "pretomanid", "linezolid", "moxifloxacin"])
    assert "Bdq" in s2 and "Lzd" in s2
    print(f"  BPaL-M: {s2}")

    # Mixed
    s3 = format_regimen_string(["isoniazid", "rifampicin", "pyrazinamide", "ethambutol", "levofloxacin"])
    assert "Lfx" in s3
    print(f"  HRZE+Lfx: {s3}")

    # Empty
    assert format_regimen_string([]) == ""
    print("  Empty: OK")

    # Partial first-line
    s4 = format_regimen_string(["rifampicin", "ethambutol", "pyrazinamide", "levofloxacin"])
    print(f"  REZ+Lfx: {s4}")


def test_format_treatment_summary():
    print("--- 6. format_treatment_summary ---")
    rec = generate_treatment_recommendation(
        TBVariant(
            variant_id="FMT_TEST",
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
        ),
        "MDR",
    )
    summary = format_treatment_summary(rec)
    assert "Treatment Recommendation" in summary
    assert "Resistance Class: MDR" in summary
    assert "Recommended Regimen:" in summary
    assert "Monitoring:" in summary
    assert "Next Steps:" in summary
    assert "Nepal NTP Compliant:" in summary
    print("  Format output contains all sections")
    # Print just first 15 lines as sample
    for line in summary.split("\n")[:15]:
        print(f"    {line}")
    print("    ...")


def test_treatment_centers():
    print("--- 7. get_treatment_centers ---")
    mdr_centers = get_treatment_centers("MDR")
    assert len(mdr_centers) > 0
    print(f"  MDR centers: {len(mdr_centers)}")
    for c in mdr_centers:
        print(f"    - {c['name']} ({c['location']})")

    xdr_centers = get_treatment_centers("XDR")
    assert len(xdr_centers) > 0
    assert len(xdr_centers) <= len(mdr_centers)
    print(f"  XDR centers: {len(xdr_centers)}")

    all_centers = get_treatment_centers("SUSCEPTIBLE")
    assert len(all_centers) >= len(mdr_centers)
    print(f"  Susceptible (all): {len(all_centers)}")


def test_legacy_determine_regimen_key():
    print("--- 8. Legacy determine_regimen_key ---")
    assert determine_regimen_key("MDR", []) == "mdr"
    assert determine_regimen_key("XDR", []) == "xdr"
    assert determine_regimen_key("Pre-XDR", []) == "pre_xdr"
    assert determine_regimen_key("LOW", ["isoniazid"]) == "isoniazid_resistant"
    assert determine_regimen_key("LOW", []) == "susceptible"
    assert determine_regimen_key("LOW", ["rifampicin", "isoniazid"]) == "mdr"
    print("  All regimen key mappings OK")


if __name__ == "__main__":
    test_get_treatment_for_known_variant()
    test_infer_treatment_from_similarity()
    test_generate_treatment_recommendation()
    test_check_drug_interactions()
    test_format_regimen_string()
    test_format_treatment_summary()
    test_treatment_centers()
    test_legacy_determine_regimen_key()
    print()
    print("ALL TESTS PASSED")
