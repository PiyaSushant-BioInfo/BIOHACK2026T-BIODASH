"""Integration test for report_generator module."""
import sys
import os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from modules.report_generator import (
    generate_doctor_report, generate_patient_report,
    generate_pdf_report, generate_comparison_chart_data,
    generate_gene_change_table,
    generate_clinical_report, generate_patient_summary,
)
from schema import (
    ClinicalReport, TBVariant, ResistanceGene, Mutation,
    DrugSensitivity, ResistanceLevel, ComparisonResult,
    RiskScore, RiskLevel, RiskColor, ConfidenceLevel,
)
from datetime import datetime


def _make_mdr_variant():
    return TBVariant(
        variant_id="TEST_RPT_001",
        name="Test MDR Variant",
        lineage="Lineage 2 (Beijing)",
        resistance_genes=[
            ResistanceGene(gene_name="rpoB", mutations=[
                Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                         is_resistance_conferring=True, drug_affected="rifampicin",
                         resistance_level=ResistanceLevel.HIGH),
                Mutation(position=514, reference_amino_acid="F", mutant_amino_acid="F",
                         is_synonymous=True, one_step_away_risk=True,
                         one_step_away_drug="rifampicin",
                         reference_codon="TTC", mutant_codon="TTT"),
            ], drug_target="rifampicin"),
            ResistanceGene(gene_name="katG", mutations=[
                Mutation(position=315, reference_amino_acid="S", mutant_amino_acid="T",
                         is_resistance_conferring=True, drug_affected="isoniazid",
                         resistance_level=ResistanceLevel.HIGH),
            ], drug_target="isoniazid"),
        ],
        drug_resistance={
            "rifampicin": DrugSensitivity.RESISTANT,
            "isoniazid": DrugSensitivity.RESISTANT,
            "ethambutol": DrugSensitivity.SENSITIVE,
        },
    )


def _make_report(variant=None):
    if variant is None:
        variant = _make_mdr_variant()
    return ClinicalReport(
        patient_id="NP-TEST-001",
        variant=variant,
        risk_score=RiskScore(
            variant_id=variant.variant_id,
            score=74.3,
            level=RiskLevel.HIGH,
            color=RiskColor.RED,
            factors=[
                "rpoB S531L: resistance to rifampicin",
                "katG S315T: resistance to isoniazid",
                "Beijing/Lineage 2: higher transmissibility",
                "1 precursor mutation(s) detected",
            ],
        ),
        comparison_result=ComparisonResult(
            new_variant_id=variant.variant_id,
            matched_variant_id="TB_VAR_002",
            weighted_final_score=87.5,
            protein_similarity_score=92.0,
            gene_similarity_score=85.0,
            confidence_level=ConfidenceLevel.HIGH,
            resistance_mutations=[
                Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                         is_resistance_conferring=True),
            ],
            novel_mutations=[
                Mutation(position=999, reference_amino_acid="A", mutant_amino_acid="V"),
            ],
            treatment_recommendation="MDR",
        ),
        treatment="=== Treatment ===\nResistance Class: MDR\nConfidence: 57%\nRegimen: Bdq + Lzd + Mfx + Cfz\n\nRecommended Regimen:\n  - bedaquiline [Available]\n  - linezolid [Available]\n  - moxifloxacin [Available]\n  - clofazimine [Available]\n\nDuration: 26 weeks\n\nMonitoring:\n  * Monthly sputum cultures\n  * ECG weekly first month",
        doctor_report="MDR-TB confirmed. Beijing lineage with precursor mutations detected. Recommend close monitoring.",
        patient_report="Your TB needs special medicines. Follow your doctor's plan carefully.",
    )


OUTPUT_DIR = Path(__file__).resolve().parent / "output" / "test_reports"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def test_generate_doctor_report():
    print("--- 1. generate_doctor_report ---")
    report = _make_report()
    path = generate_doctor_report(report, str(OUTPUT_DIR / "test_doctor.pdf"))
    assert os.path.exists(path), f"PDF not created at {path}"
    size = os.path.getsize(path)
    assert size > 1000, f"PDF too small: {size} bytes"
    print(f"  Doctor report: {path} ({size:,} bytes)")


def test_generate_patient_report():
    print("--- 2. generate_patient_report ---")
    report = _make_report()
    path = generate_patient_report(report, str(OUTPUT_DIR / "test_patient.pdf"))
    assert os.path.exists(path), f"PDF not created at {path}"
    size = os.path.getsize(path)
    assert size > 1000, f"PDF too small: {size} bytes"
    print(f"  Patient report: {path} ({size:,} bytes)")


def test_generate_pdf_report_dispatch():
    print("--- 3. generate_pdf_report dispatch ---")
    report = _make_report()

    doc_path = generate_pdf_report(report, "doctor", str(OUTPUT_DIR / "test_dispatch_doctor.pdf"))
    assert os.path.exists(doc_path)
    print(f"  doctor dispatch: OK")

    pat_path = generate_pdf_report(report, "patient", str(OUTPUT_DIR / "test_dispatch_patient.pdf"))
    assert os.path.exists(pat_path)
    print(f"  patient dispatch: OK")


def test_susceptible_variant_report():
    print("--- 4. susceptible variant report ---")
    variant = TBVariant(
        variant_id="TEST_SUSC",
        name="Susceptible Strain",
        lineage="Lineage 4",
    )
    report = ClinicalReport(
        patient_id="NP-TEST-SUSC",
        variant=variant,
        risk_score=RiskScore(
            variant_id="TEST_SUSC", score=25.0,
            level=RiskLevel.LOW, color=RiskColor.GREEN, factors=[],
        ),
    )
    doc_path = generate_doctor_report(report, str(OUTPUT_DIR / "test_susc_doctor.pdf"))
    pat_path = generate_patient_report(report, str(OUTPUT_DIR / "test_susc_patient.pdf"))
    assert os.path.exists(doc_path)
    assert os.path.exists(pat_path)
    print(f"  Susceptible doctor: {os.path.getsize(doc_path):,} bytes")
    print(f"  Susceptible patient: {os.path.getsize(pat_path):,} bytes")


def test_xdr_variant_report():
    print("--- 5. XDR variant report ---")
    variant = TBVariant(
        variant_id="TEST_XDR",
        name="XDR Strain",
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
            "bedaquiline": DrugSensitivity.RESISTANT,
        },
    )
    report = ClinicalReport(
        patient_id="NP-TEST-XDR",
        variant=variant,
        risk_score=RiskScore(
            variant_id="TEST_XDR", score=95.0,
            level=RiskLevel.CRITICAL, color=RiskColor.BLACK,
            factors=["XDR profile", "Beijing lineage"],
        ),
    )
    path = generate_patient_report(report, str(OUTPUT_DIR / "test_xdr_patient.pdf"))
    assert os.path.exists(path)
    print(f"  XDR patient report: {os.path.getsize(path):,} bytes")


def test_generate_comparison_chart_data():
    print("--- 6. generate_comparison_chart_data ---")

    results = [
        ComparisonResult(
            new_variant_id="NEW", matched_variant_id="TB_VAR_002",
            weighted_final_score=87.5, protein_similarity_score=92.0,
            gene_similarity_score=85.0, confidence_level=ConfidenceLevel.HIGH,
            resistance_mutations=[
                Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                         is_resistance_conferring=True),
            ],
            novel_mutations=[],
        ),
        ComparisonResult(
            new_variant_id="NEW", matched_variant_id="TB_VAR_003",
            weighted_final_score=72.0, protein_similarity_score=78.0,
            gene_similarity_score=70.0, confidence_level=ConfidenceLevel.MODERATE,
            resistance_mutations=[],
            novel_mutations=[
                Mutation(position=999, reference_amino_acid="A", mutant_amino_acid="V"),
            ],
        ),
    ]

    data = generate_comparison_chart_data(results)
    assert data["chart_ready"] is True
    assert len(data["labels"]) == 2
    assert data["labels"][0] == "TB_VAR_002"
    assert data["weighted_scores"][0] == 87.5
    assert data["protein_scores"][0] == 92.0
    assert data["gene_scores"][0] == 85.0
    assert data["resistance_mutation_counts"][0] == 1
    assert data["novel_mutation_counts"][1] == 1
    assert data["confidence_levels"][0] == "HIGH"
    assert data["top_match_id"] == "TB_VAR_002"
    assert data["top_match_score"] == 87.5
    print(f"  Labels: {data['labels']}")
    print(f"  Scores: {data['weighted_scores']}")
    print(f"  Top match: {data['top_match_id']} at {data['top_match_score']}%")

    empty = generate_comparison_chart_data([])
    assert empty["chart_ready"] is False
    assert len(empty["labels"]) == 0
    print("  Empty data: chart_ready=False")


def test_generate_gene_change_table():
    print("--- 7. generate_gene_change_table ---")

    mutations = [
        Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                 is_resistance_conferring=True, drug_affected="rifampicin",
                 resistance_level=ResistanceLevel.HIGH),
        Mutation(position=514, reference_amino_acid="F", mutant_amino_acid="F",
                 is_synonymous=True, one_step_away_risk=True,
                 one_step_away_drug="rifampicin"),
        Mutation(position=315, reference_amino_acid="S", mutant_amino_acid="T",
                 is_resistance_conferring=True, drug_affected="isoniazid",
                 resistance_level=ResistanceLevel.HIGH),
        Mutation(position=999, reference_amino_acid="A", mutant_amino_acid="V"),
    ]

    gene_names = {531: "rpoB", 514: "rpoB", 315: "katG"}

    table = generate_gene_change_table(mutations, gene_names)
    assert "Gene" in table
    assert "Pos" in table
    assert "Change" in table
    assert "Drug Affected" in table
    assert "Risk" in table
    assert "Type" in table
    assert "531" in table
    assert "HIGH" in table
    assert "PRECURSOR" in table
    assert "unknown" in table

    lines = table.split("\n")
    assert len(lines) == 6  # header + sep + 4 mutations
    print(f"  Table ({len(lines)} lines):")
    for line in lines:
        print(f"    {line}")

    empty_table = generate_gene_change_table([])
    assert empty_table == "No mutations detected."
    print("  Empty: OK")


def test_backward_compatible_aliases():
    print("--- 8. backward-compatible aliases ---")
    report = _make_report()

    path1 = generate_clinical_report(report, str(OUTPUT_DIR / "test_alias_clinical.pdf"))
    assert os.path.exists(path1)
    print(f"  generate_clinical_report: OK")

    path2 = generate_patient_summary(report, str(OUTPUT_DIR / "test_alias_patient.pdf"))
    assert os.path.exists(path2)
    print(f"  generate_patient_summary: OK")


if __name__ == "__main__":
    test_generate_doctor_report()
    test_generate_patient_report()
    test_generate_pdf_report_dispatch()
    test_susceptible_variant_report()
    test_xdr_variant_report()
    test_generate_comparison_chart_data()
    test_generate_gene_change_table()
    test_backward_compatible_aliases()
    print()
    print("ALL TESTS PASSED")
