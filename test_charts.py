"""Integration test for charts module."""
import sys
import os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from modules.charts import (
    similarity_radar_chart, resistance_profile_heatmap,
    mutation_map, risk_gauge, similarity_bar_chart,
    lineage_distribution_nepal, save_chart, generate_all_charts,
    resistance_profile_chart, mutation_frequency_chart,
    comparison_heatmap, treatment_availability_chart,
    risk_score_gauge,
)
from schema import (
    TBVariant, ResistanceGene, Mutation, ComparisonResult,
    RiskScore, RiskLevel, RiskColor, ClinicalReport,
    DrugSensitivity, ResistanceLevel, ConfidenceLevel,
)
import plotly.graph_objects as go


OUTPUT_DIR = Path(__file__).resolve().parent / "output" / "test_charts"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def _make_mdr_variant():
    return TBVariant(
        variant_id="CHART_TEST_001",
        name="Test MDR Variant",
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
        drug_resistance={
            "rifampicin": DrugSensitivity.RESISTANT,
            "isoniazid": DrugSensitivity.RESISTANT,
            "ethambutol": DrugSensitivity.SENSITIVE,
        },
    )


def _make_comparison():
    return ComparisonResult(
        new_variant_id="CHART_TEST_001",
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
    )


def _make_comparisons():
    return [
        _make_comparison(),
        ComparisonResult(
            new_variant_id="CHART_TEST_001",
            matched_variant_id="TB_VAR_003",
            weighted_final_score=72.0,
            protein_similarity_score=78.0,
            gene_similarity_score=70.0,
            confidence_level=ConfidenceLevel.MODERATE,
            resistance_mutations=[],
            treatment_recommendation="Pre-XDR",
        ),
        ComparisonResult(
            new_variant_id="CHART_TEST_001",
            matched_variant_id="TB_VAR_001",
            weighted_final_score=55.0,
            protein_similarity_score=60.0,
            gene_similarity_score=50.0,
            confidence_level=ConfidenceLevel.LOW,
            resistance_mutations=[],
            treatment_recommendation="MONO-RESISTANT",
        ),
    ]


def test_similarity_radar_chart():
    print("--- 1. similarity_radar_chart ---")
    comp = _make_comparison()
    fig = similarity_radar_chart(comp)
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 2
    assert fig.data[0].type == "scatterpolar"
    print(f"  Traces: {len(fig.data)}")
    print(f"  Values: {fig.data[0].r}")

    fig2 = similarity_radar_chart(comp, protein_score=95, gene_score=80,
                                   binding_score=88, resistance_score=90)
    assert fig2.data[0].r[0] == 95
    print("  Custom scores: OK")


def test_resistance_profile_heatmap():
    print("--- 2. resistance_profile_heatmap ---")
    v1 = _make_mdr_variant()
    v2 = TBVariant(
        variant_id="SUSC_001", name="Susceptible",
        drug_resistance={
            "rifampicin": DrugSensitivity.SENSITIVE,
            "isoniazid": DrugSensitivity.SENSITIVE,
            "ethambutol": DrugSensitivity.SENSITIVE,
        },
    )
    fig = resistance_profile_heatmap([v1, v2])
    assert isinstance(fig, go.Figure)
    assert fig.data[0].type == "heatmap"
    assert len(fig.data[0].y) == 2
    print(f"  Variants: {fig.data[0].y}")
    print(f"  Drugs: {list(fig.data[0].x)}")


def test_mutation_map():
    print("--- 3. mutation_map ---")
    mutations = [
        Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                 is_resistance_conferring=True, drug_affected="rifampicin",
                 resistance_level=ResistanceLevel.HIGH),
        Mutation(position=514, reference_amino_acid="F", mutant_amino_acid="F",
                 is_synonymous=True, one_step_away_risk=True),
        Mutation(position=526, reference_amino_acid="H", mutant_amino_acid="D",
                 is_resistance_conferring=True, drug_affected="rifampicin",
                 resistance_level=ResistanceLevel.HIGH),
        Mutation(position=450, reference_amino_acid="A", mutant_amino_acid="V"),
        Mutation(position=300, reference_amino_acid="G", mutant_amino_acid="G",
                 is_synonymous=True),
    ]
    fig = mutation_map(mutations, "rpoB", 1172)
    assert isinstance(fig, go.Figure)
    trace_names = [t.name for t in fig.data if t.name]
    assert "Resistance" in trace_names
    assert "Silent Precursor" in trace_names
    print(f"  Traces: {trace_names}")
    print(f"  Gene length: 1172 codons")


def test_risk_gauge():
    print("--- 4. risk_gauge ---")
    risk = RiskScore(
        variant_id="TEST", score=74.3,
        level=RiskLevel.HIGH, color=RiskColor.RED,
        factors=["test factor"],
    )
    fig = risk_gauge(risk)
    assert isinstance(fig, go.Figure)
    assert fig.data[0].type == "indicator"
    assert fig.data[0].value == 74.3
    print(f"  Value: {fig.data[0].value}")

    risk_low = RiskScore(
        variant_id="TEST_LOW", score=15.0,
        level=RiskLevel.LOW, color=RiskColor.GREEN, factors=[],
    )
    fig_low = risk_gauge(risk_low)
    assert fig_low.data[0].value == 15.0
    print(f"  Low risk: {fig_low.data[0].value}")

    risk_crit = RiskScore(
        variant_id="TEST_CRIT", score=95.0,
        level=RiskLevel.CRITICAL, color=RiskColor.BLACK, factors=[],
    )
    fig_crit = risk_gauge(risk_crit)
    assert fig_crit.data[0].value == 95.0
    print(f"  Critical risk: {fig_crit.data[0].value}")


def test_similarity_bar_chart():
    print("--- 5. similarity_bar_chart ---")
    comparisons = _make_comparisons()
    fig = similarity_bar_chart(comparisons)
    assert isinstance(fig, go.Figure)
    assert fig.data[0].orientation == "h"
    assert len(fig.data[0].y) == 3
    print(f"  Variants: {list(fig.data[0].y)}")
    print(f"  Scores: {list(fig.data[0].x)}")

    single = similarity_bar_chart([_make_comparison()])
    assert len(single.data[0].y) == 1
    print("  Single variant: OK")


def test_lineage_distribution_nepal():
    print("--- 6. lineage_distribution_nepal ---")
    variants = [_make_mdr_variant()]
    fig = lineage_distribution_nepal(variants)
    assert isinstance(fig, go.Figure)
    assert fig.data[0].type == "pie"
    assert len(fig.data[0].labels) == 4
    assert sum(fig.data[0].values) > 99
    print(f"  Labels: {list(fig.data[0].labels)}")
    print(f"  Values: {list(fig.data[0].values)}")

    pull = list(fig.data[0].pull)
    assert pull[0] == 0.1, "Beijing should be highlighted"
    print(f"  Beijing highlighted: pull={pull}")

    fig_l4 = lineage_distribution_nepal([], highlight_lineage="Lineage 4 (Euro-American)")
    pull_l4 = list(fig_l4.data[0].pull)
    assert pull_l4[2] == 0.1, "Euro-American should be highlighted"
    print(f"  Euro-American highlighted: pull={pull_l4}")

    fig_none = lineage_distribution_nepal([])
    pull_none = list(fig_none.data[0].pull)
    assert all(p == 0 for p in pull_none)
    print("  No highlight: OK")


def test_save_chart():
    print("--- 7. save_chart ---")
    fig = risk_gauge(RiskScore(
        variant_id="SAVE_TEST", score=50.0,
        level=RiskLevel.MODERATE, color=RiskColor.YELLOW, factors=[],
    ))

    html_path = save_chart(fig, str(OUTPUT_DIR / "test_save.html"), "html")
    assert os.path.exists(html_path)
    assert html_path.endswith(".html")
    size = os.path.getsize(html_path)
    assert size > 100
    print(f"  HTML: {html_path} ({size:,} bytes)")

    html_path2 = save_chart(fig, str(OUTPUT_DIR / "test_save_noext"), "html")
    assert html_path2.endswith(".html")
    print(f"  HTML auto-ext: {html_path2}")


def test_generate_all_charts():
    print("--- 8. generate_all_charts ---")
    variant = _make_mdr_variant()
    report = ClinicalReport(
        patient_id="NP-CHART-TEST",
        variant=variant,
        risk_score=RiskScore(
            variant_id=variant.variant_id, score=74.3,
            level=RiskLevel.HIGH, color=RiskColor.RED,
            factors=["test"],
        ),
        comparison_result=_make_comparison(),
    )
    comparisons = _make_comparisons()

    charts = generate_all_charts(report, comparisons, str(OUTPUT_DIR / "all"))
    assert isinstance(charts, dict)
    assert len(charts) > 0
    print(f"  Generated {len(charts)} charts:")
    for name, path in charts.items():
        exists = os.path.exists(path)
        size = os.path.getsize(path) if exists else 0
        print(f"    {name}: {path} ({size:,} bytes) {'OK' if exists else 'MISSING'}")
        assert exists, f"Chart {name} not created"

    assert "similarity_radar" in charts
    assert "similarity_bar" in charts
    assert "risk_gauge" in charts
    assert "lineage_distribution" in charts
    print("  All expected charts present")


def test_legacy_matplotlib_charts():
    print("--- 9. legacy matplotlib charts ---")
    variant = _make_mdr_variant()

    path = resistance_profile_chart(variant, str(OUTPUT_DIR / "legacy_resistance.png"))
    assert os.path.exists(path)
    print(f"  resistance_profile_chart: {os.path.getsize(path):,} bytes")

    path = mutation_frequency_chart(variant, str(OUTPUT_DIR / "legacy_mutation.png"))
    assert os.path.exists(path)
    print(f"  mutation_frequency_chart: {os.path.getsize(path):,} bytes")

    path = comparison_heatmap(_make_comparisons(), str(OUTPUT_DIR / "legacy_heatmap.png"))
    assert os.path.exists(path)
    print(f"  comparison_heatmap: {os.path.getsize(path):,} bytes")

    path = risk_score_gauge(74.3, "TEST", str(OUTPUT_DIR / "legacy_gauge.png"))
    assert os.path.exists(path)
    print(f"  risk_score_gauge: {os.path.getsize(path):,} bytes")


if __name__ == "__main__":
    test_similarity_radar_chart()
    test_resistance_profile_heatmap()
    test_mutation_map()
    test_risk_gauge()
    test_similarity_bar_chart()
    test_lineage_distribution_nepal()
    test_save_chart()
    test_generate_all_charts()
    test_legacy_matplotlib_charts()
    print()
    print("ALL TESTS PASSED")
