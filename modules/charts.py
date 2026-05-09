"""
TBAnalytica Charts Module
Generates all visualizations: interactive (plotly) and static (matplotlib).
"""

from pathlib import Path
from typing import Optional
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from schema import (
    TBVariant, ComparisonResult, RiskScore, ClinicalReport,
    Mutation, DrugSensitivity,
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import plotly.graph_objects as go

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

_CONFIDENCE_COLORS = {
    "HIGH": "#2ecc71",
    "MODERATE": "#f1c40f",
    "LOW": "#e74c3c",
}

_RISK_ZONE_COLORS = [
    (0, 20, "#2ecc71", "Low"),
    (20, 40, "#f1c40f", "Moderate"),
    (40, 60, "#f39c12", "Moderate-High"),
    (60, 80, "#e67e22", "High"),
    (80, 100, "#e74c3c", "Critical"),
]

_MUTATION_COLORS = {
    "resistance": "#e74c3c",
    "silent_precursor": "#f39c12",
    "silent_benign": "#95a5a6",
    "novel": "#8e44ad",
}

_NEPAL_LINEAGE_DISTRIBUTION = {
    "Lineage 2 (Beijing)": 48.4,
    "Lineage 3 (CAS/Delhi)": 30.7,
    "Lineage 4 (Euro-American)": 14.5,
    "Lineage 1 (Indo-Oceanic)": 6.4,
}


# ---------------------------------------------------------------------------
# 1. similarity_radar_chart  (plotly)
# ---------------------------------------------------------------------------

def similarity_radar_chart(
    comparison: ComparisonResult,
    protein_score: float | None = None,
    gene_score: float | None = None,
    binding_score: float | None = None,
    resistance_score: float | None = None,
) -> go.Figure:
    """Radar/spider chart showing similarity across 4 dimensions.

    Scores can be passed explicitly or derived from ComparisonResult.
    All values are 0-100 percentages.
    """
    ps = protein_score if protein_score is not None else comparison.protein_similarity_score
    gs = gene_score if gene_score is not None else comparison.gene_similarity_score
    bs = binding_score if binding_score is not None else comparison.protein_similarity_score
    rs = resistance_score if resistance_score is not None else comparison.gene_similarity_score

    categories = [
        "Protein Sequence",
        "Gene Sequence",
        "Drug Binding Sites",
        "Resistance Profile",
    ]
    values = [ps, gs, bs, rs]
    values_closed = values + [values[0]]
    categories_closed = categories + [categories[0]]

    fig = go.Figure()

    fig.add_trace(go.Scatterpolar(
        r=values_closed,
        theta=categories_closed,
        fill="toself",
        fillcolor="rgba(41, 128, 185, 0.25)",
        line=dict(color="#2980b9", width=2),
        name=f"vs {comparison.matched_variant_id}",
        hovertemplate="%{theta}: %{r:.1f}%<extra></extra>",
    ))

    fig.add_trace(go.Scatterpolar(
        r=[100, 100, 100, 100, 100],
        theta=categories_closed,
        line=dict(color="#bdc3c7", width=1, dash="dot"),
        name="Perfect Match",
        hoverinfo="skip",
    ))

    fig.update_layout(
        polar=dict(
            radialaxis=dict(visible=True, range=[0, 100], ticksuffix="%"),
        ),
        title=dict(
            text=f"Similarity: {comparison.new_variant_id} vs {comparison.matched_variant_id}",
            x=0.5,
        ),
        showlegend=True,
        template="plotly_white",
        width=600,
        height=500,
    )

    return fig


# ---------------------------------------------------------------------------
# 2. resistance_profile_heatmap  (plotly)
# ---------------------------------------------------------------------------

def resistance_profile_heatmap(variants: list[TBVariant]) -> go.Figure:
    """Heatmap: variants (rows) vs drugs (columns).

    Green=sensitive, yellow=intermediate, red=resistant.
    """
    all_drugs: list[str] = []
    for v in variants:
        for drug in v.drug_resistance:
            if drug not in all_drugs:
                all_drugs.append(drug)

    if not all_drugs:
        all_drugs = ["isoniazid", "rifampicin", "ethambutol", "pyrazinamide",
                     "levofloxacin", "moxifloxacin", "bedaquiline", "linezolid"]

    status_map = {"sensitive": 0, "intermediate": 0.5, "resistant": 1}
    color_scale = [
        [0, "#2ecc71"],
        [0.5, "#f1c40f"],
        [1.0, "#e74c3c"],
    ]

    variant_labels = []
    z = []
    hover_text = []

    for v in variants:
        variant_labels.append(v.name or v.variant_id)
        row = []
        hover_row = []
        for drug in all_drugs:
            status = v.drug_resistance.get(drug)
            if status:
                val = status_map.get(status.value, 0)
                row.append(val)
                hover_row.append(f"{drug}: {status.value}")
            else:
                row.append(0)
                hover_row.append(f"{drug}: not tested")
        z.append(row)
        hover_text.append(hover_row)

    fig = go.Figure(data=go.Heatmap(
        z=z,
        x=all_drugs,
        y=variant_labels,
        colorscale=color_scale,
        zmin=0,
        zmax=1,
        text=hover_text,
        hovertemplate="%{text}<extra></extra>",
        showscale=False,
    ))

    fig.update_layout(
        title=dict(text="Drug Resistance Profile Comparison", x=0.5),
        xaxis=dict(title="Drug", tickangle=45),
        yaxis=dict(title="Variant", autorange="reversed"),
        template="plotly_white",
        width=max(600, len(all_drugs) * 80),
        height=max(400, len(variants) * 60 + 150),
    )

    return fig


# ---------------------------------------------------------------------------
# 3. mutation_map  (plotly)
# ---------------------------------------------------------------------------

def mutation_map(
    mutations: list[Mutation],
    gene_name: str,
    gene_length: int,
) -> go.Figure:
    """Linear map of a gene with mutations marked by type.

    Red=resistance, Orange=silent precursor, Grey=benign silent, Purple=novel.
    """
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=[0, gene_length],
        y=[0, 0],
        mode="lines",
        line=dict(color="#bdc3c7", width=8),
        hoverinfo="skip",
        showlegend=False,
    ))

    categories = {
        "Resistance": {"color": _MUTATION_COLORS["resistance"], "muts": []},
        "Silent Precursor": {"color": _MUTATION_COLORS["silent_precursor"], "muts": []},
        "Silent (Benign)": {"color": _MUTATION_COLORS["silent_benign"], "muts": []},
        "Novel": {"color": _MUTATION_COLORS["novel"], "muts": []},
    }

    for m in mutations:
        if m.is_resistance_conferring:
            categories["Resistance"]["muts"].append(m)
        elif m.is_synonymous and m.one_step_away_risk:
            categories["Silent Precursor"]["muts"].append(m)
        elif m.is_synonymous:
            categories["Silent (Benign)"]["muts"].append(m)
        else:
            categories["Novel"]["muts"].append(m)

    for cat_name, info in categories.items():
        if not info["muts"]:
            continue
        positions = [m.position for m in info["muts"]]
        labels = [m.short_code() for m in info["muts"]]
        drugs = [m.drug_affected or "" for m in info["muts"]]
        hover = [
            f"{l}<br>Drug: {d}" if d else l
            for l, d in zip(labels, drugs)
        ]

        fig.add_trace(go.Scatter(
            x=positions,
            y=[0] * len(positions),
            mode="markers+text",
            marker=dict(
                color=info["color"],
                size=14,
                symbol="diamond",
                line=dict(width=1, color="white"),
            ),
            text=labels,
            textposition="top center",
            textfont=dict(size=9),
            hovertext=hover,
            hoverinfo="text",
            name=cat_name,
        ))

    fig.update_layout(
        title=dict(text=f"Mutation Map: {gene_name} ({gene_length} codons)", x=0.5),
        xaxis=dict(title="Codon Position", range=[-gene_length * 0.05, gene_length * 1.05]),
        yaxis=dict(visible=False, range=[-0.5, 1.5]),
        template="plotly_white",
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        width=max(700, gene_length // 2),
        height=350,
    )

    return fig


# ---------------------------------------------------------------------------
# 4. risk_gauge  (plotly)
# ---------------------------------------------------------------------------

def risk_gauge(risk_score: RiskScore) -> go.Figure:
    """Gauge/speedometer chart for risk score (0-100)."""
    steps = []
    for low, high, color, label in _RISK_ZONE_COLORS:
        steps.append(dict(range=[low, high], color=color, name=label))

    bar_color = "#2c3e50"
    for low, high, color, _ in _RISK_ZONE_COLORS:
        if low <= risk_score.score < high or (high == 100 and risk_score.score == 100):
            bar_color = color
            break

    fig = go.Figure(go.Indicator(
        mode="gauge+number+delta",
        value=risk_score.score,
        number=dict(suffix="%", font=dict(size=36)),
        title=dict(
            text=f"Risk Level: {risk_score.level.value}",
            font=dict(size=16),
        ),
        gauge=dict(
            axis=dict(range=[0, 100], tickwidth=1, tickcolor="#2c3e50"),
            bar=dict(color=bar_color, thickness=0.3),
            steps=steps,
            threshold=dict(
                line=dict(color="#2c3e50", width=4),
                thickness=0.8,
                value=risk_score.score,
            ),
        ),
    ))

    fig.update_layout(
        template="plotly_white",
        width=500,
        height=350,
    )

    return fig


# ---------------------------------------------------------------------------
# 5. similarity_bar_chart  (plotly)
# ---------------------------------------------------------------------------

def similarity_bar_chart(comparison_results: list[ComparisonResult]) -> go.Figure:
    """Horizontal bar chart of top 5 closest known variants by similarity."""
    top5 = sorted(comparison_results, key=lambda c: c.weighted_final_score, reverse=True)[:5]
    top5.reverse()

    labels = [c.matched_variant_id for c in top5]
    scores = [c.weighted_final_score for c in top5]
    colors = [_CONFIDENCE_COLORS.get(c.confidence_level.value, "#3498db") for c in top5]
    conf_labels = [c.confidence_level.value for c in top5]

    fig = go.Figure()

    fig.add_trace(go.Bar(
        y=labels,
        x=scores,
        orientation="h",
        marker=dict(color=colors, line=dict(width=1, color="white")),
        text=[f"{s:.1f}% ({cl})" for s, cl in zip(scores, conf_labels)],
        textposition="outside",
        hovertemplate="%{y}: %{x:.1f}%<extra></extra>",
    ))

    fig.update_layout(
        title=dict(text="Top 5 Closest Known Variants", x=0.5),
        xaxis=dict(title="Weighted Similarity (%)", range=[0, 110]),
        yaxis=dict(title=""),
        template="plotly_white",
        showlegend=False,
        width=700,
        height=max(350, len(top5) * 60 + 150),
    )

    return fig


# ---------------------------------------------------------------------------
# 6. lineage_distribution_nepal  (plotly)
# ---------------------------------------------------------------------------

def lineage_distribution_nepal(
    variants: list[TBVariant],
    highlight_lineage: str | None = None,
) -> go.Figure:
    """Donut chart of Nepal TB lineage distribution.

    Highlights which lineage the current patient belongs to.
    """
    labels = list(_NEPAL_LINEAGE_DISTRIBUTION.keys())
    values = list(_NEPAL_LINEAGE_DISTRIBUTION.values())

    base_colors = ["#e74c3c", "#3498db", "#2ecc71", "#f39c12"]

    if highlight_lineage is None and variants:
        for v in variants:
            if v.lineage:
                highlight_lineage = v.lineage
                break

    pull = [0] * len(labels)
    final_colors = list(base_colors)
    if highlight_lineage:
        hl = highlight_lineage.lower()
        matched_idx = None
        for i, label in enumerate(labels):
            ll = label.lower()
            if hl == ll or hl in ll or ll in hl:
                matched_idx = i
                break
        if matched_idx is None:
            hl_words = set(hl.replace("(", " ").replace(")", " ").split())
            best_overlap = 0
            for i, label in enumerate(labels):
                lw = set(label.lower().replace("(", " ").replace(")", " ").split())
                overlap = len(hl_words & lw)
                if overlap > best_overlap:
                    best_overlap = overlap
                    matched_idx = i
        if matched_idx is not None:
            pull[matched_idx] = 0.1

    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        hole=0.45,
        pull=pull,
        marker=dict(colors=final_colors, line=dict(width=2, color="white")),
        textinfo="label+percent",
        textposition="outside",
        hovertemplate="%{label}: %{value}%<extra></extra>",
    )])

    annotation_text = f"Patient: {highlight_lineage}" if highlight_lineage else "Nepal TB"
    if len(annotation_text) > 25:
        annotation_text = annotation_text[:22] + "..."

    fig.update_layout(
        title=dict(text="TB Lineage Distribution in Nepal", x=0.5),
        annotations=[dict(text=annotation_text, x=0.5, y=0.5, font_size=11, showarrow=False)],
        template="plotly_white",
        showlegend=True,
        width=650,
        height=500,
    )

    return fig


# ---------------------------------------------------------------------------
# 7. save_chart
# ---------------------------------------------------------------------------

def save_chart(
    fig: go.Figure,
    output_path: str,
    format: str = "png",
) -> str:
    """Save a plotly figure as PNG (static) or HTML (interactive).

    Returns the saved file path.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if format.lower() == "html":
        if not str(path).endswith(".html"):
            path = path.with_suffix(".html")
        fig.write_html(str(path), include_plotlyjs="cdn")
    else:
        if not str(path).endswith(".png"):
            path = path.with_suffix(".png")
        try:
            fig.write_image(str(path), scale=2)
        except (ValueError, ImportError):
            fig.write_html(str(path).replace(".png", ".html"), include_plotlyjs="cdn")
            return str(path).replace(".png", ".html")

    return str(path)


# ---------------------------------------------------------------------------
# 8. generate_all_charts
# ---------------------------------------------------------------------------

def generate_all_charts(
    report: ClinicalReport,
    comparison_results: list[ComparisonResult],
    output_dir: str,
) -> dict:
    """Generate all relevant charts for a report.

    Returns dict of chart_name: file_path.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    vid = report.variant.variant_id
    charts: dict[str, str] = {}

    if report.comparison_result and report.comparison_result.matched_variant_id:
        fig = similarity_radar_chart(report.comparison_result)
        charts["similarity_radar"] = save_chart(fig, str(out / f"{vid}_similarity_radar.html"), "html")

    if comparison_results:
        fig = similarity_bar_chart(comparison_results)
        charts["similarity_bar"] = save_chart(fig, str(out / f"{vid}_similarity_bar.html"), "html")

    variants_for_heatmap = [report.variant]
    if comparison_results:
        from modules.variant_db import get_variant_by_id
        for cr in comparison_results[:4]:
            matched = get_variant_by_id(cr.matched_variant_id)
            if matched:
                variants_for_heatmap.append(matched)
    if len(variants_for_heatmap) > 1 or report.variant.drug_resistance:
        fig = resistance_profile_heatmap(variants_for_heatmap)
        charts["resistance_heatmap"] = save_chart(fig, str(out / f"{vid}_resistance_heatmap.html"), "html")

    for rg in report.variant.resistance_genes:
        if rg.mutations:
            gene_len = 1000
            max_pos = max(m.position for m in rg.mutations)
            if max_pos > 800:
                gene_len = max_pos + 200
            fig = mutation_map(rg.mutations, rg.gene_name, gene_len)
            charts[f"mutation_map_{rg.gene_name}"] = save_chart(
                fig, str(out / f"{vid}_{rg.gene_name}_mutation_map.html"), "html"
            )

    fig = risk_gauge(report.risk_score)
    charts["risk_gauge"] = save_chart(fig, str(out / f"{vid}_risk_gauge.html"), "html")

    fig = lineage_distribution_nepal([report.variant])
    charts["lineage_distribution"] = save_chart(fig, str(out / f"{vid}_lineage_nepal.html"), "html")

    return charts


# ---------------------------------------------------------------------------
# Legacy matplotlib functions (kept for main.py backward compatibility)
# ---------------------------------------------------------------------------

def resistance_profile_chart(variant: TBVariant, output_path: Optional[str] = None) -> str:
    if not output_path:
        output_path = str(OUTPUT_DIR / f"{variant.variant_id}_resistance_profile.png")

    resistant = variant.resistant_drugs()
    susceptible = variant.susceptible_drugs()
    drugs = resistant + susceptible
    status = ["Resistant"] * len(resistant) + ["Susceptible"] * len(susceptible)
    colors = ["#e74c3c" if s == "Resistant" else "#2ecc71" for s in status]

    fig, ax = plt.subplots(figsize=(10, max(6, len(drugs) * 0.5)))
    y_pos = range(len(drugs))
    ax.barh(y_pos, [1] * len(drugs), color=colors, edgecolor="white", height=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(drugs, fontsize=10)
    ax.set_xlim(0, 1.2)
    ax.set_xticks([])
    ax.set_title(f"Drug Resistance Profile: {variant.name}", fontsize=13, fontweight="bold")

    resistant_patch = mpatches.Patch(color="#e74c3c", label="Resistant")
    susceptible_patch = mpatches.Patch(color="#2ecc71", label="Susceptible")
    ax.legend(handles=[resistant_patch, susceptible_patch], loc="lower right")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    return output_path


def mutation_frequency_chart(variant: TBVariant, output_path: Optional[str] = None) -> str:
    if not output_path:
        output_path = str(OUTPUT_DIR / f"{variant.variant_id}_mutation_freq.png")

    all_mutations = []
    for rg in variant.resistance_genes:
        for m in rg.mutations:
            all_mutations.append((rg.gene_name, m))

    if not all_mutations:
        return ""

    labels = [f"{gene}\n{m.short_code()}" for gene, m in all_mutations]
    colors = []
    for _, m in all_mutations:
        if m.is_synonymous:
            colors.append("#f39c12")
        elif m.resistance_level.value == "high":
            colors.append("#e74c3c")
        elif m.resistance_level.value == "moderate":
            colors.append("#e67e22")
        else:
            colors.append("#3498db")

    bar_values = [1.0 if m.is_resistance_conferring else 0.5 for _, m in all_mutations]

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(range(len(labels)), bar_values, color=colors, edgecolor="white", width=0.6)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=9, rotation=45, ha="right")
    ax.set_ylabel("Resistance Impact", fontsize=11)
    ax.set_title(f"Mutations: {variant.name}", fontsize=13, fontweight="bold")
    ax.set_ylim(0, 1.2)

    for bar, (_, m) in zip(bars, all_mutations):
        label = "R" if m.is_resistance_conferring else ("S" if m.is_synonymous else "?")
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                label, ha="center", va="bottom", fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    return output_path


def comparison_heatmap(comparisons: list[ComparisonResult], output_path: Optional[str] = None) -> str:
    if not output_path:
        output_path = str(OUTPUT_DIR / "comparison_heatmap.png")

    if not comparisons:
        return ""

    labels = [c.matched_variant_id for c in comparisons]
    shared_counts = [len(c.resistance_mutations) for c in comparisons]
    scores = [c.weighted_final_score for c in comparisons]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(6, len(labels) * 0.5)))

    y_pos = range(len(labels))
    ax1.barh(y_pos, shared_counts, color="#3498db", edgecolor="white", height=0.6)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(labels, fontsize=9)
    ax1.set_xlabel("Resistance Mutations")
    ax1.set_title("Resistance Mutation Overlap", fontsize=11, fontweight="bold")

    colors = ["#e74c3c" if s > 60 else "#f39c12" if s > 30 else "#2ecc71" for s in scores]
    ax2.barh(y_pos, scores, color=colors, edgecolor="white", height=0.6)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(labels, fontsize=9)
    ax2.set_xlabel("Weighted Score (%)")
    ax2.set_title("Similarity Score", fontsize=11, fontweight="bold")
    ax2.set_xlim(0, 100)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    return output_path


def treatment_availability_chart(rec: dict, output_path: Optional[str] = None) -> str:
    if not output_path:
        output_path = str(OUTPUT_DIR / f"{rec['variant_id']}_treatment_avail.png")

    all_drugs = list(dict.fromkeys(rec["recommended_regimen"] + rec["alternative_regimen"]))
    if not all_drugs:
        return ""

    availability = []
    drug_type = []
    for drug in all_drugs:
        avail = rec["nepal_availability"].get(drug, {})
        availability.append(1 if avail.get("available", False) else 0)
        drug_type.append("Primary" if drug in rec["recommended_regimen"] else "Alternative")

    fig, ax = plt.subplots(figsize=(10, max(6, len(all_drugs) * 0.4)))

    y_pos = range(len(all_drugs))
    colors = []
    for avail, dtype in zip(availability, drug_type):
        if not avail:
            colors.append("#e74c3c")
        elif dtype == "Primary":
            colors.append("#2ecc71")
        else:
            colors.append("#3498db")

    ax.barh(y_pos, [1] * len(all_drugs), color=colors, edgecolor="white", height=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(all_drugs, fontsize=10)
    ax.set_xlim(0, 1.3)
    ax.set_xticks([])

    contraindicated = set(rec["contraindicated_drugs"])
    for i, drug in enumerate(all_drugs):
        if drug in contraindicated:
            ax.text(1.05, i, "CONTRAINDICATED", fontsize=8, color="#e74c3c", fontweight="bold", va="center")

    ax.set_title(f"Treatment Availability (Nepal): {rec['variant_id']}", fontsize=12, fontweight="bold")

    patches = [
        mpatches.Patch(color="#2ecc71", label="Available (Primary)"),
        mpatches.Patch(color="#3498db", label="Available (Alternative)"),
        mpatches.Patch(color="#e74c3c", label="Not Available"),
    ]
    ax.legend(handles=patches, loc="lower right", fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    return output_path


def risk_score_gauge(risk_score: float, variant_id: str, output_path: Optional[str] = None) -> str:
    if not output_path:
        output_path = str(OUTPUT_DIR / f"{variant_id}_risk_gauge.png")

    fig, ax = plt.subplots(figsize=(8, 4), subplot_kw={"projection": "polar"})

    for i in range(100):
        angle = np.pi - (i / 100) * np.pi
        if i < 20:
            color = "#2ecc71"
        elif i < 40:
            color = "#f1c40f"
        elif i < 60:
            color = "#f39c12"
        elif i < 80:
            color = "#e67e22"
        else:
            color = "#e74c3c"
        ax.bar(angle, 1, width=np.pi/100, bottom=0.5, color=color, alpha=0.8)

    needle_angle = np.pi - (risk_score / 100) * np.pi
    ax.plot([needle_angle, needle_angle], [0, 1.3], color="black", linewidth=2)
    ax.plot(needle_angle, 1.3, "ko", markersize=5)

    ax.set_ylim(0, 1.8)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.spines["polar"].set_visible(False)
    ax.grid(False)

    ax.text(np.pi/2, -0.3, f"{risk_score:.1f}%", ha="center", va="center",
            fontsize=20, fontweight="bold", transform=ax.transAxes)
    ax.text(np.pi/2, -0.15, "Resistance Risk Score", ha="center", va="center",
            fontsize=11, transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    return output_path
