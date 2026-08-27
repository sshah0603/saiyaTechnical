"""Design tokens and Plotly figure builders.

Imported by both app.py (interactive) and make_report.py (static export)
so the two never drift apart.

The palette is taken from the fluorophore channels a cytometrist reads at
the instrument: BV421 violet, FITC green, PE amber, APC red, PerCP teal.
One channel per population, held constant across every view, so a color
means the same cell type wherever it appears.
"""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go

INK = "#16202B"
MUTE = "#5C6B7A"
GROUND = "#EEF1F4"
CARD = "#FFFFFF"
RULE = "#D5DDE4"

CHANNEL = {
    "b_cell": "#6C4FD6",      # BV421
    "cd4_t_cell": "#2FA84F",  # FITC
    "cd8_t_cell": "#E8871A",  # PE
    "nk_cell": "#C2384A",     # APC
    "monocyte": "#17869B",    # PerCP
}

RESPONSE_COLOR = {"yes": "#2FA84F", "no": "#5C6B7A"}

LABEL = {
    "b_cell": "B cell",
    "cd4_t_cell": "CD4+ T",
    "cd8_t_cell": "CD8+ T",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}

BASE_LAYOUT = dict(
    paper_bgcolor=CARD,
    plot_bgcolor=CARD,
    font=dict(family="IBM Plex Sans, sans-serif", color=INK, size=12),
    margin=dict(l=48, r=16, t=36, b=40),
)


def population_order(df: pd.DataFrame) -> list[str]:
    present = set(df["population"])
    return [p for p in CHANNEL if p in present]


def label_of(population: str) -> str:
    return LABEL.get(population, population)


def composition_strips(df: pd.DataFrame, sort_by: str | None = None) -> go.Figure:
    """Stacked 100 percent bar, one thin strip per sample.

    The signature view. A cytometrist reads composition, not raw counts, so
    the whole cohort is shown as composition in a single image.
    """
    wide = df.pivot_table(
        index="sample", columns="population", values="percentage", aggfunc="mean"
    )
    if sort_by and sort_by in wide.columns:
        wide = wide.sort_values(sort_by, ascending=False)

    fig = go.Figure()
    for pop in population_order(df):
        if pop not in wide.columns:
            continue
        fig.add_bar(
            x=wide.index,
            y=wide[pop],
            name=label_of(pop),
            marker_color=CHANNEL[pop],
            marker_line_width=0,
            hovertemplate="%{x}<br>" + label_of(pop) + " %{y:.2f}%<extra></extra>",
        )
    fig.update_layout(
        **BASE_LAYOUT,
        barmode="stack",
        bargap=0.15,
        height=330,
        legend=dict(orientation="h", y=1.12, x=0, title=None),
        yaxis=dict(title="percent of sample", range=[0, 100], gridcolor=RULE),
        xaxis=dict(showticklabels=False, title=f"{len(wide)} samples"),
    )
    return fig


def responder_boxplot(df: pd.DataFrame, order: list[str] | None = None) -> go.Figure:
    """Grouped boxplot of population frequency, responder against non-responder.

    Individual observations are overlaid because group sizes are small
    enough that a box alone can hide how thin the evidence is.
    """
    order = order or population_order(df)
    fig = go.Figure()
    for resp in ["yes", "no"]:
        sub = df[df["response"] == resp]
        fig.add_box(
            x=sub["population"].map(label_of),
            y=sub["percentage"],
            name="Responder" if resp == "yes" else "Non-responder",
            marker_color=RESPONSE_COLOR[resp],
            line_width=1.4,
            boxpoints="all",
            jitter=0.45,
            pointpos=0,
            marker=dict(size=5, opacity=0.55),
            hovertemplate="%{y:.2f}%<extra></extra>",
        )
    fig.update_layout(
        **BASE_LAYOUT,
        boxmode="group",
        height=430,
        legend=dict(orientation="h", y=1.1, x=0, title=None),
        yaxis=dict(title="percent of sample", gridcolor=RULE, zeroline=False),
        xaxis=dict(
            title=None,
            categoryorder="array",
            categoryarray=[label_of(p) for p in order],
        ),
    )
    return fig


def effect_size_bars(results: pd.DataFrame) -> go.Figure:
    """Cliff's delta per population, so effect size is legible next to p-values."""
    res = results.dropna(subset=["cliffs_delta"]).copy()
    res = res.sort_values("cliffs_delta")
    colors = [
        CHANNEL.get(p, MUTE) if sig else RULE
        for p, sig in zip(res["population"], res["significant"])
    ]
    fig = go.Figure(
        go.Bar(
            x=res["cliffs_delta"],
            y=res["population"].map(label_of),
            orientation="h",
            marker_color=colors,
            marker_line_width=0,
            hovertemplate="%{y}: delta %{x:.2f}<extra></extra>",
        )
    )
    fig.update_layout(
        **BASE_LAYOUT,
        height=260,
        showlegend=False,
        xaxis=dict(
            title="Cliff's delta (responder minus non-responder)",
            range=[-1.05, 1.05],
            gridcolor=RULE,
            zeroline=True,
            zerolinecolor=INK,
            zerolinewidth=1,
        ),
        yaxis=dict(title=None),
    )
    return fig


def headline(results: pd.DataFrame, unit_noun: str = "subjects") -> str:
    """One sentence stating what the statistics support."""
    sig = results[results["significant"]]
    n_r = int(results["n_responder"].max())
    n_nr = int(results["n_non_responder"].max())

    if sig.empty:
        return (
            f"No population separates responders from non-responders once the five "
            f"tests are corrected for. n = {n_r} responder {unit_noun}, "
            f"{n_nr} non-responder."
        )

    parts = []
    for _, r in sig.iterrows():
        direction = "higher" if r["median_difference"] > 0 else "lower"
        parts.append(
            f"{label_of(r['population'])} is {direction} in responders by "
            f"{abs(r['median_difference']):.1f} percentage points "
            f"(adjusted p = {r['p_value_adj']:.2g})"
        )
    return f"n = {n_r} responder {unit_noun}, {n_nr} non-responder. " + "; ".join(parts) + "."
