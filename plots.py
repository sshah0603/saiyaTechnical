#!/usr/bin/env python3
"""Static matplotlib plots written to disk by the pipeline.

Plotly drives the dashboard because it is interactive. Static PNGs need a
headless browser for Plotly export, which is a heavy dependency to add for
a grading run, so the on-disk plots use matplotlib instead. Both read their
colors from figures.py, so a population is the same color everywhere.

Called by run_pipeline.py. Also runnable on its own:

    python plots.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

import analysis as an
import figures as fx

ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"
DOCS = ROOT / "docs"

AVAILABLE = {f.name for f in font_manager.fontManager.ttflist}
BODY = next((f for f in ["IBM Plex Sans", "DejaVu Sans"] if f in AVAILABLE), "sans-serif")

plt.rcParams.update(
    {
        "font.family": BODY,
        "font.size": 9,
        "figure.facecolor": fx.CARD,
        "axes.facecolor": fx.CARD,
        "axes.edgecolor": fx.RULE,
        "axes.labelcolor": fx.INK,
        "text.color": fx.INK,
        "xtick.color": fx.MUTE,
        "ytick.color": fx.MUTE,
        "axes.grid": True,
        "grid.color": fx.RULE,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.facecolor": fx.CARD,
        "savefig.bbox": "tight",
        "savefig.dpi": 160,
    }
)


def composition(annotated, order, out_dir: Path):
    wide = annotated.pivot_table(
        index="sample", columns="population", values="percentage", aggfunc="mean"
    )
    # Sort by whichever population spreads the most. Sorting by a population
    # that barely varies produces a flat block that shows nothing.
    sort_key = wide[order].std().idxmax()
    wide = wide.sort_values(sort_key, ascending=False)

    fig, ax = plt.subplots(figsize=(9, 2.9))
    ax.set_axisbelow(True)
    bottom = 0
    for pop in order:
        ax.bar(
            range(len(wide)),
            wide[pop].to_numpy(),
            bottom=bottom,
            width=1.0,
            color=fx.CHANNEL[pop],
            label=fx.label_of(pop),
            linewidth=0,
        )
        bottom = bottom + wide[pop].to_numpy()

    ax.set_xlim(-0.5, len(wide) - 0.5)
    ax.set_ylim(0, 100)
    ax.set_xticks([])
    ax.set_ylabel("percent of sample")
    ax.set_xlabel(f"{len(wide)} samples, sorted by {fx.label_of(sort_key)}")
    ax.grid(axis="x", visible=False)
    ax.legend(
        ncol=5, frameon=False, loc="lower center",
        bbox_to_anchor=(0.5, 1.01), fontsize=8,
    )
    fig.savefig(out_dir / "part2_composition.png")
    plt.close(fig)


def responder_box(by_subject, order, results, out_dir: Path):
    fig, ax = plt.subplots(figsize=(9, 3.6))
    ax.set_axisbelow(True)
    width = 0.34
    sig = set(results.loc[results["significant"], "population"])

    for offset, resp in [(-width / 2, "yes"), (width / 2, "no")]:
        data = [
            by_subject.loc[
                (by_subject["population"] == p) & (by_subject["response"] == resp),
                "percentage",
            ].to_numpy()
            for p in order
        ]
        color = fx.RESPONSE_COLOR[resp]
        bp = ax.boxplot(
            data,
            positions=[i + offset for i in range(len(order))],
            widths=width * 0.86,
            patch_artist=True,
            medianprops=dict(color=fx.INK, linewidth=1.2),
            whiskerprops=dict(color=color, linewidth=1.1),
            capprops=dict(color=color, linewidth=1.1),
            flierprops=dict(markersize=3, markerfacecolor=color, markeredgecolor="none"),
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.28)
            patch.set_edgecolor(color)
            patch.set_linewidth(1.2)
        ax.plot([], [], color=color, linewidth=6, alpha=0.5,
                label="Responder" if resp == "yes" else "Non-responder")

    labels = [
        f"{fx.label_of(p)}{' *' if p in sig else ''}" for p in order
    ]
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels)
    ax.set_ylabel("percent of sample")
    ax.grid(axis="x", visible=False)
    ax.legend(frameon=False, ncol=2, loc="lower center",
              bbox_to_anchor=(0.5, 1.01), fontsize=8)
    ax.text(
        0.995, 0.02, "* adjusted p < 0.05", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=7.5, color=fx.MUTE,
    )
    fig.savefig(out_dir / "part3_responder_boxplot.png")
    plt.close(fig)


def per_population_strip(annotated, order, out_dir: Path):
    """Small multiples: one panel per population, frequency across all samples.

    The stacked view answers "what is this sample made of". This answers
    "how much does this population vary", which is the question that decides
    whether a population is worth testing at all.
    """
    fig, axes = plt.subplots(1, len(order), figsize=(10, 2.6), sharey=False)
    for ax, pop in zip(axes, order):
        vals = annotated.loc[annotated["population"] == pop, "percentage"]
        ax.hist(vals, bins=28, color=fx.CHANNEL[pop], linewidth=0)
        ax.set_title(fx.label_of(pop), fontsize=9)
        ax.set_yticks([])
        ax.grid(axis="x", visible=False)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("samples")
    fig.supxlabel("percent of sample", fontsize=9, color=fx.MUTE)
    fig.tight_layout()
    fig.savefig(out_dir / "part2_population_distributions.png")
    plt.close(fig)


def write_all(out_dirs: list[Path]) -> list[Path]:
    """Generate every static plot into each directory given."""
    conn = an.connect()
    try:
        annotated = an.annotated_frequencies(conn)
        cohort = an.responder_comparison(conn)
    finally:
        conn.close()

    order = fx.population_order(annotated)
    results = an.compare_responders(cohort, aggregate_by_subject=True)
    by_subject = cohort.groupby(
        ["subject", "response", "population"], as_index=False
    )["percentage"].mean()

    written = []
    for out_dir in out_dirs:
        out_dir.mkdir(parents=True, exist_ok=True)
        composition(annotated, order, out_dir)
        per_population_strip(annotated, order, out_dir)
        if not by_subject.empty:
            responder_box(by_subject, order, results, out_dir)
        written.extend(sorted(out_dir.glob("*.png")))
    return written


def main() -> None:
    for path in write_all([OUTPUTS, DOCS]):
        print(f"  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
