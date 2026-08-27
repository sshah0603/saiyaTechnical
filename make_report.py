#!/usr/bin/env python3
"""Write report.html, a static snapshot of the analysis.

    python make_report.py

Purpose of this is to make the same numbers and same figures as the dashboard, but to put it in a single file that
opens in any browser. This is useful for sending to someone who is not going to use the GitHub and just needs the results. 
This is also what I typially give members of my lab who don't know how to code. 
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import analysis as an
import figures as fx

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "report.html"

FONTS = (
    '<link rel="preconnect" href="https://fonts.googleapis.com">'
    '<link href="https://fonts.googleapis.com/css2?'
    "family=Fraunces:opsz,wght@9..144,400;9..144,600&"
    "family=IBM+Plex+Mono:wght@400;500&"
    'family=IBM+Plex+Sans:wght@400;500;600&display=swap" rel="stylesheet">'
)

CSS = f"""
  :root {{
    --ink: {fx.INK}; --mute: {fx.MUTE}; --ground: {fx.GROUND};
    --card: {fx.CARD}; --rule: {fx.RULE};
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; background: var(--ground); color: var(--ink);
    font-family: 'IBM Plex Sans', system-ui, sans-serif; line-height: 1.55;
  }}
  .wrap {{ max-width: 1080px; margin: 0 auto; padding: 3rem 1.5rem 5rem; }}
  .eyebrow {{
    font-family: 'IBM Plex Mono', monospace; font-size: 0.72rem;
    letter-spacing: 0.16em; text-transform: uppercase; color: var(--mute);
  }}
  h1 {{
    font-family: 'Fraunces', Georgia, serif; font-weight: 600;
    font-size: clamp(2rem, 5vw, 3rem); letter-spacing: -0.015em;
    margin: 0.2rem 0 0.4rem;
  }}
  h2 {{
    font-family: 'Fraunces', Georgia, serif; font-weight: 600;
    font-size: 1.5rem; margin: 3rem 0 0.3rem;
  }}
  h3 {{ font-size: 1rem; font-weight: 600; margin: 1.8rem 0 0.5rem; }}
  p.lede {{ color: var(--mute); max-width: 64ch; margin: 0 0 1.2rem; }}
  .panel {{
    background: var(--card); border: 1px solid var(--rule);
    border-radius: 3px; padding: 1.1rem 1.25rem; margin: 1rem 0;
  }}
  .stats {{ display: flex; gap: 0.75rem; flex-wrap: wrap; margin: 1.2rem 0; }}
  .stats .panel {{ flex: 1 1 150px; margin: 0; }}
  .stat {{ font-family: 'IBM Plex Mono', monospace; font-size: 1.7rem; font-weight: 500; }}
  .stat-label {{
    font-family: 'IBM Plex Mono', monospace; font-size: 0.68rem;
    letter-spacing: 0.14em; text-transform: uppercase; color: var(--mute);
  }}
  .cols {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 0.75rem; }}
  table {{
    border-collapse: collapse; width: 100%;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.8rem;
  }}
  th {{
    text-align: left; font-weight: 500; color: var(--mute);
    text-transform: uppercase; letter-spacing: 0.08em; font-size: 0.66rem;
    border-bottom: 1px solid var(--rule); padding: 0.5rem 0.6rem;
  }}
  td {{ padding: 0.45rem 0.6rem; border-bottom: 1px solid var(--ground); }}
  tr:last-child td {{ border-bottom: none; }}
  .sig td {{ background: rgba(47, 168, 79, 0.07); }}
  footer {{
    margin-top: 4rem; padding-top: 1.2rem; border-top: 1px solid var(--rule);
    color: var(--mute); font-size: 0.82rem;
  }}
  @media (prefers-reduced-motion: reduce) {{ * {{ transition: none !important; }} }}
"""


def table_html(df: pd.DataFrame, highlight_col: str | None = None) -> str:
    head = "".join(f"<th>{c.replace('_', ' ')}</th>" for c in df.columns)
    rows = []
    for _, row in df.iterrows():
        cls = " class='sig'" if highlight_col and bool(row.get(highlight_col)) else ""
        cells = "".join(f"<td>{'' if pd.isna(v) else v}</td>" for v in row)
        rows.append(f"<tr{cls}>{cells}</tr>")
    return f"<div class='panel'><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"


def stats_html(items) -> str:
    blocks = "".join(
        f"<div class='panel'><div class='stat-label'>{k}</div><div class='stat'>{v}</div></div>"
        for k, v in items
    )
    return f"<div class='stats'>{blocks}</div>"


def main() -> None:
    conn = an.connect()
    try:
        freq = an.frequency_table(conn)
        annotated = an.annotated_frequencies(conn)
        cohort = an.responder_comparison(conn)
        base = an.baseline_cohort(conn)
    finally:
        conn.close()

    results = an.compare_responders(cohort, aggregate_by_subject=True)
    by_subject = (
        cohort.groupby(["subject", "response", "population"], as_index=False)["percentage"].mean()
    )
    cuts = an.baseline_breakdown(base)
    order = fx.population_order(annotated)

    fig_strips = fx.composition_strips(annotated)
    fig_box = fx.responder_boxplot(by_subject, order=order)
    fig_effect = fx.effect_size_bars(results)

    plots = [
        fig_strips.to_html(full_html=False, include_plotlyjs="cdn", config={"displaylogo": False}),
        fig_box.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False}),
        fig_effect.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False}),
    ]

    means = (
        annotated.groupby("population")["percentage"]
        .agg(["mean", "std", "min", "max"])
        .reindex(order)
        .round(2)
        .reset_index()
    )
    means["population"] = means["population"].map(fx.label_of)
    means.columns = ["population", "mean %", "sd", "min %", "max %"]

    stat_tbl = results.copy()
    stat_tbl["population"] = stat_tbl["population"].map(fx.label_of)
    for col in ["median_responder", "median_non_responder", "median_difference", "cliffs_delta"]:
        stat_tbl[col] = stat_tbl[col].round(2)
    for col in ["p_value", "p_value_adj", "welch_p_value"]:
        stat_tbl[col] = stat_tbl[col].map(lambda v: f"{v:.3g}" if pd.notna(v) else "")
    display_cols = [
        "population", "n_responder", "n_non_responder",
        "median_responder", "median_non_responder", "median_difference",
        "p_value", "p_value_adj", "welch_p_value", "cliffs_delta", "effect_size_label",
    ]

    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cell population frequencies</title>
{FONTS}<style>{CSS}</style></head><body><div class="wrap">

<div class="eyebrow">Loblaw Bio · clinical immunophenotyping</div>
<h1>Cell population frequencies</h1>
<p class="lede">Relative abundance of five immune populations across every sample in
the trial, the responder comparison for miraclib in melanoma, and the baseline
cohort breakdown. Generated from cell_counts.db.</p>

{stats_html([
    ("Samples", f"{annotated['sample'].nunique():,}"),
    ("Subjects", f"{annotated['subject'].nunique():,}"),
    ("Projects", f"{annotated['project'].nunique():,}"),
    ("Cells counted", f"{int(annotated.drop_duplicates('sample')['total_count'].sum()):,}"),
])}

<h2>Composition of every sample</h2>
<p class="lede">Each vertical strip is one sample, stacked to 100 percent.</p>
<div class="panel">{plots[0]}</div>
{table_html(means)}

<h2>Responders versus non-responders</h2>
<p class="lede">Melanoma, miraclib, PBMC only. Each subject's repeat samples are
averaged before testing, so one subject counts once.</p>
<div class="panel">{plots[1]}</div>
<div class="panel"><strong>Finding.</strong> {fx.headline(results)}</div>
{table_html(stat_tbl[display_cols].assign(significant=results['significant'].values), 'significant')}
<div class="panel">
<p class="lede" style="margin:0">Primary test is a two sided Mann-Whitney U, which
ranks the groups instead of assuming normality. p_value_adj applies a
Benjamini-Hochberg correction across the five populations and is the figure to
quote. The Welch t-test is a cross-check; agreement between the two is
reassuring, and disagreement means the result rests on few observations.
Cliff's delta reports how far the distributions are shifted, so a small p-value
paired with a negligible delta can be flagged as thin.</p>
</div>
<div class="panel">{plots[2]}</div>

<h2>Baseline cohort</h2>
<p class="lede">Melanoma PBMC samples at time zero from miraclib-treated patients.</p>
{stats_html([
    ("Samples", f"{len(base):,}"),
    ("Subjects", f"{base['subject'].nunique():,}"),
    ("Projects", f"{base['project'].nunique():,}"),
    ("Median age", f"{base.drop_duplicates('subject')['age'].median():.0f}"),
])}
<div class="cols">
  <div><h3>Samples per project</h3>{table_html(cuts['by_project'])}</div>
  <div><h3>Subjects by response</h3>{table_html(cuts['by_response'])}</div>
  <div><h3>Subjects by sex</h3>{table_html(cuts['by_sex'])}</div>
</div>

<footer>
Part 2 summary table holds {len(freq):,} rows, one per population per sample.
Run <code>streamlit run app.py</code> for the interactive version with filters and CSV export.
</footer>
</div></body></html>"""

    OUT.write_text(html, encoding="utf-8")
    print(f"Wrote {OUT.name} ({OUT.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
