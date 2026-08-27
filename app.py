"""Interactive dashboard for the Loblaw Bio cell-count analysis.

    streamlit run app.py

Reads cell_counts.db. Run `python load_data.py` first.
"""

from __future__ import annotations

import io

import pandas as pd
import streamlit as st

import analysis as an

import figures as fx
from figures import CARD, CHANNEL, GROUND, INK, LABEL, MUTE, RULE

st.set_page_config(page_title="Cell populations", layout="wide", page_icon="◍")

st.markdown(
    f"""
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,400;9..144,600&family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&display=swap" rel="stylesheet">
    <style>
      .stApp {{ background: {GROUND}; }}
      html, body, [class*="css"] {{
          font-family: 'IBM Plex Sans', system-ui, sans-serif;
          color: {INK};
      }}
      h1, h2, h3 {{ font-family: 'Fraunces', Georgia, serif; font-weight: 600; letter-spacing: -0.01em; }}
      h1 {{ font-size: 2.4rem; margin-bottom: 0.1rem; }}
      .eyebrow {{
          font-family: 'IBM Plex Mono', monospace; font-size: 0.72rem;
          letter-spacing: 0.16em; text-transform: uppercase; color: {MUTE};
      }}
      .lede {{ color: {MUTE}; font-size: 0.95rem; max-width: 62ch; }}
      .panel {{
          background: {CARD}; border: 1px solid {RULE}; border-radius: 3px;
          padding: 1.1rem 1.25rem; margin-bottom: 0.9rem;
      }}
      .stat {{ font-family: 'IBM Plex Mono', monospace; font-size: 1.7rem; font-weight: 500; }}
      .stat-label {{
          font-family: 'IBM Plex Mono', monospace; font-size: 0.68rem;
          letter-spacing: 0.14em; text-transform: uppercase; color: {MUTE};
      }}
      .stDataFrame {{ font-family: 'IBM Plex Mono', monospace; }}
      .stTabs [data-baseweb="tab"] {{
          font-family: 'IBM Plex Mono', monospace; font-size: 0.78rem;
          letter-spacing: 0.1em; text-transform: uppercase;
      }}
      hr {{ border-color: {RULE}; }}
    </style>
    """,
    unsafe_allow_html=True,
)



# ---------------------------------------------------------------- data

@st.cache_data(show_spinner=False)
def load_frames():
    conn = an.connect()
    try:
        return an.frequency_table(conn), an.annotated_frequencies(conn)
    finally:
        conn.close()


@st.cache_data(show_spinner=False)
def load_cohort(condition, treatment, sample_type):
    conn = an.connect()
    try:
        return an.responder_comparison(conn, condition, treatment, sample_type)
    finally:
        conn.close()


@st.cache_data(show_spinner=False)
def load_baseline(condition, treatment, sample_type, timepoint):
    conn = an.connect()
    try:
        return an.baseline_cohort(conn, condition, treatment, sample_type, timepoint)
    finally:
        conn.close()


def download(df: pd.DataFrame, label: str, filename: str, key: str):
    buf = io.StringIO()
    df.to_csv(buf, index=False)
    st.download_button(label, buf.getvalue(), file_name=filename, mime="text/csv", key=key)


def stat_block(items):
    cols = st.columns(len(items))
    for col, (label, value) in zip(cols, items):
        col.markdown(
            f"<div class='panel'><div class='stat-label'>{label}</div>"
            f"<div class='stat'>{value}</div></div>",
            unsafe_allow_html=True,
        )


try:
    freq, annotated = load_frames()
except FileNotFoundError as exc:
    st.error(str(exc))
    st.stop()

ORDER = [p for p in CHANNEL if p in set(annotated["population"])]


# ---------------------------------------------------------------- header

st.markdown("<div class='eyebrow'>Loblaw Bio · clinical immunophenotyping</div>",
            unsafe_allow_html=True)
st.markdown("# Cell population frequencies")
st.markdown(
    "<p class='lede'>Relative abundance of five immune populations across every "
    "sample in the trial, with the responder comparison for miraclib in melanoma "
    "and the baseline cohort breakdown.</p>",
    unsafe_allow_html=True,
)

tab_overview, tab_response, tab_baseline, tab_data = st.tabs(
    ["Overview", "Response", "Baseline cohort", "Data"]
)


# ---------------------------------------------------------------- Part 2

with tab_overview:
    st.markdown("### Composition of every sample")
    st.markdown(
        "<p class='lede'>Each vertical strip is one sample, stacked to 100 percent. "
        "Sort by a population to see where the cohort separates.</p>",
        unsafe_allow_html=True,
    )

    left, mid, right = st.columns([1, 1, 1])
    conditions = sorted(annotated["condition"].dropna().unique())
    sample_types = sorted(annotated["sample_type"].dropna().unique())

    pick_conditions = left.multiselect("Condition", conditions, default=conditions)
    pick_types = mid.multiselect("Sample type", sample_types, default=sample_types)
    sort_by = right.selectbox(
        "Sort strips by", ["sample id"] + [LABEL[p] for p in ORDER], index=0
    )

    view = annotated[
        annotated["condition"].isin(pick_conditions)
        & annotated["sample_type"].isin(pick_types)
    ]

    if view.empty:
        st.info("No samples match those filters. Widen the selection.")
    else:
        sort_key = None
        if sort_by != "sample id":
            sort_key = {v: k for k, v in LABEL.items()}[sort_by]
        fig = fx.composition_strips(view, sort_by=sort_key)
        st.plotly_chart(fig, width="stretch")

        stat_block(
            [
                ("Samples", f"{view['sample'].nunique():,}"),
                ("Subjects", f"{view['subject'].nunique():,}"),
                ("Projects", f"{view['project'].nunique():,}"),
                ("Cells counted", f"{int(view.drop_duplicates('sample')['total_count'].sum()):,}"),
            ]
        )

        st.markdown("#### Mean frequency by population")
        means = (
            view.groupby("population")["percentage"]
            .agg(["mean", "std", "min", "max"])
            .reindex(ORDER)
            .round(2)
            .reset_index()
        )
        means["population"] = means["population"].map(LABEL)
        means.columns = ["population", "mean %", "sd", "min %", "max %"]
        st.dataframe(means, width="stretch", hide_index=True)


# ---------------------------------------------------------------- Part 3

with tab_response:
    st.markdown("### Responders versus non-responders")

    c1, c2, c3, c4 = st.columns([1, 1, 1, 1.3])
    cond = c1.selectbox("Condition", conditions,
                        index=conditions.index("melanoma") if "melanoma" in conditions else 0)
    treatments = sorted(annotated["treatment"].dropna().unique())
    treat = c2.selectbox("Treatment", treatments,
                         index=treatments.index("miraclib") if "miraclib" in treatments else 0)
    stype = c3.selectbox("Sample type", sample_types,
                         index=sample_types.index("PBMC") if "PBMC" in sample_types else 0)
    unit = c4.radio(
        "Unit of analysis",
        ["Subject mean", "Each sample"],
        horizontal=True,
        help=(
            "Subjects contribute several samples across timepoints and replicates. "
            "Treating each sample as independent inflates n and overstates significance."
        ),
    )

    cohort = load_cohort(cond, treat, stype)

    if cohort.empty:
        st.info("No samples with a recorded response match that combination.")
    else:
        by_subject = unit == "Subject mean"
        results = an.compare_responders(cohort, aggregate_by_subject=by_subject)

        plot_df = (
            cohort.groupby(["subject", "response", "population"], as_index=False)["percentage"].mean()
            if by_subject
            else cohort[["sample", "response", "population", "percentage"]]
        )

        fig = fx.responder_boxplot(plot_df, order=ORDER)
        st.plotly_chart(fig, width="stretch")

        noun = "subjects" if by_subject else "samples"
        headline = fx.headline(results, unit_noun=noun)
        st.markdown(f"<div class='panel'>{headline}</div>", unsafe_allow_html=True)

        st.plotly_chart(fx.effect_size_bars(results), width="stretch")

        show = results.copy()
        show["population"] = show["population"].map(LABEL).fillna(show["population"])
        for col in ["median_responder", "median_non_responder", "median_difference", "cliffs_delta"]:
            show[col] = show[col].round(3)
        for col in ["p_value", "p_value_adj", "welch_p_value"]:
            show[col] = show[col].map(lambda v: f"{v:.3g}" if pd.notna(v) else "")
        st.dataframe(show, width="stretch", hide_index=True)

        with st.expander("How to read this table"):
            st.markdown(
                """
- **p_value** is a two sided Mann-Whitney U test. It ranks the two groups rather
  than assuming normal distributions, which suits frequency data at these group sizes.
- **p_value_adj** applies a Benjamini-Hochberg correction across the five
  populations. Five tests at 0.05 gives roughly a one in four chance of at least
  one false positive if uncorrected, so the adjusted column is the one to quote.
- **welch_p_value** is a Welch t-test run as a cross-check. Agreement between the
  two tests is reassuring. Disagreement means the result rests on a few
  observations and should be described as preliminary.
- **cliffs_delta** measures how far the distributions are shifted, from -1 to 1.
  A small p-value with a negligible delta means a real but clinically thin difference.
"""
            )

        download(results, "Download statistics (CSV)", "responder_statistics.csv", "dl_stats")


# ---------------------------------------------------------------- Part 4

with tab_baseline:
    st.markdown("### Baseline cohort")
    st.markdown(
        "<p class='lede'>Samples drawn before dosing, for the selected condition, "
        "treatment and sample type.</p>",
        unsafe_allow_html=True,
    )

    b1, b2, b3, b4 = st.columns(4)
    b_cond = b1.selectbox("Condition ", conditions,
                          index=conditions.index("melanoma") if "melanoma" in conditions else 0)
    b_treat = b2.selectbox("Treatment ", treatments,
                           index=treatments.index("miraclib") if "miraclib" in treatments else 0)
    b_type = b3.selectbox("Sample type ", sample_types,
                          index=sample_types.index("PBMC") if "PBMC" in sample_types else 0)
    timepoints = sorted(annotated["time_from_treatment_start"].dropna().unique().tolist())
    b_time = b4.selectbox("Timepoint (days)", timepoints,
                          index=timepoints.index(0) if 0 in timepoints else 0)

    base = load_baseline(b_cond, b_treat, b_type, int(b_time))

    if base.empty:
        st.info("No samples match that combination.")
    else:
        cuts = an.baseline_breakdown(base)
        stat_block(
            [
                ("Samples", f"{len(base):,}"),
                ("Subjects", f"{base['subject'].nunique():,}"),
                ("Projects", f"{base['project'].nunique():,}"),
                ("Median age", f"{base.drop_duplicates('subject')['age'].median():.0f}"),
            ]
        )

        c1, c2, c3 = st.columns(3)
        for col, (key, title) in zip(
            (c1, c2, c3),
            [("by_project", "Samples per project"),
             ("by_response", "Subjects by response"),
             ("by_sex", "Subjects by sex")],
        ):
            col.markdown(f"**{title}**")
            col.dataframe(cuts[key], width="stretch", hide_index=True)

        st.markdown("#### Cohort samples")
        st.dataframe(base, width="stretch", hide_index=True, height=300)
        download(base, "Download cohort (CSV)", "baseline_cohort.csv", "dl_base")


# ---------------------------------------------------------------- data

with tab_data:
    st.markdown("### Part 2 summary table")
    st.markdown(
        "<p class='lede'>One row per population per sample: total count, count and "
        "relative frequency. This is the table every other view is built from.</p>",
        unsafe_allow_html=True,
    )
    search = st.text_input("Filter by sample id", placeholder="e.g. s0012")
    table = freq[freq["sample"].str.contains(search, case=False, na=False)] if search else freq
    st.dataframe(table, width="stretch", hide_index=True, height=480)
    st.caption(f"{len(table):,} rows")
    download(freq, "Download summary table (CSV)", "cell_frequencies.csv", "dl_freq")
