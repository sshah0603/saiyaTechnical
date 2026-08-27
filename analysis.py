"""Analysis layer.

Parts map to functions as follows:
    Part 2  ->  frequency_table
    Part 3  ->  responder_comparison, compare_responders
    Part 4  ->  baseline_cohort, baseline_breakdown
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

DB_PATH = Path(__file__).resolve().parent / "cell_counts.db"


def connect(db_path: str | Path = DB_PATH) -> sqlite3.Connection:
    db_path = Path(db_path)
    if not db_path.exists():
        raise FileNotFoundError(
            f"{db_path.name} not found. Run `python load_data.py` first."
        )
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


# Part 2

FREQUENCY_SQL = """
SELECT sample, total_count, population, count, percentage
FROM sample_population_frequency
ORDER BY sample, population
"""


def frequency_table(conn: sqlite3.Connection) -> pd.DataFrame:
    """One row per (sample, population) with relative frequency as a percent."""
    return pd.read_sql(FREQUENCY_SQL, conn)


ANNOTATED_SQL = """
SELECT
    f.sample,
    f.total_count,
    f.population,
    f.count,
    f.percentage,
    pr.project_code               AS project,
    sb.subject_code               AS subject,
    sb.condition,
    sb.age,
    sb.sex,
    sm.treatment,
    sm.response,
    sm.sample_type,
    sm.time_from_treatment_start
FROM sample_population_frequency f
JOIN sample  sm ON sm.sample_code = f.sample
JOIN subject sb ON sb.subject_id  = sm.subject_id
JOIN project pr ON pr.project_id  = sb.project_id
"""


def annotated_frequencies(conn: sqlite3.Connection) -> pd.DataFrame:
    """Frequency table joined to sample and subject metadata.

    This is the working frame for the dashboard. Keeping it separate from
    frequency_table means the Part 2 deliverable stays exactly as specified.
    """
    return pd.read_sql(ANNOTATED_SQL, conn)


# Part 3

RESPONDER_SQL = ANNOTATED_SQL + """
WHERE sb.condition = :condition
  AND sm.treatment = :treatment
  AND sm.sample_type = :sample_type
  AND sm.response IN ('yes', 'no')
ORDER BY f.sample, f.population
"""


def responder_comparison(
    conn: sqlite3.Connection,
    condition: str = "melanoma",
    treatment: str = "miraclib",
    sample_type: str = "PBMC",
) -> pd.DataFrame:
    """Melanoma miraclib PBMC samples with a recorded response, by default."""
    return pd.read_sql(
        RESPONDER_SQL,
        conn,
        params={
            "condition": condition,
            "treatment": treatment,
            "sample_type": sample_type,
        },
    )


def _cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta: P(a > b) - P(a < b). Range -1 to 1, 0 means no shift."""
    if len(a) == 0 or len(b) == 0:
        return np.nan
    diff = np.sign(a[:, None] - b[None, :])
    return float(diff.mean())


def _bh_adjust(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR adjustment, implemented locally to avoid a dep."""
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    if ok.sum() == 0:
        return out.tolist()
    sub = p[ok]
    n = sub.size
    order = np.argsort(sub)
    ranked = sub[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    restored = np.empty(n)
    restored[order] = adj
    out[ok] = restored
    return out.tolist()


def compare_responders(
    df: pd.DataFrame,
    aggregate_by_subject: bool = True,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Test responder vs non-responder frequency for each population.

    aggregate_by_subject collapses a subject's repeat samples to their mean
    before testing. Subjects contribute multiple samples across timepoints
    and replicates, so treating each sample as independent inflates n and
    understates the p-values. Set it False to test at the sample level.

    Primary test is Mann-Whitney U, two sided. It makes no normality
    assumption, which suits frequency data at these group sizes. Welch's
    t-test is reported alongside as a cross-check; when the two disagree,
    the result is fragile and should be read as such. p-values are adjusted
    across the five populations using Benjamini-Hochberg.
    """
    unit = "subject" if aggregate_by_subject else "sample"

    if aggregate_by_subject:
        work = (
            df.groupby(["subject", "response", "population"], as_index=False)["percentage"]
            .mean()
        )
    else:
        work = df[["sample", "response", "population", "percentage"]].copy()

    rows = []
    for population, grp in work.groupby("population", sort=False):
        resp = grp.loc[grp["response"] == "yes", "percentage"].to_numpy()
        nonresp = grp.loc[grp["response"] == "no", "percentage"].to_numpy()

        if len(resp) < 3 or len(nonresp) < 3:
            rows.append(
                {
                    "population": population,
                    "n_responder": len(resp),
                    "n_non_responder": len(nonresp),
                    "median_responder": np.median(resp) if len(resp) else np.nan,
                    "median_non_responder": np.median(nonresp) if len(nonresp) else np.nan,
                    "median_difference": np.nan,
                    "mannwhitney_u": np.nan,
                    "p_value": np.nan,
                    "welch_p_value": np.nan,
                    "cliffs_delta": np.nan,
                    "note": "too few observations to test",
                }
            )
            continue

        u_stat, p_mw = stats.mannwhitneyu(resp, nonresp, alternative="two-sided")
        _, p_welch = stats.ttest_ind(resp, nonresp, equal_var=False)

        rows.append(
            {
                "population": population,
                "n_responder": len(resp),
                "n_non_responder": len(nonresp),
                "median_responder": float(np.median(resp)),
                "median_non_responder": float(np.median(nonresp)),
                "median_difference": float(np.median(resp) - np.median(nonresp)),
                "mannwhitney_u": float(u_stat),
                "p_value": float(p_mw),
                "welch_p_value": float(p_welch),
                "cliffs_delta": _cliffs_delta(resp, nonresp),
                "note": "",
            }
        )

    out = pd.DataFrame(rows)
    out["p_value_adj"] = _bh_adjust(out["p_value"].tolist())
    out["significant"] = out["p_value_adj"] < alpha
    out["effect_size_label"] = out["cliffs_delta"].abs().map(_delta_label)
    out.attrs["unit_of_analysis"] = unit
    out.attrs["alpha"] = alpha

    cols = [
        "population",
        "n_responder",
        "n_non_responder",
        "median_responder",
        "median_non_responder",
        "median_difference",
        "mannwhitney_u",
        "p_value",
        "p_value_adj",
        "significant",
        "welch_p_value",
        "cliffs_delta",
        "effect_size_label",
        "note",
    ]
    return out[cols].sort_values("p_value_adj", na_position="last").reset_index(drop=True)


def _delta_label(value: float) -> str:
    if pd.isna(value):
        return ""
    v = abs(value)
    if v < 0.147:
        return "negligible"
    if v < 0.33:
        return "small"
    if v < 0.474:
        return "medium"
    return "large"


# Part 4

BASELINE_SQL = """
SELECT
    sm.sample_code                AS sample,
    pr.project_code               AS project,
    sb.subject_code               AS subject,
    sb.condition,
    sb.age,
    sb.sex,
    sm.treatment,
    sm.response,
    sm.sample_type,
    sm.time_from_treatment_start,
    t.total_count
FROM sample sm
JOIN subject       sb ON sb.subject_id = sm.subject_id
JOIN project       pr ON pr.project_id = sb.project_id
JOIN sample_totals t  ON t.sample_id   = sm.sample_id
WHERE sb.condition = :condition
  AND sm.treatment = :treatment
  AND sm.sample_type = :sample_type
  AND sm.time_from_treatment_start = :timepoint
ORDER BY pr.project_code, sb.subject_code, sm.sample_code
"""


def baseline_cohort(
    conn: sqlite3.Connection,
    condition: str = "melanoma",
    treatment: str = "miraclib",
    sample_type: str = "PBMC",
    timepoint: int = 0,
) -> pd.DataFrame:
    """Melanoma PBMC samples at baseline from miraclib-treated patients."""
    return pd.read_sql(
        BASELINE_SQL,
        conn,
        params={
            "condition": condition,
            "treatment": treatment,
            "sample_type": sample_type,
            "timepoint": timepoint,
        },
    )


def baseline_breakdown(cohort: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Counts by project, response and sex for a baseline cohort.

    The brief asks for samples per project and subjects per response and
    sex. Those are different denominators, since a subject can contribute
    more than one sample at baseline. Both are reported for each cut so the
    difference is visible rather than buried.
    """
    if cohort.empty:
        empty = pd.DataFrame()
        return {"by_project": empty, "by_response": empty, "by_sex": empty}

    subjects = cohort.drop_duplicates(subset="subject")

    def cut(column: str, label: str) -> pd.DataFrame:
        samples = (
            cohort.groupby(column, dropna=False)
            .size()
            .rename("samples")
            .reset_index()
        )
        subs = (
            subjects.groupby(column, dropna=False)
            .size()
            .rename("subjects")
            .reset_index()
        )
        out = samples.merge(subs, on=column, how="outer").fillna(0)
        out[["samples", "subjects"]] = out[["samples", "subjects"]].astype(int)
        return out.rename(columns={column: label}).sort_values(label).reset_index(drop=True)

    return {
        "by_project": cut("project", "project"),
        "by_response": cut("response", "response"),
        "by_sex": cut("sex", "sex"),
    }


# ---------------------------------------------------------------- CLI

def _print(title: str, df: pd.DataFrame) -> None:
    print(f"\n{title}")
    print("-" * len(title))
    if df.empty:
        print("(no rows)")
    else:
        print(df.to_string(index=False))


def main() -> None:
    conn = connect()
    try:
        freq = frequency_table(conn)
        _print(f"Part 2: relative frequency table ({len(freq)} rows, first 10 shown)",
               freq.head(10))

        cohort_df = responder_comparison(conn)
        results = compare_responders(cohort_df)
        _print("Part 3: responders vs non-responders, melanoma miraclib PBMC", results)

        base = baseline_cohort(conn)
        _print(f"Part 4: baseline cohort ({len(base)} samples, first 10 shown)",
               base.head(10))
        for name, table in baseline_breakdown(base).items():
            _print(f"Part 4: {name}", table)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
