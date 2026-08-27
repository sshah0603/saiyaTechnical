"""Tests for the cell-count pipeline.

    python -m pytest test_pipeline.py -q

Covers the load, the Part 2 arithmetic, the Part 3 statistics and the
Part 4 filters. 
"""

from __future__ import annotations

import sqlite3
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import analysis as an

ROOT = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def conn():
    if not (ROOT / "cell_counts.db").exists():
        subprocess.run([sys.executable, "load_data.py"], cwd=ROOT, check=True)
    c = an.connect()
    yield c
    c.close()


# ---------------------------------------------------------------- load

def test_load_is_idempotent():
    for _ in range(2):
        subprocess.run([sys.executable, "load_data.py"], cwd=ROOT, check=True)
    c = sqlite3.connect(ROOT / "cell_counts.db")
    try:
        dupes = c.execute(
            "SELECT COUNT(*) FROM (SELECT sample_code FROM sample GROUP BY sample_code HAVING COUNT(*) > 1)"
        ).fetchone()[0]
    finally:
        c.close()
    assert dupes == 0


def test_every_sample_has_every_population(conn):
    n_pops = conn.execute("SELECT COUNT(*) FROM population").fetchone()[0]
    bad = conn.execute(
        "SELECT COUNT(*) FROM (SELECT sample_id FROM cell_count GROUP BY sample_id HAVING COUNT(*) != ?)",
        (n_pops,),
    ).fetchone()[0]
    assert bad == 0


def test_row_count_matches_csv(conn):
    csv_rows = len(pd.read_csv(ROOT / "cell-count.csv"))
    db_rows = conn.execute("SELECT COUNT(*) FROM sample").fetchone()[0]
    assert db_rows == csv_rows


def test_foreign_keys_hold(conn):
    assert conn.execute("PRAGMA foreign_key_check").fetchall() == []


# ---------------------------------------------------------------- Part 2

def test_frequency_table_shape(conn):
    freq = an.frequency_table(conn)
    assert list(freq.columns) == ["sample", "total_count", "population", "count", "percentage"]
    n_samples = conn.execute("SELECT COUNT(*) FROM sample").fetchone()[0]
    n_pops = conn.execute("SELECT COUNT(*) FROM population").fetchone()[0]
    assert len(freq) == n_samples * n_pops


def test_percentages_sum_to_100(conn):
    freq = an.frequency_table(conn)
    sums = freq.groupby("sample")["percentage"].sum()
    assert np.allclose(sums, 100.0, atol=0.01)


def test_total_count_equals_sum_of_counts(conn):
    freq = an.frequency_table(conn)
    per_sample = freq.groupby("sample").agg(
        stated=("total_count", "first"), summed=("count", "sum")
    )
    assert (per_sample["stated"] == per_sample["summed"]).all()


def test_percentage_matches_count_over_total(conn):
    freq = an.frequency_table(conn)
    expected = 100.0 * freq["count"] / freq["total_count"]
    assert np.allclose(freq["percentage"], expected, atol=0.001)


# ---------------------------------------------------------------- Part 3

def test_responder_filter_is_tight(conn):
    df = an.responder_comparison(conn)
    assert set(df["condition"]) == {"melanoma"}
    assert set(df["treatment"]) == {"miraclib"}
    assert set(df["sample_type"]) == {"PBMC"}
    assert set(df["response"]) <= {"yes", "no"}


def test_subject_aggregation_reduces_n(conn):
    df = an.responder_comparison(conn)
    per_sample = an.compare_responders(df, aggregate_by_subject=False)
    per_subject = an.compare_responders(df, aggregate_by_subject=True)
    assert per_subject["n_responder"].max() <= per_sample["n_responder"].max()


def test_adjusted_p_is_never_below_raw(conn):
    res = an.compare_responders(an.responder_comparison(conn))
    ok = res["p_value"].notna()
    assert (res.loc[ok, "p_value_adj"] >= res.loc[ok, "p_value"] - 1e-12).all()


def test_no_difference_yields_no_significance():
    rng = np.random.default_rng(0)
    rows = []
    for i in range(40):
        for pop in ["b_cell", "cd4_t_cell", "monocyte"]:
            rows.append(
                {
                    "subject": f"s{i}",
                    "sample": f"x{i}",
                    "response": "yes" if i % 2 else "no",
                    "population": pop,
                    "percentage": rng.normal(20, 3),
                }
            )
    res = an.compare_responders(pd.DataFrame(rows))
    assert not res["significant"].any()


def test_planted_difference_is_detected():
    rows = []
    for i in range(24):
        responder = i % 2 == 0
        rows.append(
            {
                "subject": f"s{i}",
                "sample": f"x{i}",
                "response": "yes" if responder else "no",
                "population": "cd4_t_cell",
                "percentage": 42.0 + i * 0.05 if responder else 30.0 + i * 0.05,
            }
        )
        rows.append(
            {
                "subject": f"s{i}",
                "sample": f"x{i}",
                "response": "yes" if responder else "no",
                "population": "nk_cell",
                "percentage": 9.0 + (i % 5) * 0.1,
            }
        )
    res = an.compare_responders(pd.DataFrame(rows)).set_index("population")
    assert res.loc["cd4_t_cell", "significant"]
    assert not res.loc["nk_cell", "significant"]


def test_cliffs_delta_bounds():
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([10.0, 11.0, 12.0])
    assert an._cliffs_delta(a, b) == pytest.approx(-1.0)
    assert an._cliffs_delta(b, a) == pytest.approx(1.0)
    assert an._cliffs_delta(a, a) == pytest.approx(0.0)


def test_small_group_is_flagged_not_tested():
    rows = [
        {"subject": f"s{i}", "sample": f"x{i}", "response": "yes" if i < 2 else "no",
         "population": "b_cell", "percentage": 10.0 + i}
        for i in range(8)
    ]
    res = an.compare_responders(pd.DataFrame(rows))
    assert res.loc[0, "note"] == "too few observations to test"
    assert pd.isna(res.loc[0, "p_value"])


# ---------------------------------------------------------------- Part 4

def test_baseline_filter(conn):
    base = an.baseline_cohort(conn)
    assert set(base["condition"]) == {"melanoma"}
    assert set(base["treatment"]) == {"miraclib"}
    assert set(base["sample_type"]) == {"PBMC"}
    assert set(base["time_from_treatment_start"]) == {0}


def test_breakdown_totals_reconcile(conn):
    base = an.baseline_cohort(conn)
    cuts = an.baseline_breakdown(base)
    for cut in cuts.values():
        assert cut["samples"].sum() == len(base)
        assert cut["subjects"].sum() == base["subject"].nunique()


def test_breakdown_handles_empty():
    cuts = an.baseline_breakdown(pd.DataFrame())
    assert all(df.empty for df in cuts.values())
