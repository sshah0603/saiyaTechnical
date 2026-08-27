#!/usr/bin/env python3
"""Initialize the SQLite database and load cell-count.csv into it.

Run from the repository root:

    python load_data.py

Creates cell_counts.db in the repository root. Safe to re-run: the schema
is dropped and rebuilt each time, so the load is idempotent.
"""

import csv
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
SCHEMA_PATH = ROOT / "schema.sql"

# Columns that are population measurements. Everything else is metadata.
# Adding a population here is the only change needed to ingest a sixth column.
POPULATION_COLUMNS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

DISPLAY_NAMES = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8+ T cell",
    "cd4_t_cell": "CD4+ T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}


def clean(value):
    """Normalize a CSV field: strip whitespace, treat blanks as NULL."""
    if value is None:
        return None
    value = value.strip()
    return value if value else None


def to_int(value):
    value = clean(value)
    if value is None:
        return None
    return int(float(value))


def read_rows(csv_path):
    with open(csv_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in POPULATION_COLUMNS if c not in reader.fieldnames]
        if missing:
            raise SystemExit(
                f"cell-count.csv is missing expected population columns: {missing}"
            )
        return list(reader)


def create_schema(conn):
    conn.executescript(SCHEMA_PATH.read_text())


def load(conn, rows):
    cur = conn.cursor()

    # Populations first. Order here fixes the display order downstream.
    for name in POPULATION_COLUMNS:
        cur.execute(
            "INSERT INTO population (population, display_name) VALUES (?, ?)",
            (name, DISPLAY_NAMES.get(name, name)),
        )
    pop_ids = dict(cur.execute("SELECT population, population_id FROM population"))

    project_ids = {}
    subject_ids = {}
    subject_seen = {}
    warnings = []

    for line_no, row in enumerate(rows, start=2):
        project_code = clean(row.get("project"))
        subject_code = clean(row.get("subject"))
        sample_code = clean(row.get("sample"))

        if not sample_code:
            warnings.append(f"line {line_no}: no sample id, row skipped")
            continue
        if not subject_code:
            warnings.append(f"line {line_no}: no subject id, row skipped")
            continue

        if project_code not in project_ids:
            cur.execute(
                "INSERT OR IGNORE INTO project (project_code) VALUES (?)",
                (project_code,),
            )
            project_ids[project_code] = cur.execute(
                "SELECT project_id FROM project WHERE project_code = ?",
                (project_code,),
            ).fetchone()[0]

        subject_attrs = (
            project_ids[project_code],
            clean(row.get("condition")),
            to_int(row.get("age")),
            clean(row.get("sex")),
        )
        if subject_code not in subject_ids:
            cur.execute(
                """INSERT INTO subject (subject_code, project_id, condition, age, sex)
                   VALUES (?, ?, ?, ?, ?)""",
                (subject_code, *subject_attrs),
            )
            subject_ids[subject_code] = cur.lastrowid
            subject_seen[subject_code] = subject_attrs
        elif subject_seen[subject_code] != subject_attrs:
            warnings.append(
                f"line {line_no}: subject {subject_code} has conflicting "
                f"demographics across rows; first values kept"
            )

        try:
            cur.execute(
                """INSERT INTO sample
                       (sample_code, subject_id, treatment, response,
                        sample_type, time_from_treatment_start)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    sample_code,
                    subject_ids[subject_code],
                    clean(row.get("treatment")),
                    clean(row.get("response")),
                    clean(row.get("sample_type")),
                    to_int(row.get("time_from_treatment_start")),
                ),
            )
        except sqlite3.IntegrityError:
            warnings.append(f"line {line_no}: duplicate sample id {sample_code}, skipped")
            continue

        sample_id = cur.lastrowid
        for pop in POPULATION_COLUMNS:
            count = to_int(row.get(pop))
            if count is None:
                warnings.append(
                    f"line {line_no}: sample {sample_code} has no {pop} count, stored as 0"
                )
                count = 0
            cur.execute(
                "INSERT INTO cell_count (sample_id, population_id, count) VALUES (?, ?, ?)",
                (sample_id, pop_ids[pop], count),
            )

    conn.commit()
    return warnings


def summarize(conn):
    q = lambda sql: conn.execute(sql).fetchone()[0]
    return {
        "projects": q("SELECT COUNT(*) FROM project"),
        "subjects": q("SELECT COUNT(*) FROM subject"),
        "samples": q("SELECT COUNT(*) FROM sample"),
        "populations": q("SELECT COUNT(*) FROM population"),
        "cell_counts": q("SELECT COUNT(*) FROM cell_count"),
    }


def main():
    if not CSV_PATH.exists():
        raise SystemExit(f"Expected {CSV_PATH.name} in the repository root. Not found.")
    if not SCHEMA_PATH.exists():
        raise SystemExit(f"Expected {SCHEMA_PATH.name} in the repository root. Not found.")

    if DB_PATH.exists():
        DB_PATH.unlink()

    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        create_schema(conn)
        rows = read_rows(CSV_PATH)
        warnings = load(conn, rows)
        counts = summarize(conn)
    finally:
        conn.close()

    print(f"Created {DB_PATH.name}")
    for label, value in counts.items():
        print(f"  {label:<13} {value}")

    if warnings:
        print(f"\n{len(warnings)} data warning(s):")
        for w in warnings[:20]:
            print(f"  {w}")
        if len(warnings) > 20:
            print(f"  ... and {len(warnings) - 20} more")

    return 0


if __name__ == "__main__":
    sys.exit(main())
