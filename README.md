# Cell population analysis

Loads `cell-count.csv` into SQLite, computes the relative frequency of five
immune populations in every sample, tests which populations separate miraclib
responders from non-responders in melanoma, and serves the result as an
interactive dashboard.


Both open in a browser with nothing installed. 

---

## Running it

Three targets, as specified.

```bash
make setup       # install dependencies from requirements.txt
make pipeline    # build the database, generate every table and plot
make dashboard   # start the dashboard on port 8501
```

`make setup` then `make pipeline` runs start to finish with no manual step.
`make pipeline` prints each stage and the file it wrote.

Two extra targets, not required but useful:

```bash
make test        # 26 tests
make clean       # remove the database and generated outputs
```


### Outputs

`make pipeline` writes everything into `outputs/`, which is committed so the
results are readable without running anything.

| File | Contents |
| --- | --- |
| `part2_cell_frequencies.csv` | The Part 2 summary table. One row per population per sample. |
| `part2_population_summary.csv` | Mean, sd and range for each population. |
| `part2_composition.png` | Stacked composition of every sample. |
| `part2_population_distributions.png` | Frequency spread per population. |
| `part3_responder_statistics.csv` | The significance table. |
| `part3_responder_frequencies.csv` | The rows the tests actually ran on. |
| `part3_responder_boxplot.png` | The Part 3 boxplot. |
| `part4_baseline_cohort.csv` | Baseline melanoma PBMC samples, miraclib. |
| `part4_counts_by_project.csv` | Samples and subjects per project. |
| `part4_counts_by_response.csv` | Responders and non-responders. |
| `part4_counts_by_sex.csv` | Males and females. |
| `summary.md` | The findings in prose, with every table inline. |
| `report.html` (repo root) | Static snapshot of the whole dashboard. |

---

## Database schema

Five tables and two views.

```
project ──< subject ──< sample ──< cell_count >── population
```

```sql
project     (project_id, project_code)
subject     (subject_id, subject_code, project_id, condition, age, sex)
sample      (sample_id, sample_code, subject_id, treatment, response,
             sample_type, time_from_treatment_start)
population  (population_id, population, display_name)
cell_count  (sample_id, population_id, count)   -- PK (sample_id, population_id)
```

### Why this shape

**The measurements are normalized out of the wide CSV.** The source file puts
the five populations in five columns. That makes adding a sixth population a
schema migration, and it makes the central question of this exercise, the
frequency of each population, an awkward unpivot in every query that needs it.
One row per `(sample, population)` turns a new population into an insert. The
composite primary key on `cell_count` makes a duplicate count for the same
population in the same sample impossible rather than merely unlikely.

The cost is honest: the table is five times longer than the CSV, and reading a
single sample's full profile needs a pivot. For an analytics workload where the
common question is "this population across many samples" rather than "all
populations for one sample", that trade goes the right way.

**Attributes sit at the level they actually vary.** Demographics and disease
condition belong to the subject and do not change between draws, so they live
on `subject` and are stored once. Sample type, timepoint, treatment and
response belong to the sample.

Response is the debatable one. It is logically a fact about a subject under a
treatment, not about an individual draw, so a stricter design would put it on a
`subject_treatment` table. I kept it on `sample` because it tolerates a subject
crossing into a second treatment arm without a schema change, which is a
realistic thing for a trial to do. The cost is that the value repeats across a
subject's samples and nothing in the schema forces those copies to agree, so
`load_data.py` warns when they conflict instead of silently keeping one. If
crossover never happens in this program, promoting response to its own table is
the better call and the loader change is small.

**The arithmetic lives in views**

- `sample_totals` sums counts per sample across whatever populations exist.
- `sample_population_frequency` is the Part 2 answer, with `percentage`
  computed in SQL.

Every consumer, the dashboard, the pipeline, the CLI, the static report, reads
the same definition. Nothing recomputes a percentage locally, so no two views
of the data can disagree about what a frequency is.

### Scaling to hundreds of projects and thousands of samples

At that size the shape of the schema does not need to change. Four things do.

**Move to Postgres.** SQLite has one writer and no real concurrency. The DDL is
close to portable already; the changes are identity columns in place of
`INTEGER PRIMARY KEY`, and a connection string instead of a file path.
`analysis.py` uses parameterized SQL throughout, so the queries carry over.

**Index for the filters that actually get used.** The current composite index on
`sample(sample_type, treatment, time_from_treatment_start)` covers the Part 3
and Part 4 predicates. At scale, add `cell_count(population_id, sample_id)` to
support "one population across everything", and consider partitioning
`cell_count` by project once it passes a few hundred million rows.

**Materialize the frequency view.** `sample_population_frequency` recomputes a
division per row per query. That is free today and will not be at a hundred
million rows. Make it a materialized view refreshed at load time, or a real
table written by the loader. The definition stays in one place either way,
which is the property worth protecting.

**Split the write path from the read path.** Loading is append-heavy and
analytics are scan-heavy. Once they contend, the normalized tables stay as the
system of record and a columnar copy, Parquet on object storage queried through
DuckDB or a warehouse table, serves the analytics. The normalized layout is
already the right source for that copy.

---

## Code structure

| File | Job |
| --- | --- |
| `load_data.py` | Part 1. Builds the schema and loads the CSV. Root, no arguments. |
| `schema.sql` | The DDL as readable, commented SQL. |
| `analysis.py` | Parts 2, 3 and 4 as functions returning DataFrames. Also a CLI. |
| `figures.py` | Design tokens and Plotly figure builders. |
| `plots.py` | Static matplotlib plots written to disk. |
| `run_pipeline.py` | Orchestrates the pipeline and writes `outputs/`. |
| `app.py` | The Streamlit dashboard. |
| `make_report.py` | Static HTML snapshot of the analysis. |
| `test_pipeline.py` | 26 tests. |

### Why it is arranged this way

**One layer owns the analysis.** Every question in the brief is a function in
`analysis.py` that takes a connection and returns a DataFrame. The dashboard,
the pipeline, the static report and the CLI all call the same functions. A
correction to the statistics changes one place and propagates everywhere. The
alternative, where the dashboard recomputes something the pipeline already
computed, is how two views of the same data start disagreeing.

**Filtering happens in SQL, statistics happen in Python.** The database is good
at predicates and aggregation and has no rank tests. Splitting on that line
keeps the filters reusable and the statistics testable.

**Figures are separated from the code that displays them.** `figures.py` holds
the color tokens and the Plotly builders. The dashboard and the HTML report
both import from it, so a population is the same color in both. `plots.py` does
the same job in matplotlib for the on-disk PNGs. Two plotting libraries is a
real cost, taken deliberately: Plotly static export needs a headless Chrome
download, which is a poor thing to put in a grader's path.

**The database is a build artifact, not source.** `cell_counts.db` is
gitignored and regenerated by `load_data.py`. `app.py` builds it on first run
if it is missing, so a hosted deployment works on a cold clone.

**The loader is self-contained.** The DDL is embedded in `load_data.py` rather
than read from `schema.sql` at runtime, so the script works when copied into an
empty directory with only the CSV beside it. `schema.sql` stays because
commented SQL is easier to read than a Python string, and a test fails if the
two diverge.

**Tests cover the failures that produce a wrong number rather than an
exception.** Percentages summing to 100, the stated total matching the sum of
counts, the adjusted p never falling below the raw p, a planted difference
being detected, and a synthetic null producing no significance. Six more lock
in the Part 1 file requirements, the three Makefile targets, and the exact
column names Part 2 specifies.

---

## Part 3: how significance is decided

Five steps in `compare_responders()`.

1. Filter to melanoma, miraclib, PBMC, response recorded as yes or no.
2. Average each subject's repeat samples, so a subject counts once.
3. Two-sided Mann-Whitney U per population, with a Welch t-test alongside as a
   cross-check.
4. Benjamini-Hochberg correction across the five tests.
5. `significant = p_value_adj < 0.05`.

Mann-Whitney is the primary test because it ranks rather than assuming
normality, which matters at these group sizes where normality cannot be checked
with any confidence. Five tests at 0.05 carry roughly a 23 percent chance of at
least one false positive, so the adjusted column is the one to quote. Cliff's
delta sits beside every p-value, because a small p paired with a negligible
delta is a real but clinically thin difference and the table should say so.

A population with fewer than three observations in either arm is skipped and
returns `NaN`, which can never satisfy the threshold. Untestable and
non-significant stay distinguishable.

The dashboard exposes a toggle between subject-level and sample-level analysis.
Sample-level roughly quintuples n by counting the same subjects repeatedly and
drives the p-values down by many orders of magnitude. It is there for
inspection. The subject-level figure is the one to report.


