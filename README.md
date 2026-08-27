# Cell population analysis

Loads `cell-count.csv` into SQLite, computes relative frequencies for five
immune populations, compares responders against non-responders for miraclib
in melanoma, and serves the results as an interactive dashboard.

## Run it

```bash
pip install -r requirements.txt
python load_data.py        # creates cell_counts.db in the repo root
streamlit run app.py       # dashboard at http://localhost:8501
```

Two extras:

```bash
python analysis.py         # prints Parts 2, 3 and 4 to the terminal
python make_report.py      # writes report.html, a static snapshot to email around
python -m pytest -q        # 18 tests
```

## Schema

Five tables and two views. `schema.sql` carries the full rationale in comments.

```
project ──< subject ──< sample ──< cell_count >── population
```

The source CSV is one wide row per sample with the populations as columns.
That layout makes a sixth population a schema migration and makes "frequency
of each population" an awkward unpivot. The schema normalizes to one row per
`(sample, population)`, so a new population is inserted as data. The composite
primary key on `cell_count` prevents two counts for the same population in the
same sample.

Attribute placement follows what changes and when:

- **subject** holds demographics and disease condition. These do not change
  between draws for the same subject.
- **sample** holds treatment, response, sample type and timepoint. Response is
  logically a subject-per-treatment fact. Holding it on the sample tolerates a
  subject who crosses to a second treatment arm without a schema change. The
  tradeoff is that the value repeats across a subject's samples, so the loader
  has to keep it consistent. `load_data.py` warns on conflicting values rather
  than silently picking one.

Two views do the arithmetic in the database:

- `sample_totals` sums counts per sample across whatever populations exist.
- `sample_population_frequency` is the Part 2 answer, one row per population
  per sample with the percentage already computed.

Putting the percentage in a view means every consumer, the dashboard, the CLI,
the static report, reads the same definition. Nothing recomputes it locally.

## Scaling

The current file is small enough that SQLite and pandas are the right tools.
Two things change as this grows toward hundreds of projects and thousands of
samples:

- **Storage.** The normalized layout is already the shape you want. Moving to
  Postgres is a connection-string change plus swapping the two views over. The
  natural next step is a partial index on `sample(sample_type, treatment,
  time_from_treatment_start)` and, if analytics dominate, a materialized
  version of `sample_population_frequency` refreshed on load.
- **Analytics.** The frequency view scales fine. The statistics do not scale in
  the sense that matters. Adding samples does not fix the underlying problem
  that subjects contribute multiple correlated samples. See below.

Adding a new analysis is a function in `analysis.py` returning a DataFrame plus
a figure builder in `figures.py`. Both the Streamlit app and the static report
pick it up.

## What the statistics do and do not say

The responder comparison uses a two sided Mann-Whitney U test on each of the
five populations. It ranks the two groups rather than assuming normality, which
suits frequency data at these group sizes. A Welch t-test runs alongside as a
cross-check. Agreement between the two is reassuring. Disagreement means the
result rests on a handful of observations and should be reported as
preliminary.

Three deliberate choices worth flagging in review:

- **Multiple testing.** Five populations tested at 0.05 gives roughly a one in
  four chance of at least one false positive. p-values are adjusted with
  Benjamini-Hochberg and `p_value_adj` is the column to quote.
- **Unit of analysis.** Subjects contribute several samples across timepoints
  and replicates. Treating each sample as independent inflates n and overstates
  significance. The default collapses a subject's samples to their mean before
  testing. The dashboard exposes a toggle so the difference is visible rather
  than buried, and the sample-level view is there for inspection, not for
  quoting.
- **Effect size.** Cliff's delta accompanies every p-value. A small p-value
  paired with a negligible delta is a real but clinically thin difference, and
  the table says so in words.

Two caveats that no amount of code fixes. The five populations are
compositional, so they sum to 100 percent by construction and one going up
forces others down. A shift in monocytes and a shift in CD4 T cells are not
independent findings. And this is an observational comparison within a treated
arm, so a population that separates responders from non-responders is a
candidate biomarker rather than a mechanism. Confirming either point needs a
prespecified analysis on a held-out cohort.

## Files

| File | Purpose |
| --- | --- |
| `load_data.py` | Builds the schema and loads the CSV. No arguments. |
| `schema.sql` | Tables, indexes, views, with design rationale in comments. |
| `analysis.py` | Parts 2, 3 and 4 as functions returning DataFrames. Also a CLI. |
| `figures.py` | Design tokens and Plotly figure builders, shared by the app and the report. |
| `app.py` | Streamlit dashboard. |
| `make_report.py` | Static HTML export of the same analysis. |
| `test_pipeline.py` | Tests for the load, the arithmetic, the statistics and the filters. |
| `make_sample_csv.py` | Generates a stand-in `cell-count.csv`. Delete once the real export is in place. |

## Note on the data file

`cell-count.csv` in this repo is generated by `make_sample_csv.py`, not the real
export. It matches the expected column layout and carries a planted responder
signal so the statistical path is exercised end to end. Replace the file with
the production export and every script runs unchanged.
