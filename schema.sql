-- Schema for the Loblaw Bio cell-count database.
--
-- Design notes
--
-- The source CSV is one wide row per sample, with the five populations as
-- columns. That layout makes adding a sixth population a schema migration
-- and makes "frequency of each population" an awkward unpivot. The schema
-- below normalizes to one row per (sample, population), so new populations
-- are inserted as data.
--
-- Attribute placement:
--   * project        -> its own table. Projects gain samples over time.
--   * subject        -> demographics and disease condition live here. These
--                       do not change between draws for the same subject.
--   * sample         -> treatment, response, sample_type and timepoint live
--                       here. Response is logically a subject-per-treatment
--                       fact, but holding it on the sample tolerates
--                       subjects who cross over to a second treatment arm
--                       without a schema change. The cost is that the value
--                       repeats across a subject's samples, so loaders must
--                       keep it consistent.
--   * cell_count     -> the measurement. Composite PK prevents a duplicate
--                       count for the same population in the same sample.

PRAGMA foreign_keys = ON;

DROP VIEW  IF EXISTS sample_population_frequency;
DROP VIEW  IF EXISTS sample_totals;
DROP TABLE IF EXISTS cell_count;
DROP TABLE IF EXISTS sample;
DROP TABLE IF EXISTS subject;
DROP TABLE IF EXISTS project;
DROP TABLE IF EXISTS population;

CREATE TABLE project (
    project_id   INTEGER PRIMARY KEY,
    project_code TEXT NOT NULL UNIQUE
);

CREATE TABLE subject (
    subject_id   INTEGER PRIMARY KEY,
    subject_code TEXT NOT NULL UNIQUE,
    project_id   INTEGER NOT NULL REFERENCES project(project_id),
    condition    TEXT,
    age          INTEGER CHECK (age IS NULL OR age BETWEEN 0 AND 120),
    sex          TEXT CHECK (sex IS NULL OR sex IN ('M', 'F'))
);

CREATE TABLE population (
    population_id INTEGER PRIMARY KEY,
    population    TEXT NOT NULL UNIQUE,
    display_name  TEXT
);

CREATE TABLE sample (
    sample_id                 INTEGER PRIMARY KEY,
    sample_code               TEXT NOT NULL UNIQUE,
    subject_id                INTEGER NOT NULL REFERENCES subject(subject_id),
    treatment                 TEXT,
    response                  TEXT CHECK (response IS NULL OR response IN ('yes', 'no')),
    sample_type               TEXT,
    time_from_treatment_start INTEGER
);

CREATE TABLE cell_count (
    sample_id     INTEGER NOT NULL REFERENCES sample(sample_id) ON DELETE CASCADE,
    population_id INTEGER NOT NULL REFERENCES population(population_id),
    count         INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id)
);

CREATE INDEX idx_subject_project      ON subject(project_id);
CREATE INDEX idx_sample_subject       ON sample(subject_id);
CREATE INDEX idx_sample_filter        ON sample(sample_type, treatment, time_from_treatment_start);
CREATE INDEX idx_cell_count_pop       ON cell_count(population_id);

-- Total cells per sample, summed across every population present.
CREATE VIEW sample_totals AS
SELECT
    s.sample_id,
    s.sample_code,
    SUM(cc.count) AS total_count
FROM sample s
JOIN cell_count cc ON cc.sample_id = s.sample_id
GROUP BY s.sample_id, s.sample_code;

-- Part 2 answer, expressed in the database rather than in pandas.
CREATE VIEW sample_population_frequency AS
SELECT
    s.sample_code                              AS sample,
    t.total_count                              AS total_count,
    p.population                               AS population,
    cc.count                                   AS count,
    ROUND(100.0 * cc.count / t.total_count, 4) AS percentage
FROM cell_count cc
JOIN sample        s ON s.sample_id     = cc.sample_id
JOIN population    p ON p.population_id = cc.population_id
JOIN sample_totals t ON t.sample_id     = cc.sample_id;
