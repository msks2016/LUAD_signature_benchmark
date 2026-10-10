# OSF Update — 2026-10-10 — Coding error correction

## What happened
The pipeline used in the initial transportability and direction-consistency
analyses applied survival::concordance() without reverse = TRUE for a risk
score in which higher values indicate higher mortality risk. This inverted
the C-index for affected signatures. The prior 0.323 value for Yao 2021
in CPTAC-3 was 1 - 0.677; the true value is 0.677.

## Consequences
- The earlier direction-inversion finding does not survive correction
- The Yao 2021 post-hoc disclosure is no longer relevant as a reversal
- The primary hypotheses (H1, H2, H3) remain intact and are re-evaluated
  in the corrected pipeline

## Actions taken
- reverse = TRUE added to every formula-based concordance() call
- All downstream analyses re-run
- Earlier results preserved in the repository's git history and marked
  superseded in the README
- A deviations log records the correction with date and reason

Nothing has been deleted. The correction is documented for transparency.
