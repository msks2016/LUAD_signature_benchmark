# Project Handover

## What the project is
A systematic benchmark of 15 published LUAD lncRNA prognostic signatures
across one reference cohort (TCGA-LUAD) and four independent cohorts
(CPTAC-3, GSE31210, GSE50081, GSE30219).

## Current state
- Full pipeline exists and runs end-to-end.
- Manuscript is drafted, revised per a detailed review, and near submission.
- One known open issue: TCGA-LUAD log-transform correction was applied on
  2026-10-13 and the manuscript must be updated to match the corrected
  pipeline (see "Outstanding work").

## Where to start reading
1. README.md — how to install and run the pipeline
2. docs/deviations_log.md — every correction made and why
3. results/tables/ — all derived tables
4. figures/ — all figures
5. Manuscript draft in the repo root

## What was corrected most recently (2026-10-13)
1. TCGA-LUAD expression matrix was on raw TPM scale (0 to 136,166) while the
   four independent cohorts were on log scale (0 to ~17). The TCGA-LUAD v3
   matrix was log2-transformed and re-run. See `scripts/R/06b_prepare_TCGA_v3.R`.
2. `scripts/R/07b_score_signatures_relaxed.R` had three bugs:
   - referenced `surv$event` / `surv$time_months` (columns do not exist;
     they are `OS` / `OS_time`)
   - applied `log2(expr + 1)` a second time after the v3 matrix was already
     log-transformed
   - used a non-deterministic rule to pick one primary-tumor aliquot per
     patient
   All three were fixed.
3. `scripts/R/09_stage_lrt.R` had the same survival column bug and the same
   double-log bug. The survival column bug is fixed. The double-log bug is
   fixed. A remaining issue with the C-index scale is unresolved (see
   "Outstanding work").

## Outstanding work
1. `scripts/R/09_stage_lrt.R` writes C-index values in the thousands
   (24,849 and 33,000–35,000) instead of the 0–1 C-index scale. The value
   looks like a concordance count rather than a concordance fraction. The
   line that computes the C-index must be checked: `concordance(...)$count`
   should be `concordance(...)$concordance`. Regenerate
   `results/tables/stage_adjusted_LRT.csv` after the fix, then update the
   four ΔC-index values in the manuscript.
2. `scripts/22_extra_tables.R` reports a stale leave-one-out pooled HR range
   (1.68–1.8). The corrected four-cohort range is 1.380–1.496. Check what
   file the script reads for this value.
3. Manuscript patient count needs to be updated from 505 to 504 to match
   the corrected pipeline.
4. Figure 3 TCGA-LUAD reference HR needs updating from 1.635 (1.419–1.883)
   to 1.801 (1.560–2.079).
5. Figure 2, Supplementary Table S9, Supplementary Table S10 need
   regeneration with the corrected TCGA-LUAD metrics.
6. Supplementary Table S28 (Yao gene × cohort crosswalk) revealed that the
   Yao 2021 signature has fewer genes present per cohort than the manuscript
   claims. This requires a decision: verify the coverage accounting against
   the actual pipeline, or correct the manuscript.

## Known open reviewer issues
The manuscript has been reviewed and nine issues were raised. Issues 1–7 are
fully addressed. Issue 8 is partially addressed. Issue 9 has a critical open
finding (see #6 above). Full details in the review file and in
docs/deviations_log.md.

## Reproducing the pipeline
cd "/Users/apple/Desktop/Lab Works/BenchMark3/BenchMark3"
R
renv::restore()
source("scripts/bash/run_all.sh")
