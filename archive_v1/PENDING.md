# Pending Tasks

## Priority 1 — Fix pipeline bugs
- [ ] Fix C-index scale in scripts/R/09_stage_lrt.R (use $concordance not $count)
- [ ] Rerun 09_stage_lrt.R and save output
- [ ] Fix stale leave-one-out range in scripts/22_extra_tables.R
- [ ] Rerun 22_extra_tables.R and save output

## Priority 2 — Update manuscript to match corrected pipeline
- [ ] 505 -> 504 patients everywhere
- [ ] Figure 3 TCGA-LUAD HR 1.635 -> 1.801 (95% CI 1.560-2.079)
- [ ] Regenerate Figure 2 with log-transformed TCGA-LUAD C-indexes
- [ ] Regenerate Figure 3 with new TCGA-LUAD reference row
- [ ] Regenerate Supplementary Table S9
- [ ] Regenerate Supplementary Table S10
- [ ] Update the four ΔC-index values in the stage-adjusted section if the
      regenerated stage_adjusted_LRT.csv changes them

## Priority 3 — Resolve the Yao coverage finding
- [ ] Cross-check the manuscript's per-cohort Yao gene counts against
      Supplementary_Table_S28_Yao_gene_cohort_crosswalk.csv
- [ ] If the pipeline scored fewer genes than claimed, correct the
      manuscript or the pipeline

## Priority 4 — Additional reviewer items
- [ ] Add null-test specification to Methods (Issue 8)
- [ ] Add sample-selection rule to Methods (Issue 9)
- [ ] Add S26 (recheck crosswalk) and S27 (conflict quotations)
- [ ] Shorten Abstract to <=400 words
- [ ] Fix "404 identified, 404 screened" duplicate-count wording
- [ ] Verify TRIPOD and reporting checklist files are supplied

## Priority 5 — Final submission
- [ ] Regenerate all figures at 300 dpi
- [ ] Assemble submission package
- [ ] Cover letter
- [ ] Final commit and push
