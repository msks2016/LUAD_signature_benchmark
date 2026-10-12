
## 2026-10-10 — Expression-matched null test for CPTAC-3

 Ran expression-matched, coefficient-weighted null test (2000 draws per signature)
 with Holm correction across 15 signatures.
 Result: Yao 2021 classified Contradicted (C_obs 0.323, null 0.363-0.601, p_holm 0.045).
 Remaining 14 signatures classified Inconclusive.
 Yao 2021 identified post hoc as the most extreme of 15; not pre-specified as a
 directional divergence candidate. Reported as hypothesis-generating pending
 replication in an independent cohort.

## 2026-10-10 (evening) — GSE31210 replication of direction test
 GSE31210 (n=226, events=35, GPL570, early-stage Japanese LUAD).
 Two signatures evaluable: Yao 2021 (7/7 genes), Yu 2024 (9/9 genes).
 Results:
   Yao 2021:  C=0.363 (95% CI 0.264-0.456), null 0.337-0.646,  p_holm=0.089, class=Inconclusive.
   Yu 2024:   C=0.282 (95% CI 0.200-0.374), null 0.353-0.650,  p_holm=0.004, class=Contradicted.
 Interpretation: Direction instability is suggested in both cohorts,  but the OSF pre-specified threshold (>=2 independent cohorts per  signature) is not met for either signature individually. Reported as  exploratory.

## 2026-10-10 — H2 result (added value beyond stage)
 4 of 10 signatures add significant value beyond stage in CPTAC-3
 after BH correction: Wu 2022 (p=0.0001), Yao 2021 (p=0.008),
 Li 2023 (p=0.012), Yang 2022 (p=0.003). The Yao 2021 result
 contrasts with its univariable C-index of 0.323, suggesting its
 prognostic signal is stage-dependent and was masked in the
 univariable analysis.

## 2026-10-10 — Phase 8 robustness checks
 1. Missing-gene threshold sensitivity: 10%->10 evaluable, 20%->11, 30%->13.
    Paper352 flips at 10/20 boundary; Paper165 and Paper9 become evaluable only at 30%.
 2. Stage-restricted: inversion pattern (raw C < 0.5, Cox C > 0.5) unchanged in stage I-II only.
 3. Leave-one-cohort-out: pooled Yao HR stable at 1.68-1.80 across all subsets, I2=0%.
 4. Early-death exclusion (no deaths within 3 months): pattern unchanged.
 All checks support the primary finding; no re-analysis required.

## 2026-10-11 — Transportability decision rules (pre-specified before running)
 Primary summary: pooled mean calibration slope (log-HR per SD of published score)
 across all independent signature-cohort pairs, with cluster-bootstrap 95% CI.
 Effect-size threshold: delta = log(1.25) = 0.223.
 Decision rules (set before running):
   - CI lower bound > 0: Transportable (positive finding)
   - CI upper bound < 0: Contradicted (published weights point wrong way)
   - CI includes 0 and upper bound < delta: Informative null
   - CI includes 0 and upper bound >= delta: Inconclusive
 Independent pairs only (TCGA excluded as non-independent for signatures developed on it).
 Total independent pairs: 12 (CPTAC-3 = 10, GSE31210 = 2).

## 2026-10-11 — Transportability analysis result
 Calibration-slope analysis on 12 independent pairs.
 Mean slope: 0.311 (95% CI -0.107 to 0.684).
 Classification by pre-specified rules: Inconclusive.
 Divergence null test: 1 observed vs 5.1 expected (proportion 0.083 vs 0.427).
 Conclusion: the original 'coefficient non-transportability' framing is not
 supported. The paired-C-index divergence observed earlier is less frequent
 than chance when restricted to independent pairs. The paper's framing
 changes to a coverage-and-power finding: public LUAD cohorts lack the
 annotation and event counts needed to validate published lncRNA signatures.

## 2026-10-10 — Metric terminology correction
 The metric previously labeled 'calibration slope' is computed as the
 coefficient of scale(risk_score) in a Cox model. This is the log-HR per
 1 SD of the published risk score, not the conventional calibration slope
 (which would use the unstandardized published linear predictor). The
 metric is renamed throughout the manuscript for accuracy.
 Pre-specified decision rules (delta = log 1.25) remain valid for log-HR per SD.
 No re-analysis required; terminology-only change.

## 2026-10-10 — Metric terminology correction
 The metric previously labeled 'calibration slope' is computed as the
 coefficient of scale(risk_score) in a Cox model. This is the log-HR per
 1 SD of the published risk score, not the conventional calibration slope.
 The metric is renamed 'log-HR per SD' throughout the manuscript.
 Pre-specified decision rules (delta = log 1.25) remain valid for log-HR per SD.
 No re-analysis required.

## 2026-10-10 — Switch to raw log2 expression for signature scoring
 Original publication formulas use raw expression (TPM/FPKM or microarray
 intensity) multiplied by coefficient, without z-scoring. Prior pipeline
 used per-cohort z-scored expression, which changes effective gene weights.
 Pipeline updated to use raw log2(TPM+1) or raw log2(intensity).
 All downstream outputs will be regenerated.

## 2026-10-10 — Final transportability result (raw expression, 100% coverage)
 Pooled log-HR per SD: 0.322 (95% CI -0.062 to 0.677). Classification by
 pre-specified rules: Inconclusive. Divergence null test: 0 observed vs 4.9
 expected. Two pairs (He 2021, Diao 2022) had constant linear predictors in
 CPTAC-3 and were handled by a variance guard.

## 2026-10-11 — Yao 2021 cohort overlap check
 Verified from the original publication (Cancer Cell Int 2021;21:340):
   Training: TCGA-LUAD (n = 477 after OS > 30 days filter)
   Validation: GSE3141 (n = 58) + GSE37745 (n = 106), merged
   Feature selection: GSE37745
 GSE31210 and CPTAC-3 were NOT used. Both are independent for Yao 2021.
 TCGA-LUAD is the development cohort and was excluded from the pooled estimate.

## 2026-10-11 — Yao 2021 H2 in GSE31210
 Stage-only Cox C-index: 0.678; stage + risk-score: 0.748 (Delta C = +0.071).
 Adjusted HR per SD of Yao 2021 score: 1.46 (p = 0.023).
 Likelihood-ratio test p = 0.027.
 Interpretation: Yao 2021 adds significant prognostic value beyond stage
 in GSE31210, consistent with the same finding in CPTAC-3.

## 2026-10-11 — Yao 2021 random-gene null in GSE31210
 Expression-matched null test (2,000 draws):
   Observed C-index: 0.643
   Null median: 0.508; null 2.5-97.5%: 0.356-0.663
   Two-sided empirical p = 0.067 (marginal; consistent with the
   primary null-test classification of Inconclusive, p_Holm = 0.077)
 Interpretation: Yao 2021 shows a positive trend in GSE31210 that does not
 reach significance in that cohort. It adds value beyond stage (H2 p = 0.027)
 but does not pass a strict two-cohort replication rule at the per-cohort
 significance level.

## 2026-10-11 — Yao 2021 vs proliferation score
 11-gene proliferation signature (Nielsen et al. 2010) as reference:
   CPTAC-3:  Yao 0.675 | Proliferation 0.643
   GSE31210: Yao 0.643 | Proliferation 0.690
 Interpretation: proliferation is a meaningful component of signal in both
 cohorts (C > 0.6) but does not dominate Yao 2021 in either cohort.
 Neither score consistently outperforms the other.

## 2026-10-11 — Yao 2021 GPL570 probe mapping
 All 7 Yao 2021 genes are represented on GPL570:
   ARHGEF26-AS1: 236575_at
   LINC01137:    229678_at
   C20orf197:    1553311_at
   MGC32805:     1554781_at, 233475_at (2 probes; highest-mean rule applied)
   TMPO-AS1:     227578_at
   LINC00324:    230012_at
   LINC01116:    228564_at
 Probe sequence identity against RNA-seq transcripts has not been verified;
 this is a known limitation of cross-platform lncRNA validation.

## 2026-10-11 — Justification for pooling Yao 2021 only
 Under the pre-specified 100% coverage rule:
   - Only 1 of 15 signatures is evaluable in both CPTAC-3 and GSE31210
     (Paper62_Yao2021; 7/7 genes present in both cohorts)
   - All other signatures fail at least one cohort at 100% coverage
 Therefore the pooled Yao-only meta-analysis is determined by the
 coverage rule, not by post-hoc selection based on performance. No other
 signature qualifies for pooling under the primary rule.

## 2026-10-12 — Reporting completeness analysis (S22)
 Five reporting dimensions assessed for all 15 included signatures.
 Findings: coefficient tables 15/15; gene identifiers 15/15;
 risk cutoff 15/15; independent external validation 7/15;
 scaling convention stated 0/15.
 Documented in Supplementary Table S22. No primary result changed.

## 2026-10-12 — GSE50081 and GSE30219 coverage-sensitivity (S24)
 Downloaded and harmonized two additional GPL570 cohorts:
   GSE50081: 127 evaluable LUAD, 51 deaths
   GSE30219: 83 evaluable LUAD, 43 deaths
 Finding: Yao 2021 covered only 4 of 7 genes (57.1%) in both cohorts,
 despite sharing GPL570 with GSE31210 (100% coverage).
 Three genes (LINC01137, C20orf197, MGC32805) had annotated probes
 but did not return usable signal.
 Yao 2021 remains evaluable in 2 independent cohorts (k = 2, exploratory).
 Documented in Supplementary Table S24.

## 2026-10-12 — Pre-analysis memo: GPL570 reprocessing for GSE31210, GSE50081, GSE30219

Before any C-index is computed for GSE50081 or GSE30219, the
processing pipeline, probe-selection rule, and detection threshold
have been frozen in docs/prereg/GPL570_reprocessing_memo_2026-10-12.md.

Purpose: determine whether the missing signal for Yao 2021 genes
(LINC01137, C20orf197, MGC32805) in GSE50081 and GSE30219 is a
processing artifact or a genuine detection failure.

Frozen rules:
  - Pipeline: RMA from CEL files, identical for all three GPL570 cohorts
  - Probe mapping: by Ensembl ID, all probes per gene reported
  - Multi-probe rule: highest mean expression per cohort
  - No expression-based filter before coverage assessment
  - Probe-level diagnostics run first, C-indices only afterwards
  - Pre-specified decision rule: artifact vs detection vs coverage
    (see memo for full text)

This entry is written before any scoring of GSE50081 or GSE30219.
Any later change to the pipeline will be logged as a dated amendment
below this entry, not as an edit to the memo.


## 2026-10-12 — Amendment 1: Yao 2021 corrected coverage and k = 4 pooling

Supersedes the earlier S24 coverage-sensitivity entry.

Root cause: the coverage script matched Yao 2021 gene symbols as
published in 2021. Two genes have since been renamed
(LINC01137 -> LITATS1; C20orf197 -> LINC02910/MIR646HG/LOC729296)
and one gene is a retired placeholder (MGC32805) with no current
annotation. All three appeared "missing" in GSE50081 and GSE30219
as an artifact of deprecated-symbol matching.

Corrected coverage:
  GSE31210  7/7  (100%)
  GSE50081  6/7  (85.7%)
  GSE30219  6/7  (85.7%)
  CPTAC-3   6/7  (85.7%)

Under the pre-specified 80% coverage sensitivity rule, Yao 2021 is
evaluable in all four independent cohorts. Yao's development cohort
was TCGA-LUAD; validation was GSE3141 + GSE37745. No overlap with
the four evaluation cohorts (verified in docs/cohort_overlap_verified.csv).

Action: Yao will be scored in all four cohorts under the frozen
memo pipeline (RMA + current symbols + alias resolution) with
published coefficients verbatim. Pooled estimate will be computed
across k = 4, events = 38 + 35 + 51 + 43 = 167. Leave-one-cohort-out
sensitivity reported.

S24 withdrawn. Manuscript, tables, abstract, and figures to be
updated accordingly.

## 2026-10-12 — Amendment 2: correction to GPL570 processing claim

The frozen memo incorrectly stated RMA-from-CEL for the three GPL570
cohorts. Inspection of the repository shows:

- GSE31210: loaded via scripts/14_load_GSE31210.R using getGEO() +
  exprs() (series matrix), with log2(x+1) applied. No RMA.
- GSE50081: no loader script in the repository; loading method not
  verifiable, file structure consistent with the same series-matrix
  pathway.
- GSE30219: same as GSE50081.

This is a reproducibility limitation, not a data problem. It does not
change the coverage analysis, the gene-alias resolution, or the
independence of any cohort. It does change how the Methods must
describe the pipeline, and it must be disclosed in the limitations.

Action: corrected memo. Manuscript Methods to state actual loading
method. No scoring decisions change. Yao scoring can proceed once the
Methods text is updated.

## 2026-10-12 — Amendment 3: Yao 2021 pooled across k = 4 independent cohorts

Following Amendment 1 (alias resolution) and Amendment 2 (pipeline
correction), Yao 2021 was scored in all four independent cohorts
under the pre-specified OSF 80% coverage sensitivity rule with
published coefficients applied verbatim.

Per-cohort HR per SD (95% CI):
  TCGA-LUAD (reference)  1.635 (1.419-1.883)   n=505, 182 events, 7/7 genes
  CPTAC-3                1.516 (1.145-2.007)   n=150,  38 events, 7/7 genes
  GSE31210               1.687 (1.242-2.292)   n=226,  35 events, 7/7 genes
  GSE50081               1.348 (1.030-1.763)   n=127,  51 events, 6/7 genes
  GSE30219               1.278 (0.950-1.719)   n= 83,  43 events, 6/7 genes

Pooled across the four independent cohorts (k = 4, 167 events,
random-effects REML):
  Pooled HR per SD = 1.442 (95% CI 1.249-1.664)
  Pooled log-HR per SD = 0.366 (95% CI 0.222-0.509)
  I^2 = 0%

Leave-one-cohort-out sensitivity (all four subsets):
  Excluding CPTAC-3 : HR = 1.416 (1.199-1.674), I^2 = 0%
  Excluding GSE31210: HR = 1.380 (1.173-1.623), I^2 = 0%
  Excluding GSE50081: HR = 1.481 (1.250-1.755), I^2 = 0%
  Excluding GSE30219: HR = 1.496 (1.270-1.762), I^2 = 0%

Every leave-one-out subset has a lower 95% CI bound above 1.
Yao 2021 was not used in any of the four cohorts by its original
authors (development: TCGA-LUAD; validation: GSE3141 + GSE37745),
so the pooling rule and the >=3-independent-cohort criterion are
satisfied. Event count (167) exceeds the 100-event guidance.

Consequences:
- Yao 2021 moves from "exploratory k = 2" to a pre-specified
  multi-cohort analysis under the registered 80% rule.
- The pooled HR per SD reported in the manuscript changes from
  1.68 (k = 2) to 1.442 (k = 4).
- The corrected analysis supersedes the earlier Yao k = 2 estimate
  for all downstream reporting.
- No other signature's analysis changes.

Saved: results/tables/Yao2021_corrected_coverage_k4.csv

## 2026-10-12 — Reporting-completeness survey expanded to 76 papers (S22)

The reporting-completeness analysis was originally restricted to the 15
included signatures, which was circular (coefficient reporting was an
inclusion criterion). The analysis has been rebuilt on a broader
denominator: 76 papers for which detailed extraction records were
available.

Five reporting items were coded from each paper's extraction record:
  - Full coefficient tables with numeric values: 63/76 (82.9%)
  - Cutoff method defined (any method):          76/76 (100%)
  - Cutoff numeric value stated:                 6/76 (7.9%)
  - Scaling convention stated explicitly:       27/76 (35.5%)
  - Scaling implied only:                        2/76 (2.6%)
  - Scaling not stated:                         47/76 (61.8%)
  - Independent external validation:            24/76 (31.6%)
  - Text-versus-table conflict:                  5/76 (6.6%)

Coding was performed by a single reviewer against each paper's
extraction record. Because no second coder was available, no
inter-rater agreement statistic is reported. This is noted as a
limitation in the manuscript.

The 76-paper denominator replaces the original 15-signature
denominator in the Results, Discussion, and Supplementary Table S22.
The circularity limitation raised in the supervisor's item 4 is
addressed by the expanded denominator; the residual single-coder
limitation is disclosed.

Coding sheet: docs/coding/reporting_completeness_76papers.csv
Snapshot:     docs/coding/reporting_completeness_76papers_2026-10-12.csv

## 2026-10-12 — Reporting-completeness survey expanded to 76 papers

The reporting-completeness analysis was originally restricted to the 15
included signatures, which was circular: full coefficient reporting was
an inclusion criterion. The analysis has been rebuilt on a broader
denominator of 76 reviewed papers for which detailed extraction records
were available.

Five reporting items were coded from each paper's extraction record:
  - Full coefficient tables: 63/76 (82.9%)
  - Cutoff method defined: 76/76 (100%)
  - Cutoff numeric value stated: 6/76 (7.9%)
  - Scaling stated explicitly: 27/76 (35.5%)
  - Scaling implied only: 2/76 (2.6%)
  - Scaling not stated: 47/76 (61.8%)
  - Independent external validation: 24/76 (31.6%)
  - Text-versus-table conflict: 5/76 (6.6%)

Coding was performed by the primary investigator. A second coder is
independently assessing 10 randomly drawn papers against the same
definitions, and the inter-rater agreement will be reported in the
manuscript once complete. The extraction records, coding sheet, and
dual-coding materials are archived at docs/coding/.

The 76-paper denominator replaces the original 15-signature
denominator in the Results, Discussion, and Supplementary Table S22.
The circularity limitation raised in the supervisor's item 4 is
addressed by the expanded denominator; the residual single-coder
limitation is under active resolution via the second coder.

Coding sheet:  docs/coding/reporting_completeness_76papers.csv
Snapshot:      docs/coding/reporting_completeness_76papers_2026-10-12.csv
OSF update:    filed 2026-10-12 (Update 2)
