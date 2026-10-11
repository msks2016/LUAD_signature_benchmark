
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
