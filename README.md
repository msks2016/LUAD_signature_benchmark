# LUAD Signature Benchmark

Benchmarking published lncRNA/miRNA prognostic signatures in lung adenocarcinoma (LUAD) across independent cohorts.

## Overview

Many lncRNA and miRNA prognostic signatures have been published for LUAD, but few are tested outside their original cohort. This study asks whether they predict survival better than chance, whether they add anything beyond tumor stage, and whether performance drops from the original cohort to independent cohorts.

## Hypotheses

- **H1:** Published signatures perform better than size-matched random gene sets.
- **H2:** Signatures add prognostic value beyond stage (and age/sex where available).
- **H3:** Performance drops from the original cohort to independent cohorts, and the drop differs between signatures.

## Data Sources

| Source | Description |
|--------|-------------|
| TCGA-LUAD | Reference cohort (RNA-seq + clinical data) |
| GSE31210 | Validation cohort |
| GSE50081 | Validation cohort |
| GSE72094 | Validation cohort |
| GSE30219 | Validation cohort |

## Analysis Plan

1. Systematic PubMed search (2018–2025) for LUAD lncRNA/miRNA prognostic signatures
2. Extract signatures with full gene lists and coefficients
3. Download and harmonize cohorts (TCGA + GEO)
4. Apply each signature's published formula unchanged
5. Calculate C-index, time-dependent AUC (1/3/5 years), log-rank test, hazard ratio
6. Compare against 1,000 random signatures of the same size
7. Multivariable Cox + likelihood-ratio test for added value beyond stage
8. Random-effects meta-analysis for cross-cohort heterogeneity
9. Generate figures and tables

## Project Structure
