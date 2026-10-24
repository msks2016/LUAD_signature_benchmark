# Deviations log, v2 pipeline

Registration: osf.io/qckvr (filed 2026-10-06). Updates: 1 (2026-10-10), 2 (2026-10-12), 3 (2026-10-20).
The v0.3 deviations log (docs/deviations_log.md) records everything before the v2 rebuild; the results seen with the v0.3 code are listed in Update 3 ("Explanation of foreknowledge"). This log records what happened after Update 3 was filed. Each entry gives the date, what changed, why, and whether any v2 outcome had been seen.

Freeze point: git tag `v2-preanalysis` (commit 871d199, 2026-10-24), set after the gene mapping and coverage were final and before any v2 outcome analysis (scripts 29-37) was run.

## A. Decisions after Update 3 and before any v2 outcome analysis (2026-10-20 to 2026-10-24)

| # | Date | Change | Reason | Outcome seen? |
|---|------|--------|--------|---------------|
| A1 | 2026-10-20 | Independent double extraction of all 102 coefficients completed (two extractors, blind). 86 identical, 16 differed only by rounding (one extractor rounded the printed value), 0 real disagreements. The full-precision printed value is used. | Planned in Update 3. | No |
| A2 | 2026-10-20 | Gene-name release per signature fixed by a check of which GENCODE table knows the most of each paper's gene names: GENCODE v22 for Yu 2024 and Yao 2021, v36 for the other 13. | Update 3 change 1 (map each gene with the annotation its paper used). | No |
| A3 | 2026-10-20 | Six gene-signature rows mapped by documented manual override (config/gene_overrides.csv): AL034397.3 (Jin 2020, Song 2021) -> ENSG00000274536; HSPC324 (Jin 2020) -> ENSG00000228401; AC026355.1 (Zheng 2021, Zhang 2025) -> ENSG00000228561; AC010980.2 (Song 2021) -> ENSG00000267034. Evidence per row in the file (GENCODE v29-v36 name history, NCBI Gene aliases). The last two corrected matches to a snoRNA and a pseudogene that reused the clone-based names in other GENCODE releases. | A name audit across GENCODE v22, v29, v31, v32 and v36 (scripts 22b, 22c) found names that point to different genes in different releases. | No |
| A4 | 2026-10-20 | Successor Ensembl IDs (config/gene_id_alternatives.csv, 7 genes) from the official Ensembl ID history: used in a cohort only where the original ID is not measurable there and the successor is. | Retired or merged IDs made present genes look absent (e.g. LINC00941, C20orf197/MIR646HG). Not in the registration. | No |
| A5 | 2026-10-20 / 24 | Predecessor check for GSE81089 (Ensembl 73): older IDs whose official history leads to a missing signature gene. Decision rule fixed on 2026-10-20 before any outcome: accept only if measurable, the history passes through the signature gene's own ID, ends in one current gene, and is the only such ID. Result (2026-10-24): 63,054 IDs queried, 10,024 retired, none accepted. | Same rule as A4 in the other direction. | No |
| A6 | 2026-10-19 / 20 | Probe and Entrez maps restricted to primary-assembly genes (chromosomes 1-22, X, Y, MT). | Alternative-haplotype copies made probes look multi-mapping and dropped them. | No |
| A7 | 2026-10-20 | GSE72094 built from the GEO series matrix; Entrez IDs mapped to Ensembl (release 114). 398 patients, 113 deaths. No signature reached 80% coverage. | Update 3: analysed for any signature reaching 80% coverage. | No |
| A8 | 2026-10-20 | Independence table completed by two extractors from the full texts: Geng 2021 used GSE31210; Yu 2024 used GSE31210, GSE30219 and GSE50081 (with GSE37745); Yang 2022 used GSE50081 (single-gene validation; classed as used). No paper used Kaplan-Meier Plotter or a similar portal for survival validation (Yao 2021 used GEPIA for expression correlation only). | Update 3 clarification 12. | No |

## B. After the freeze (2026-10-24)

| # | Date | Change | Effect on numbers |
|---|------|--------|-------------------|
| B1 | 2026-10-24 | Code fix in the power simulation: exp(beta x) overflowed inside integrate() at a true C-index of 0.70; rewritten as plogis(beta x - log mu), the same quantity. | None for C = 0.60 and 0.65 (identical); C = 0.70 could not run before. |
| B2 | 2026-10-24 | Speed changes in the transportability bootstrap: signature genes are subset before resampled patients; bootstrap replicates use survival::coxph.fit with the inputs coxph() passes (aeqSurv time fix, zero offset, Efron ties, nocenter = c(-1, 0, 1)). Progress messages added. | None: point estimates identical to scripts 30/32 (difference 0); equality with coxph() tested. |
| B3 | 2026-10-24 | Pooling-interval sensitivity (script 37): common-effect, REML Wald and modified Hartung-Knapp (q* = max(1, q)) intervals beside the registered REML + Hartung-Knapp interval. Added after the registered interval for Yao 2021 was seen to be narrower than the common-effect interval (q = Q/(k-1) = 0.19). | Reported as a post hoc sensitivity analysis; the registered interval stays primary. Yao 2021: HR 1.41 per SD; registered 1.30-1.52; modified 1.17-1.68; common effect 1.24-1.60. |
| B4 | 2026-10-24 | Manuscript outputs (script 40): median follow-up by the reverse Kaplan-Meier method; Kaplan-Meier curves at the median score. | Descriptive only. |
| B5 | 2026-10-24 | v1 files merged in from the old repository moved to archive_v1/; none was used by any v2 script. | None. |

## C. To be filed

OSF Update 4 (post-analysis disclosure): sections A and B of this log, stating that A1-A8 were decided before any v2 outcome analysis (tag v2-preanalysis) and that B1-B5 came after it, with B3 reported as post hoc.
