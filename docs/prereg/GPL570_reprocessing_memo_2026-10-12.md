# Pre-analysis memo: GPL570 reprocessing for GSE31210, GSE50081, GSE30219

Date written: 2026-10-12
Date any result was seen: (none — this memo precedes all scoring)
Author: [your name]

## Purpose

Resolve whether the missing signal for Yao 2021 genes (LINC01137,
C20orf197, MGC32805) in GSE50081 and GSE30219 is a processing
artifact or a genuine detection failure.

## Frozen processing pipeline

All three GPL570 cohorts will be processed identically:

1. Source: raw CEL files from GEO, if available; otherwise the
   submitter-supplied series matrix.
2. Normalization: RMA via affy::rma() using the GPL570 CDF
   (hgu133plus2cdf). If CEL files are unavailable, log that
   fact and use the series matrix unchanged.
3. Log2 transformation: RMA output is already log2; series matrix
   values will be used as provided (log2 scale).
4. Probe-to-gene mapping: match by Ensembl ID via the GPL570
   annotation (hgu133plus2.db), not by old HGNC aliases.
5. Multi-probe rule: for each gene, retain the probe with the
   highest mean expression across samples in that cohort.
6. Detection / filtering: apply NO expression-based filter
   before coverage assessment. Detection status will be assessed
   separately by intensity distribution.

## Frozen probe-level diagnostics (to run before any C-index)

For each of the three genes (LINC01137, C20orf197, MGC32805) in
each of the three cohorts:
- Confirm the probe ID(s) exist as rows in the processed matrix.
- Count missing values per probe per cohort.
- Report mean, median, SD, and min-max of each probe.
- Compare distributions across cohorts (boxplot saved as figure).
- For genes with multiple probes, report ALL probes, not only the
  highest-expressing one.

## Variance guard audit

Count how many genes were removed by the variance guard
(SD < 1e-6) in each cohort. Confirm none of the three target
genes were among them.

## Decision rule (pre-specified)

- If the three probes are present, non-missing, and their
  intensities in GSE50081/GSE30219 overlap the GSE31210
  distribution: the gap is a PROCESSING artifact; reprocess
  identically and re-score Yao.
- If the three probes are present but their intensities in
  GSE50081/GSE30219 are at or near background while GSE31210
  is well above: the gap is a DETECTION finding; report with
  intensity plots; do not call it a coverage ceiling.
- If the probes are absent as rows: the gap is a COVERAGE
  finding; report as such.

## Deviation log

This memo and the decision rule above will be entered into
docs/deviations_log.md before any C-index is computed for
GSE50081 or GSE30219.

---

## Amendment 1 — 2026-10-12

Written after the original frozen memo, before any Yao scoring
of GSE50081 or GSE30219.

New evidence:

1. Signature verification (data/processed/signature_coefficients.csv):
   Paper62 (Yao 2021) stored coefficients match the published
   Table 2 exactly for all 7 genes (ARHGEF26-AS1 -0.310,
   LINC01137 +0.364, C20orf197 -0.270, MGC32805 -0.263,
   TMPO-AS1 +0.331, LINC00324 -0.404, LINC01116 +0.222).

2. Gene-symbol renaming since Yao 2021 publication:
   - LINC01137 = LITATS1 (Ensembl ENSG00000233621)
   - C20orf197 = LINC02910 / MIR646HG / LOC729296
                 (Ensembl ENSG00000228340)
   Both genes are present in GSE50081 and GSE30219 under their
   current symbols.

3. MGC32805 is a retired placeholder symbol with no current
   HGNC/Ensembl entry and no probe in current GPL570 annotation.
   Present as a row name only in GSE31210 (older annotation).
   Absent from GSE50081, GSE30219, and CPTAC-3.

4. Cohort independence (docs/cohort_overlap_verified.csv):
   Paper62 development cohort = TCGA-LUAD.
   Author validation cohorts  = GSE3141 + GSE37745.
   No overlap with GSE31210, GSE50081, GSE30219, or CPTAC-3.
   All four cohorts are independent for Yao 2021.

Decision (pre-specified under the OSF 80% coverage sensitivity rule,
registered before any cohort was scored):

Yao 2021 will be evaluated in all four independent cohorts under
the 80% coverage rule with published coefficients applied verbatim.
The MGC32805 term is omitted in cohorts where the gene is absent
(GSE50081, GSE30219, CPTAC-3); the remaining 6 genes use published
coefficients unchanged. This is a coverage-rule application, not
a signature modification.

Yao's pooled log-HR per SD and HR per SD will be computed across
all four independent cohorts, with a leave-one-cohort-out
sensitivity. The pooled estimate supersedes the exploratory k = 2
estimate in the current manuscript.

The S24 "within-platform coverage ceiling" claim is withdrawn.
The apparent 57.1% coverage in GSE50081 and GSE30219 was an
artifact of matching by deprecated HGNC symbols, not a real
platform or detection limitation.

---

## Amendment 2 — 2026-10-12

Amendment 1 corrected gene-symbol handling. This amendment corrects the
processing-pipeline claim in the original memo.

The original memo stated: "Pipeline: RMA from CEL files, identical for
all three GPL570 cohorts." Inspection of the repository shows this is
inaccurate.

Actual state, per script and file inspection:

- GSE31210: loaded via scripts/14_load_GSE31210.R using
  getGEO("GSE31210", GSEMatrix = TRUE) and exprs(). The submitter-
  normalized series matrix was used. log2(x + 1) was applied after
  exprs(). No RMA and no CEL files were used.

- GSE50081: no loader script exists in the repository. The file
  data/processed/GSE50081_expr_mapped.rds exists and was used in
  analysis. Its row names are gene symbols with an Affymetrix probe-ID
  attribute, matching the mapping pattern of 14_load_GSE31210.R, but
  the exact loading and mapping steps cannot be verified from the
  repository.

- GSE30219: same as GSE50081. No loader script; unverifiable from the
  repository.

Consequences:

1. The three GPL570 cohorts were not reprocessed with RMA from CEL
   files. They were taken from the GEO series matrices (confirmed for
   GSE31210; highly likely for the other two, per file structure).
2. The three matrices may therefore differ in preprocessing, batch
   effects, and probe-level summarization in ways that are not
   controlled. This is a limitation to state explicitly.
3. No cross-cohort pooled analysis should be described as "RMA-
   processed uniformly" because that is not what was done.

Recommendation for the manuscript: state the actual loading method
(GEO series matrices via getGEO) in the Methods, disclose that
GSE50081 and GSE30219 loaders are not in the repository, and note
RMA-from-CEL as a stronger pipeline for future work.
