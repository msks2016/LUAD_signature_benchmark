# Reporting checklist for future LUAD lncRNA prognostic signature studies

Based on the findings of this benchmark and Supplementary Table S22, we
recommend that future studies reporting lncRNA prognostic signatures for
lung adenocarcinoma (LUAD) include the following items in their main text
or supplementary material.

## Signature reporting

- [ ] Full gene list with HGNC symbols (or unambiguous stable identifiers)
- [ ] Complete coefficient table (β value for every gene, with sign)
- [ ] Expression platform used in development (RNA-seq vs microarray,
      with GEO/GDC accession)
- [ ] Explicit statement of the scaling convention applied to expression
      values before coefficient multiplication (raw log2, per-gene z-scoring,
      TPM, or another transformation)
- [ ] Risk-score cutoff or the rule used to define high- vs low-risk groups

## Validation reporting

- [ ] At least one external validation cohort, distinct from the development
      cohort
- [ ] Event count in each validation cohort, alongside sample size
- [ ] lncRNA coverage of the validation platform: fraction of signature
      genes present in the cohort's expression matrix
- [ ] If the coverage rule excludes any signature, state it explicitly
- [ ] C-index and log-HR per SD reported for each validation cohort
- [ ] If pooling across cohorts, state whether fixed- or random-effects
      and report I² with the number of cohorts (k)

## Reporting hygiene

- [ ] Text and table values for coefficients and HRs agree
- [ ] Any post-hoc analysis decisions are disclosed (with date and reason)
- [ ] Pre-registration of the analysis plan, if any, is cited

## Minimal abstract requirements

- [ ] Report the lncRNA coverage of each validation platform
- [ ] Report the event count of each validation cohort
- [ ] Report the scaling convention applied to expression values

---

Source: this checklist is derived directly from the findings of the
reporting-completeness analysis (Supplementary Table S22). Fifteen published
LUAD lncRNA prognostic signatures were assessed. Full coefficients and
gene identifiers were present in 15/15 (100%), a risk cutoff was defined in
15/15 (100%), an independent validation cohort was reported in 7/15 (47%),
and the scaling convention was stated in 0/15 (0%).
