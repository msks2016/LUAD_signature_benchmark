# Reporting checklist for lncRNA prognostic signature studies in lung adenocarcinoma

Supplementary File 2. Items that a reader needs to apply a published signature to new data, and that a validation study needs to report so that its result can be interpreted. Derived from the reporting-completeness analysis of 76 published papers (Supplementary Table S19) and from the problems met when applying 15 signatures to seven cohorts in this benchmark.

## Signature reporting

- [ ] Every signature gene with a stable identifier: Ensembl gene ID with version, plus the gene symbol
- [ ] The gene annotation release used (e.g. GENCODE v36), because clone-based lncRNA names are reassigned between releases
- [ ] The coefficient of every gene, with sign, to the precision used in the risk score; the same values in text, formula and tables
- [ ] The expression scale to which the coefficients apply (e.g. log2(TPM + 1), log2 FPKM, z-scored within cohort) and the quantification pipeline
- [ ] The numeric risk-score cutoff, if risk groups are reported, and the rule used to choose it
- [ ] The development cohort, with its accession, platform and library type (poly(A) or total RNA)

## Validation reporting

- [ ] Every validation cohort, with its accession, and whether it was used for gene selection or tuning
- [ ] The number of patients and of events in each validation cohort
- [ ] Gene coverage in each cohort: how many signature genes were measurable, and how missing genes were handled
- [ ] Discrimination (C-index with 95% CI) and the hazard ratio per SD of the score in each cohort
- [ ] A comparison with random gene sets of the same size, scored with the same coefficients, in addition to a test against no association
- [ ] Added value beyond stage and other clinical predictors
- [ ] If cohorts are pooled: the model, the number of cohorts (k), heterogeneity, and the interval method

## Reporting hygiene

- [ ] Pre-registration of the analysis plan, if any, with its updates
- [ ] Deviations from the plan, with dates and reasons
- [ ] Code and gene mappings available, so that the score can be recomputed

## Abstract

- [ ] Number of events and gene coverage in each validation cohort
- [ ] The expression scale of the coefficients
