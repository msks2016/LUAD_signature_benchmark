# Protocol Deviations Log

## Deviation 1: Missing-gene rule relaxation

Original protocol: 20 percent missing-gene rule.
Deviation applied: 25 percent missing-coefficient-mass rule.

Reason: Papers #151 and #305 became not evaluable under 20 percent rule.
Resolution: Applied both rules. Strict (13/15). Relaxed (15/15).

## Deviation 2: Gene name corrections by Ensembl ID

Paper AL034397.3 (ENSG00000272630) = current AL731563.3, NOT AL034397.1
Paper AC010980.2 (ENSG00000272631) = current AC067750.1, NOT AC010980.1
Built v3 matrix with correct mappings. Rescored affected signatures.
Commit: v2.2 (cd8e30a)

## Deviation 3: GEO microarray coverage fallback (2026-10-09)

Trigger: Phase 2 GEO coverage check on GPL570.
Result: Only 2 of 15 signatures (Paper305, Paper62) had 100% coverage.
Average coverage: ~28%. Five signatures had 0%.
Decision: pivot to RNA-seq external cohorts (GSE81089, GSE72094, CPTAC-3).
Table: results/tables/geo_coverage.csv

## Deviation 3: GEO microarray coverage fallback (2026-10-09)

Trigger: Phase 2 GEO coverage check on GPL570.
Result: Only 2 of 15 signatures (Paper305, Paper62) had 100% coverage.
Average coverage: ~28%. Five signatures had 0%.
Decision: pivot to RNA-seq external cohorts (GSE81089, GSE72094, CPTAC-3).
Table: results/tables/geo_coverage.csv

## Deviation 3: GEO microarray coverage fallback (2026-10-09)

Trigger: Phase 2 GEO coverage check on GPL570.
Result: Only 2 of 15 signatures (Paper305, Paper62) had 100% coverage.
Average coverage: ~28%. Five signatures had 0%.
Decision: pivot to RNA-seq external cohorts (GSE81089, GSE72094, CPTAC-3).
Table: results/tables/geo_coverage.csv
