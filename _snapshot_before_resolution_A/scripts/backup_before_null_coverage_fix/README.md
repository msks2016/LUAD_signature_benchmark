# Backup scripts — historical snapshots

These folders contain the state of pipeline scripts immediately before specific
corrections were applied. They are kept for audit transparency and to document
the deviations described in the deviations log (see docs/deviations_log.md).

- backup_before_reverse_fix/ — before adding reverse = TRUE to concordance() calls
- backup_before_raw_expression_fix/ — before switching from z-scored to raw log2 expression
- backup_before_null_coverage_fix/ — before enforcing the 100% coverage rule in the null test
- backup_before_score_refit/ — before the (subsequently reverted) score-level refit attempt
- backup_before_transport_raw_fix/ — before switching the transportability script to raw expression

All primary scripts in scripts/ reflect the final, corrected pipeline.
