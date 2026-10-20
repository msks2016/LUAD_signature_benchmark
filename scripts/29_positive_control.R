# 29_positive_control.R - stage positive control in every cohort (registered; Update 3, D1).
# A cohort is valid for signature evaluation if the Cox model on stage (I = reference) shows
# HR(III vs I) > 1 with p < 0.05 and a monotonic HR across stages. In a cohort with fewer than
# 10 stage III-IV patients that criterion cannot be evaluated, and HR(II vs I) > 1 with
# p < 0.05 is used instead. A cohort that fails is still reported; its signature metrics are
# flagged exploratory in scripts 30-32.
# Output: results/tables/positive_control.csv

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
cohort_names <- c(cfg$cohorts$reference, cfg$cohorts$independent)
cohorts <- load_cohorts(cohort_names)

res <- do.call(rbind, lapply(names(cohorts), function(n) {
  co <- cohorts[[n]]
  r <- positive_control(co$time, co$event, co$clin$stage,
                        min_iii_iv = cfg$positive_control$min_stage_iii_iv, alpha = cfg$alpha)
  data.frame(cohort = n, role = if (n %in% cfg$cohorts$reference) "reference" else "independent",
             r, stringsAsFactors = FALSE)
}))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, "results/tables/positive_control.csv", row.names = FALSE)
for (i in seq_len(nrow(res))) {
  record_number(paste0("positive_control_", res$cohort[i]), if (res$pass[i]) "pass" else "fail",
                note = res$rule[i])
}
print(res[, c("cohort", "n_by_stage", "rule", "hr_ii", "p_ii", "hr_iii", "p_iii", "monotonic", "pass", "c_stage")])
if (any(!res$pass)) message("Cohorts failing the stage control are kept and flagged exploratory: ",
                            paste(res$cohort[!res$pass], collapse = ", "))
