## Acute check s5 (RUN_RECORD choice 6): TE share of TEcount-assigned counts per library, lm(share ~ cohort + treatment) per tissue;
## treatment coefficient in percentage points with a 95 % t interval. Descriptive. Usage: R --file=acute_share.R --args ANALYSIS_DIR
a <- commandArgs(trailingOnly = TRUE); O <- a[1]; out <- NULL
for (t in c("brain", "gill")) {
  d <- read.delim(file.path(O, sprintf("%s_te_share_input.tsv", t))); d$share <- 100 * d$te_set / d$assigned
  d$cohort <- factor(d$cohort); d$treatment <- factor(d$treatment, levels = c("Control", "CO2"))
  f <- lm(share ~ cohort + treatment, d); ci <- confint(f)["treatmentCO2", ]
  out <- rbind(out, data.frame(tissue = t, n = nrow(d), mean_control = mean(d$share[d$treatment == "Control"]), mean_co2 = mean(d$share[d$treatment == "CO2"]),
                               effect_pp = coef(f)["treatmentCO2"], ci_lo = ci[1], ci_hi = ci[2], p = summary(f)$coefficients["treatmentCO2", 4]))
  write.table(d, file.path(O, sprintf("%s_te_share.tsv", t)), sep = "\t", quote = FALSE, row.names = FALSE)
}
write.table(out, file.path(O, "te_share_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(out, digits = 3)
