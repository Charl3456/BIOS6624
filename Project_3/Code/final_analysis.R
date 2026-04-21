library(survival)
library(survminer)
library(gtsummary)   
library(gt)          
library(flextable)   
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
library(xtable)

dir.create("../Results/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("../Results/Tables",  recursive = TRUE, showWarnings = FALSE)

fig_dir <- "../Results/Figures/"
tab_dir <- "../Results/Tables/"

frmgham <- read.csv("../RawData/frmgham2.csv")
cat("Full dataset:", nrow(frmgham), "observations,",
    length(unique(frmgham$RANDID)), "unique participants\n")
ALPHA_SELECT <- 0.10
# Helper: save a standard data frame table
save_table <- function(df, filename, caption = "", digits = 3) {
  
  # .tex
  tex_path <- paste0(tab_dir, filename, ".tex")
  xt <- xtable(df, caption = caption, digits = digits)
  print(xt, file = tex_path,
        include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = "top",
        sanitize.text.function = identity)
  cat("  Saved:", tex_path, "\n")
  
  png_path <- paste0(tab_dir, filename, ".png")
  tt <- ttheme_minimal(
    core    = list(fg_params = list(fontsize = 9),
                   padding   = unit(c(4, 4), "mm")),
    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"),
                   padding   = unit(c(4, 4), "mm"))
  )
  tbl_grob <- tableGrob(df, rows = NULL, theme = tt)
  
  w     <- convertWidth( sum(tbl_grob$widths),  "in", valueOnly = TRUE)
  h_tbl <- convertHeight(sum(tbl_grob$heights), "in", valueOnly = TRUE)
  
  h_cap   <- if (nchar(caption) > 0) 0.5 else 0
  h_total <- max(h_tbl + h_cap + 0.3, 2)
  w_total <- max(w + 0.8, 8)
  
  png(png_path, width = w_total, height = h_total,
      units = "in", res = 200, bg = "white")
  grid.newpage()
  if (nchar(caption) > 0) {
    grid.text(caption,
              x    = unit(0.4, "in"),
              y    = unit(h_total - 0.28, "in"),
              gp   = gpar(fontsize = 11, fontface = "bold"),
              just = "left",
              default.units = "in")
  }
  
  # Draw table below caption, horizontally centered
  x_off <- (w_total - w) / 2
  pushViewport(viewport(
    x      = unit(x_off, "in"),
    y      = unit(0.15, "in"),
    width  = unit(w,     "in"),
    height = unit(h_tbl, "in"),
    just   = c("left", "bottom")
  ))
  grid.draw(tbl_grob)
  popViewport()
  
  dev.off()
  cat("  Saved:", png_path, "\n")
}

# Helper: save a ggplot figure
save_fig <- function(plot_obj, filename, w = 8, h = 5) {
  path <- paste0(fig_dir, filename, ".png")
  ggsave(path, plot_obj, width = w, height = h, dpi = 200, bg = "white")
  cat("  Saved:", path, "\n")
}

# Define Analytic Cohort
baseline <- frmgham %>%
  filter(PERIOD == 1, PREVSTRK == 0)

TEN_YEARS <- 3652

baseline <- baseline %>%
  mutate(
    time_stroke_10yr = pmin(TIMESTRK, TEN_YEARS),
    stroke_10yr      = ifelse(STROKE == 1 & TIMESTRK <= TEN_YEARS, 1, 0),
    time_death_10yr  = pmin(TIMEDTH, TEN_YEARS),
    death_10yr       = ifelse(DEATH == 1 & TIMEDTH <= TEN_YEARS, 1, 0),
    # Competing event: death without prior stroke
    death_no_stroke  = ifelse(death_10yr == 1 & stroke_10yr == 0, 1, 0),
    time_years       = time_stroke_10yr / 365.25,
    sex_label        = factor(SEX, levels = c(1, 2), labels = c("Men", "Women")),
    diabetes_label   = factor(DIABETES, levels = c(0, 1),
                              labels = c("No Diabetes", "Diabetes")),
    # Labeled factor versions for clean gtsummary display
    diabetes_f = factor(DIABETES,    levels = c(0, 1), labels = c("No", "Yes")),
    bpmeds_f   = factor(BPMEDS,      levels = c(0, 1), labels = c("No", "Yes")),
    cursmoke_f = factor(CURSMOKE,    levels = c(0, 1), labels = c("No", "Yes")),
    prevchd_f  = factor(PREVCHD,     levels = c(0, 1), labels = c("No", "Yes")),
    stroke_f   = factor(stroke_10yr, levels = c(0, 1), labels = c("No", "Yes")),
    death_f    = factor(death_10yr,  levels = c(0, 1), labels = c("No", "Yes"))
  )

tab1 <- baseline %>%
  select(
    sex_label,
    AGE, SYSBP, TOTCHOL, BMI,
    diabetes_f, bpmeds_f, cursmoke_f, prevchd_f,
    stroke_f, death_f
  ) %>%
  tbl_summary(
    by = sex_label,
    label = list(
      AGE        ~ "Age (years)",
      SYSBP      ~ "Systolic BP (mmHg)",
      TOTCHOL    ~ "Total Cholesterol (mg/dL)",
      BMI        ~ "BMI (kg/m\u00b2)",
      diabetes_f ~ "Diabetes",
      bpmeds_f   ~ "Antihypertensive Medication",
      cursmoke_f ~ "Current Smoker",
      prevchd_f  ~ "Prevalent CHD",
      stroke_f   ~ "Stroke within 10 years",
      death_f    ~ "Death within 10 years"
    ),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    # show missing count/% for any variable with missing values (e.g. BPMEDS)
    missing          = "ifany",
    missing_text     = "Missing, n (%)",
    missing_stat     = "{N_miss} ({p_miss}%)"
  ) %>%
  add_overall(last = FALSE) %>%
  modify_caption("**Table 1. Baseline Characteristics by Sex**") %>%
  modify_spanning_header(c(stat_1, stat_2) ~ "**Sex**") %>%
  bold_labels()

# Print to console
tab1

# Save Table 1 as .tex
tab1_tex <- tab1 %>%
  as_kable_extra(format = "latex", booktabs = TRUE,
                 caption = "Baseline Characteristics by Sex")
writeLines(tab1_tex, paste0(tab_dir, "Table1_Baseline.tex"))
cat("  Saved:", paste0(tab_dir, "Table1_Baseline.tex"), "\n")

# Save Table 1 as .png via flextable
tab1 %>%
  as_flex_table() %>%
  flextable::bg(bg = "white", part = "all") %>%
  flextable::save_as_image(
    path = paste0(tab_dir, "Table1_Baseline.png"),
    res  = 200
  )

# Kaplan-Meier Curves
# KM by Sex
km_sex <- survfit(Surv(time_years, stroke_10yr) ~ sex_label, data = baseline)
p_km_sex <- ggsurvplot(
  km_sex, data = baseline,
  risk.table = TRUE, pval = TRUE, conf.int = TRUE,
  xlim = c(0, 10), ylim = c(0.90, 1.0), break.time.by = 2,
  palette = c("#2171b5", "#cb181d"),
  xlab = "Time (years)", ylab = "Stroke-Free Survival Probability",
  title = "Kaplan-Meier: Stroke-Free Survival by Sex",
  legend.labs = c("Men", "Women"),
  risk.table.height = 0.25, ggtheme = theme_minimal()
)
png(paste0(fig_dir, "Figure1_KM_by_Sex.png"),
    width = 8, height = 6, units = "in", res = 200, bg = "white")
print(p_km_sex)
dev.off()
cat("  Saved:", paste0(fig_dir, "Figure1_KM_by_Sex.png"), "\n")

# KM by Diabetes, STRATIFIED BY SEX 
km_dm_sex <- survfit(Surv(time_years, stroke_10yr) ~ diabetes_label + sex_label,
                     data = baseline)
sdf_dm <- surv_summary(km_dm_sex, data = baseline)
sdf_dm$group <- sdf_dm$diabetes_label
p_km_dm_sex <- ggplot(sdf_dm, aes(x = time, y = surv, color = group, fill = group)) +
  geom_step(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  facet_wrap(~ sex_label) +
  coord_cartesian(xlim = c(0,10), ylim = c(0.80, 1.0)) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_color_manual(values = c("#2ca02c", "#d62728")) +
  scale_fill_manual (values = c("#2ca02c", "#d62728")) +
  labs(x = "Time (years)", y = "Stroke-Free Survival Probability",
       color = "Diabetes Status", fill = "Diabetes Status",
       title = "Kaplan-Meier: Stroke-Free Survival by Diabetes, Stratified by Sex") +
  theme_minimal() + theme(legend.position = "bottom")
save_fig(p_km_dm_sex, "Figure2_KM_by_Diabetes_bySex", w = 10, h = 5)

# Additional KM curves for other candidate risk factors
# All stratified by sex to match the primary analysis structure.
baseline <- baseline %>%
  mutate(
    smoke_label   = factor(CURSMOKE, levels = c(0,1),
                           labels = c("Non-smoker","Current smoker")),
    prevchd_label = factor(PREVCHD,  levels = c(0,1),
                           labels = c("No prior CHD","Prevalent CHD")),
    bpmeds_label  = factor(BPMEDS,   levels = c(0,1),
                           labels = c("No BP meds","On BP meds")),
    # dichotomize age and SBP at clinically meaningful cut-points
    age_cat       = factor(ifelse(AGE   >= 55,  "Age 55+",  "Age under 55"),
                           levels = c("Age under 55", "Age 55+")),
    sbp_cat       = factor(ifelse(SYSBP >= 140, "SBP 140+", "SBP under 140"),
                           levels = c("SBP under 140", "SBP 140+"))
  )

km_extra <- function(var, ylim_low, fname, title_lab) {
  # Build survfit with the variable passed by name
  f   <- reformulate(c(var, "sex_label"), response = quote(Surv(time_years, stroke_10yr)))
  fit <- survfit(f, data = baseline)
  
  # Extract tidy survival curve data and plot manually with facet_wrap by sex.
  sdf <- surv_summary(fit, data = baseline)
  # surv_summary returns columns for each stratum variable; rename for plotting
  sdf$group <- sdf[[var]]
  p <- ggplot(sdf, aes(x = time, y = surv, color = group, fill = group)) +
    geom_step(linewidth = 0.7) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    facet_wrap(~ sex_label) +
    coord_cartesian(xlim = c(0, 10), ylim = c(ylim_low, 1.0)) +
    scale_x_continuous(breaks = seq(0, 10, 2)) +
    labs(x = "Time (years)", y = "Stroke-Free Survival Probability",
         color = NULL, fill = NULL,
         title = paste0("Kaplan-Meier: ", title_lab, ", Stratified by Sex")) +
    theme_minimal() + theme(legend.position = "bottom")
  save_fig(p, fname, w = 10, h = 5)
}
km_extra("smoke_label",   0.85, "Figure2b_KM_by_Smoking_bySex",   "Stroke-Free Survival by Smoking Status")
km_extra("prevchd_label", 0.70, "Figure2c_KM_by_PrevCHD_bySex",   "Stroke-Free Survival by Prevalent CHD")
km_extra("bpmeds_label",  0.80, "Figure2d_KM_by_BPMeds_bySex",    "Stroke-Free Survival by BP Medication Use")
km_extra("age_cat",       0.80, "Figure2e_KM_by_Age_bySex",       "Stroke-Free Survival by Age (<55 vs >=55)")
km_extra("sbp_cat",       0.85, "Figure2f_KM_by_SBP_bySex",       "Stroke-Free Survival by SBP (<140 vs >=140)")

# Cox PH — Base Models
sure_vars  <- c("AGE", "DIABETES", "SYSBP")
quest_vars <- c("PREVCHD", "BPMEDS", "CURSMOKE", "TOTCHOL", "BMI")

f_base <- as.formula(
  paste("Surv(time_years, stroke_10yr) ~", paste(sure_vars, collapse = " + "))
)

men   <- baseline %>% filter(SEX == 1)
women <- baseline %>% filter(SEX == 2)

cox_men_base   <- coxph(f_base, data = men)
cox_women_base <- coxph(f_base, data = women)

# Helper: extract Cox model results
var_labels <- c(
  AGE      = "Age (years)",
  DIABETES = "Diabetes",
  SYSBP    = "Systolic BP (mmHg)",
  CURSMOKE = "Current Smoker",
  PREVCHD  = "Prevalent CHD",
  BPMEDS   = "BP Medications",
  TOTCHOL  = "Total Cholesterol (mg/dL)",
  BMI      = "BMI (kg/m\u00b2)"
)

tidy_cox <- function(model) {
  s    <- summary(model)
  vars <- rownames(s$conf.int)
  data.frame(
    Variable = ifelse(vars %in% names(var_labels), var_labels[vars], vars),
    HR       = sprintf("%.3f", s$conf.int[, "exp(coef)"]),
    `95% CI` = sprintf("(%.3f\u2013%.3f)",
                       s$conf.int[, "lower .95"],
                       s$conf.int[, "upper .95"]),
    `Wald p` = ifelse(s$coefficients[, "Pr(>|z|)"] < 0.001, "<0.001",
                      sprintf("%.4f", s$coefficients[, "Pr(>|z|)"])),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

save_table(tidy_cox(cox_men_base), "Table2_BaseModel_Men",
           caption = sprintf("Table 2. Base Cox Model - Men (C = %.3f)",
                             summary(cox_men_base)$concordance["C"]))
save_table(tidy_cox(cox_women_base), "Table3_BaseModel_Women",
           caption = sprintf("Table 3. Base Cox Model - Women (C = %.3f)",
                             summary(cox_women_base)$concordance["C"]))

# Individual Testing of Candidate Risk Factors
individual_test <- function(data, sure_vars, quest_vars, sex_label) {
  
  results <- data.frame(
    Variable = character(), HR = numeric(),
    CI_lower = numeric(), CI_upper = numeric(),
    Wald_p = numeric(), LRT_p = numeric(),
    Concordance = numeric(), stringsAsFactors = FALSE
  )
  
  for (v in quest_vars) {
    all_needed  <- c(sure_vars, v, "time_years", "stroke_10yr")
    cc_data     <- data[complete.cases(data[, all_needed]), ]
    
    f_base_cc <- as.formula(paste("Surv(time_years, stroke_10yr) ~",
                                  paste(sure_vars, collapse = " + ")))
    cox_base_cc <- coxph(f_base_cc, data = cc_data)
    
    f_test    <- as.formula(paste("Surv(time_years, stroke_10yr) ~",
                                  paste(c(sure_vars, v), collapse = " + ")))
    cox_test  <- coxph(f_test, data = cc_data)
    
    lr      <- anova(cox_base_cc, cox_test, test = "LRT")
    lr_pval <- lr[["Pr(>|Chi|)"]][2]
    
    s   <- summary(cox_test)
    idx <- which(rownames(s$conf.int) == v)
    
    results <- rbind(results, data.frame(
      Variable = v,
      HR       = s$conf.int[idx, "exp(coef)"],
      CI_lower = s$conf.int[idx, "lower .95"],
      CI_upper = s$conf.int[idx, "upper .95"],
      Wald_p   = s$coefficients[idx, "Pr(>|z|)"],
      LRT_p    = lr_pval,
      Concordance = s$concordance["C"],
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("  Base + %-10s  n=%d  HR=%.3f (%.3f-%.3f)  Wald p=%.4f  LRT p=%.4f  C=%.4f %s\n",
                v, nrow(cc_data),
                s$conf.int[idx, "exp(coef)"],
                s$conf.int[idx, "lower .95"],
                s$conf.int[idx, "upper .95"],
                s$coefficients[idx, "Pr(>|z|)"],
                lr_pval, s$concordance["C"],
                ifelse(lr_pval < ALPHA_SELECT, " *", "")))
  }
  
  base_cc_all <- data[complete.cases(data[, c(sure_vars, "time_years", "stroke_10yr")]), ]
  base_c <- summary(coxph(f_base_cc, data = base_cc_all))$concordance["C"]
  cat(sprintf("  Base model C-index: %.4f\n", base_c))
  
  return(results)
}

results_men   <- individual_test(men,   sure_vars, quest_vars, "MEN")
results_women <- individual_test(women, sure_vars, quest_vars, "WOMEN")

format_indiv <- function(results) {
  results %>%
    mutate(
      # Replace ALL CAPS variable names with readable labels
      Variable  = ifelse(Variable %in% names(var_labels),
                         var_labels[Variable], Variable),
      HR_CI     = sprintf("%.3f (%.3f\u2013%.3f)", HR, CI_lower, CI_upper),
      `Wald p`  = ifelse(Wald_p < 0.001, "<0.001", sprintf("%.4f", Wald_p)),
      `LRT p`   = ifelse(LRT_p  < 0.001, "<0.001", sprintf("%.4f", LRT_p)),
      `C-index` = sprintf("%.4f", Concordance),
      Sig       = ifelse(LRT_p < ALPHA_SELECT, "*", "")
    ) %>%
    select(Variable, `HR (95% CI)` = HR_CI, `Wald p`, `LRT p`, `C-index`, Sig)
}

save_table(format_indiv(results_men), "Table4_IndivTest_Men",
           caption = sprintf("Table 4. Individual Variable Testing - Men (* LRT p < %.2f)",
                             ALPHA_SELECT))
save_table(format_indiv(results_women), "Table5_IndivTest_Women",
           caption = sprintf("Table 5. Individual Variable Testing - Women (* LRT p < %.2f)",
                             ALPHA_SELECT))

# SENSITIVITY: joint full model with all candidate variables 
full_formula <- as.formula(
  paste("Surv(time_years, stroke_10yr) ~",
        paste(c(sure_vars, quest_vars), collapse = " + ")))
cat("\n--- SENSITIVITY: Joint full model, MEN ---\n")
print(summary(coxph(full_formula, data = men)))
cat("\n--- SENSITIVITY: Joint full model, WOMEN ---\n")
print(summary(coxph(full_formula, data = women)))

# Build shared complete-case dataset so AIC is comparable across steps
all_vars <- c(sure_vars, quest_vars)
men_cc   <- men[complete.cases(men[, c(all_vars, "time_years", "stroke_10yr")]), ]
women_cc <- women[complete.cases(women[, c(all_vars, "time_years", "stroke_10yr")]), ]
cat(sprintf("\nComplete-case N: Men = %d, Women = %d\n", nrow(men_cc), nrow(women_cc)))

f_lower <- as.formula(paste("Surv(time_years, stroke_10yr) ~",
                            paste(sure_vars, collapse = " + ")))
f_upper <- as.formula(paste("Surv(time_years, stroke_10yr) ~",
                            paste(all_vars,  collapse = " + ")))

run_stepwise <- function(cc_data, sex_tag) {
  cat("\n==================== ", sex_tag, " ====================\n", sep = "")
  
  # Backward: start at full model, remove toward lower (sure_vars)
  full_fit <- coxph(f_upper, data = cc_data)
  cat("\n--- Backward elimination (", sex_tag, ") ---\n", sep = "")
  back_fit <- step(full_fit, scope = list(lower = f_lower, upper = f_upper),
                   direction = "backward", trace = 1)
  
  # Forward: start at base model, add toward upper
  base_fit <- coxph(f_lower, data = cc_data)
  cat("\n--- Forward selection (", sex_tag, ") ---\n", sep = "")
  fwd_fit  <- step(base_fit, scope = list(lower = f_lower, upper = f_upper),
                   direction = "forward", trace = 1)
  
  list(backward = back_fit, forward = fwd_fit, full = full_fit)
}

sel_men   <- run_stepwise(men_cc,   "MEN")
sel_women <- run_stepwise(women_cc, "WOMEN")

# Compare selected variable sets
vars_selected <- function(fit, tag) {
  v <- names(coef(fit))
  cat(sprintf("  %-25s: %s\n", tag, paste(v, collapse = ", ")))
  v
}
cat("\n--- Selected variables (MEN) ---\n")
v_m_back <- vars_selected(sel_men$backward,   "Backward (MEN)")
v_m_fwd  <- vars_selected(sel_men$forward,    "Forward  (MEN)")
cat("  Agreement: ",
    ifelse(setequal(v_m_back, v_m_fwd), "YES", "NO (see tables)"), "\n")

cat("\n--- Selected variables (WOMEN) ---\n")
v_w_back <- vars_selected(sel_women$backward, "Backward (WOMEN)")
v_w_fwd  <- vars_selected(sel_women$forward,  "Forward  (WOMEN)")
cat("  Agreement: ",
    ifelse(setequal(v_w_back, v_w_fwd), "YES", "NO (see tables)"), "\n")

# Primary final model = backward
cox_men_final   <- sel_men$backward
cox_women_final <- sel_women$backward

cat("\n--- MEN FINAL MODEL (Backward, AIC) ---\n");  print(summary(cox_men_final))
ph_men <- cox.zph(cox_men_final)
cat("--- PH Test (Men) ---\n"); print(ph_men)

cat("\n--- WOMEN FINAL MODEL (Backward, AIC) ---\n"); print(summary(cox_women_final))
ph_women <- cox.zph(cox_women_final)
cat("--- PH Test (Women) ---\n"); print(ph_women)

# Side-by-side comparison table: three selection methods x two sexes
sel_compare <- data.frame(
  Variable = all_vars,
  `Indiv LRT Men`  = ifelse(all_vars %in% results_men$Variable[results_men$LRT_p   < ALPHA_SELECT], "X", ""),
  `Backward Men`   = ifelse(all_vars %in% v_m_back, "X", ""),
  `Forward Men`    = ifelse(all_vars %in% v_m_fwd,  "X", ""),
  `Indiv LRT Women`= ifelse(all_vars %in% results_women$Variable[results_women$LRT_p < ALPHA_SELECT], "X", ""),
  `Backward Women` = ifelse(all_vars %in% v_w_back, "X", ""),
  `Forward Women`  = ifelse(all_vars %in% v_w_fwd,  "X", ""),
  check.names = FALSE, stringsAsFactors = FALSE)
# Mark force-included a priori variables
sel_compare$Variable <- ifelse(sel_compare$Variable %in% sure_vars,
                               paste0(sel_compare$Variable, " (forced)"),
                               sel_compare$Variable)
save_table(sel_compare, "Table5b_SelectionComparison",
           caption = "Table 5b. Variable Selection Comparison: Individual LRT vs Backward vs Forward (AIC). Primary model = Backward.")

save_table(tidy_cox(cox_men_final), "Table6_FinalModel_Men",
           caption = sprintf("Table 6. Final Cox Model - Men (C = %.3f)",
                             summary(cox_men_final)$concordance["C"]))
save_table(tidy_cox(cox_women_final), "Table7_FinalModel_Women",
           caption = sprintf("Table 7. Final Cox Model - Women (C = %.3f)",
                             summary(cox_women_final)$concordance["C"]))

ph_to_df <- function(ph_obj) {
  df <- as.data.frame(ph_obj$table)
  df$Variable <- rownames(df)
  df[df$Variable != "GLOBAL", ] %>%
    mutate(
      Variable = ifelse(Variable %in% names(var_labels),
                        var_labels[Variable], Variable),
      chisq = sprintf("%.3f", chisq),
      p     = sprintf("%.4f", p)
    ) %>%
    select(Variable, chisq, df, p)
}

save_table(ph_to_df(ph_men),   "Table8_PHTest_Men",
           caption = "Table 8. Proportional Hazards Test - Men")
save_table(ph_to_df(ph_women), "Table9_PHTest_Women",
           caption = "Table 9. Proportional Hazards Test - Women")

# Forest Plot
extract_hr <- function(model, sex_label) {
  s    <- summary(model)
  vars <- rownames(s$conf.int)
  data.frame(
    Variable = ifelse(vars %in% names(var_labels), var_labels[vars], vars),
    HR    = s$conf.int[, "exp(coef)"],
    Lower = s$conf.int[, "lower .95"],
    Upper = s$conf.int[, "upper .95"],
    Sex   = sex_label,
    stringsAsFactors = FALSE)
}

hr_all <- rbind(extract_hr(cox_men_final, "Men"),
                extract_hr(cox_women_final, "Women"))

p_forest <- ggplot(hr_all, aes(x = HR, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper, color = Sex),
                  position = position_dodge(width = 0.5), size = 0.6) +
  scale_x_log10() +
  scale_color_manual(values = c("Men" = "#2171b5", "Women" = "#cb181d")) +
  labs(x = "Hazard Ratio (95% CI, log scale)", y = NULL,
       title = "Forest Plot: Final Models by Sex", color = "Sex") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

save_fig(p_forest, "Figure3_Forest_Plot", w = 8, h = 5)

# Time-Varying Covariates
tv_data <- frmgham %>%
  filter(PREVSTRK == 0) %>%
  mutate(sex_label = ifelse(SEX == 1, "Men", "Women"))

# (A) Balanced cohort: participants with a record in all 3 periods
balanced_ids <- tv_data %>%
  group_by(RANDID) %>%
  summarise(np = n_distinct(PERIOD), .groups = "drop") %>%
  filter(np == 3) %>% pull(RANDID)
tv_bal <- tv_data %>% filter(RANDID %in% balanced_ids)
cat("Balanced cohort (present in all 3 periods): n =", length(balanced_ids), "\n")

cont_summary <- tv_bal %>%
  group_by(PERIOD, sex_label) %>%
  summarise(
    `SBP mean (SD)`  = sprintf("%.1f (%.1f)", mean(SYSBP,   na.rm=TRUE), sd(SYSBP,   na.rm=TRUE)),
    `Chol mean (SD)` = sprintf("%.1f (%.1f)", mean(TOTCHOL, na.rm=TRUE), sd(TOTCHOL, na.rm=TRUE)),
    `BMI mean (SD)`  = sprintf("%.1f (%.1f)", mean(BMI,     na.rm=TRUE), sd(BMI,     na.rm=TRUE)),
    .groups = "drop") %>% arrange(sex_label, PERIOD)

bin_summary <- tv_bal %>%
  group_by(PERIOD, sex_label) %>%
  summarise(
    `Diabetes (%)` = sprintf("%.1f", 100 * mean(DIABETES, na.rm=TRUE)),
    `Smoking (%)`  = sprintf("%.1f", 100 * mean(CURSMOKE, na.rm=TRUE)),
    `BP Meds (%)`  = sprintf("%.1f", 100 * mean(BPMEDS,   na.rm=TRUE)),
    .groups = "drop") %>% arrange(sex_label, PERIOD)

tv_combined <- left_join(cont_summary, bin_summary, by = c("PERIOD", "sex_label")) %>%
  rename(Period = PERIOD, Sex = sex_label)
save_table(as.data.frame(tv_combined), "Table10_TimeVarying_Balanced",
           caption = "Table 10. Risk Factor Values by Period (Balanced Cohort, n in all 3 periods)")

# (B) Within-individual change: paired deltas (Period 3 - Period 1)
wide <- tv_bal %>%
  select(RANDID, sex_label, PERIOD, SYSBP, TOTCHOL, BMI, DIABETES, CURSMOKE, BPMEDS) %>%
  pivot_wider(names_from = PERIOD,
              values_from = c(SYSBP, TOTCHOL, BMI, DIABETES, CURSMOKE, BPMEDS),
              names_sep = "_P")

within_change <- wide %>%
  mutate(dSBP = SYSBP_P3 - SYSBP_P1,
         dChol = TOTCHOL_P3 - TOTCHOL_P1,
         dBMI  = BMI_P3    - BMI_P1,
         incident_DM     = ifelse(DIABETES_P1 == 0 & DIABETES_P3 == 1, 1, 0),
         quit_smoking    = ifelse(CURSMOKE_P1 == 1 & CURSMOKE_P3 == 0, 1, 0),
         started_BPmeds  = ifelse(BPMEDS_P1   == 0 & BPMEDS_P3   == 1, 1, 0)) %>%
  group_by(Sex = sex_label) %>%
  summarise(
    N = n(),
    `Delta SBP mean (SD)`  = sprintf("%.1f (%.1f)", mean(dSBP,  na.rm=TRUE), sd(dSBP,  na.rm=TRUE)),
    `Delta Chol mean (SD)` = sprintf("%.1f (%.1f)", mean(dChol, na.rm=TRUE), sd(dChol, na.rm=TRUE)),
    `Delta BMI mean (SD)`  = sprintf("%.2f (%.2f)", mean(dBMI,  na.rm=TRUE), sd(dBMI,  na.rm=TRUE)),
    `Incident diabetes (%)` = sprintf("%.1f", 100*mean(incident_DM,    na.rm=TRUE)),
    `Quit smoking (%)`      = sprintf("%.1f", 100*mean(quit_smoking,   na.rm=TRUE)),
    `Started BP meds (%)`   = sprintf("%.1f", 100*mean(started_BPmeds, na.rm=TRUE)),
    .groups = "drop")
save_table(as.data.frame(within_change), "Table10b_WithinChange",
           caption = "Table 10b. Within-Individual Change, Period 1 to Period 3 (Balanced Cohort)")


# Competing Risk Summary
cr_summary <- baseline %>%
  group_by(Sex = sex_label) %>%
  summarise(
    N = n(),
    `Stroke n (%)`           = sprintf("%d (%.1f%%)", sum(stroke_10yr),      100 * mean(stroke_10yr)),
    `Death w/o stroke n (%)` = sprintf("%d (%.1f%%)", sum(death_no_stroke),  100 * mean(death_no_stroke)),
    `Censored n (%)`         = sprintf("%d (%.1f%%)",
                                       N - sum(stroke_10yr) - sum(death_no_stroke),
                                       100 * (N - sum(stroke_10yr) - sum(death_no_stroke)) / N),
    .groups = "drop")

print(as.data.frame(cr_summary))
save_table(as.data.frame(cr_summary), "Table11_CompetingRisk",
           caption = "Table 11. Competing Risk: 10-Year Outcomes by Sex")

# Example Predicted Probabilities
model_vars <- names(cox_men_final$coefficients)

avg_vals  <- list(AGE = 50, DIABETES = 0, SYSBP = 132, CURSMOKE = 0,
                  PREVCHD = 0, BPMEDS = 0, TOTCHOL = 234, BMI = 26)
high_vals <- list(AGE = 70, DIABETES = 1, SYSBP = 160, CURSMOKE = 1,
                  PREVCHD = 1, BPMEDS = 1, TOTCHOL = 260, BMI = 30)

surv_avg  <- survfit(cox_men_final, newdata = as.data.frame(avg_vals[model_vars]))
surv_high <- survfit(cox_men_final, newdata = as.data.frame(high_vals[model_vars]))