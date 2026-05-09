# =============================================================================
# Project 4: Variable Selection Simulation Study
# Course:    BIOS 6624
# Author:    Yunji
#
# Place this file in Project_4/Code/ and run with Project_4/ as the working
# directory (e.g. via RStudio's "Set as Working Directory" on the project root,
# or Session > Set Working Directory > To Source File Location if the Rmd is
# open). Do NOT hard-code an absolute path; relative paths are used throughout.
#
# This is a pure simulation study — no external data file is required.
# All data are generated programmatically within this script.
#
# Output directories (created automatically if absent):
#   Results/Figures/      PNG figures
#   Results/Tables/       PNG + .tex tables
#   Results/sessionInfo.txt   R session information for reproducibility
#   Results/simulation_results.csv   Full condition-level summary
#   Results/simulation_results.rds   Full results as R object
# =============================================================================

# --- Reproducibility ---------------------------------------------------------
# Seed 2026 is set as the very first executable line so every downstream call
# to any RNG (data generation, cv.glmnet fold assignment, etc.) is determined.
set.seed(2026)

# --- Packages ----------------------------------------------------------------
# Data generation: the assignment suggests hdrm::gen_data() (github.com/
# drjohnson1/hdrm). We implement the same exchangeable-Sigma MVN directly
# using MASS::mvrnorm(), which requires no additional GitHub installation and
# produces identical draws for this covariance structure.
#
# Install missing CRAN packages before loading:
#   install.packages(c("MASS","glmnet","ggplot2","dplyr",
#                      "reshape2","grid","gridExtra","xtable"))
library(MASS)       # mvrnorm(): multivariate normal data generation
library(glmnet)     # cv.glmnet() / glmnet(): LASSO and Elastic Net
library(ggplot2)    # all figures
library(dplyr)      # data wrangling
library(reshape2)   # melt(): heatmap reshape
library(grid)       # grid.newpage(), viewport(), grid.text()
library(gridExtra)  # tableGrob(): data frame → table graphic
library(xtable)     # xtable(): LaTeX table export

# --- Output directories ------------------------------------------------------
for (d in c("Results/Figures", "Results/Tables"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
fig_dir <- "Results/Figures/"
tab_dir <- "Results/Tables/"

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================
# B = 1 000 replications.  Justification: for a binomial proportion p the
# Monte Carlo standard error is sqrt(p(1-p)/B) <= sqrt(0.25/1000) ~ 0.016.
# This gives ~+/-2 percentage-point precision on every TPR/FPR cell, which is
# sufficient to detect practically meaningful differences between methods.
B         <- 1000
P         <- 20          # total predictors
P_TRUE    <- 5           # signal predictors: X1-X5
P_NULL    <- 15          # null predictors:   X6-X20
SIGMA     <- 1           # error SD (sigma^2 = 1)

TRUE_VARS <- paste0("X", 1:P_TRUE)
NULL_VARS <- paste0("X", (P_TRUE + 1):P)
ALL_VARS  <- paste0("X", 1:P)

# True coefficients as specified in the project description:
# X1-X5 get 0.5/3, 1/3, 1.5/3, 2/3, 2.5/3; X6-X20 get 0
BETA_TRUE <- setNames(c(1/6, 2/6, 3/6, 4/6, 5/6, rep(0, P_NULL)), ALL_VARS)

# Sample sizes and correlation levels: both specified by class decisions
N_VALUES   <- c(250, 500)
RHO_VALUES <- c(0, 0.35, 0.70)   # exchangeable (compound symmetry); class decision

# Backward elimination threshold (class decision)
ALPHA_BACK <- 0.10

# glmnet cross-validation folds (class decision)
CV_FOLDS   <- 10

# Elastic Net mixing parameter: alpha = 0.5 balances L1 (LASSO) and L2 (ridge)
# equally; chosen as a symmetric default since no prior information favours
# either penalty.  alpha = 1 → pure LASSO; alpha = 0 → pure ridge.
ALPHA_ENET <- 0.5

# Full-model formula (used by stepwise methods)
full_formula <- as.formula(paste("Y ~", paste(ALL_VARS, collapse = " + ")))

# Method labels (order matches columns 1-7 throughout)
METHOD_LABELS <- c(
  "Backward (p-val)",
  "AIC",
  "BIC",
  "LASSO (lambda.min)",
  "LASSO (lambda.1se)",
  "Elastic Net (lambda.min)",
  "Elastic Net (lambda.1se)"
)

# =============================================================================
# FITTING FUNCTIONS  (stepwise only; penalized fitted inline in the loop)
# =============================================================================

# --- Backward elimination by p-value -----------------------------------------
# The project description specifies F-tests; for single-df OLS tests F = t^2,
# so t-test p-values from summary.lm() are numerically identical.
# deviation from step() default: we implement this manually (not via step())
# so the removal threshold is exact p >= ALPHA_BACK rather than AIC/BIC.
# Robustness: if ALL predictors are eliminated (intercept-only model), returns
# lm(Y ~ 1) rather than erroring on an empty formula.
fit_backward <- function(df) {
  remaining <- ALL_VARS
  repeat {
    if (length(remaining) == 0) break          # guard: nothing left
    form  <- as.formula(paste("Y ~", paste(remaining, collapse = " + ")))
    fit   <- lm(form, data = df)
    coef_tbl <- summary(fit)$coefficients
    # drop intercept row; guard against intercept-only (zero non-intercept rows)
    pvals <- coef_tbl[rownames(coef_tbl) != "(Intercept)", "Pr(>|t|)", drop = TRUE]
    if (length(pvals) == 0 || max(pvals) < ALPHA_BACK) break
    remaining <- remaining[remaining != names(which.max(pvals))]
  }
  if (length(remaining) == 0)
    lm(Y ~ 1, data = df)           # intercept-only: no predictors survived
  else
    lm(as.formula(paste("Y ~", paste(remaining, collapse = " + "))), data = df)
}

# --- Stepwise AIC: step() with direction="backward", k=2 --------------------
# deviation from step() default: direction="backward" (default is "both")
fit_aic <- function(df)
  step(lm(full_formula, data = df), direction = "backward", trace = 0, k = 2)

# --- Stepwise BIC: step() with direction="backward", k=log(N) ---------------
# deviation from step() default: direction="backward"; k=log(N) instead of k=2
fit_bic <- function(df)
  step(lm(full_formula, data = df), direction = "backward", trace = 0,
       k = log(nrow(df)))

# =============================================================================
# EXTRACTION FUNCTIONS
# Each returns list(vars = character vector of selected names,
#                   coefs = named numeric vector of estimated coefficients)
# =============================================================================

extract_lm <- function(mod) {
  cf   <- coef(mod)
  vars <- setdiff(names(cf), "(Intercept)")
  list(vars = vars, coefs = cf[vars])
}

extract_glmnet <- function(fit_obj) {
  cm  <- as.matrix(coef(fit_obj$fit))[-1, 1]   # fit_obj is list(fit = glmnet(...)); extract $fit first
  sel <- names(cm)[cm != 0]
  list(vars = sel, coefs = cm[sel])
}

# =============================================================================
# EVALUATE_REPLICATION
# Returns per-replication TPR, FPR, bias vector, coverage vector, and the
# number of signal variables that contributed to coverage (the denominator).
#
# NON-SELECTION HANDLING
# ── Bias: computed UNCONDITIONALLY.
#    Signal variables not selected receive estimated coefficient = 0, so their
#    contribution to bias is (0 − beta_j) = −beta_j.  This conservative choice
#    captures the full cost of missing a true signal.
#    For penalized methods the penalized coefficient (not an OLS refit) is used
#    as the estimate for selected variables, reflecting actual shrinkage.
#    For stepwise methods the OLS coefficient from the fitted model is used.
#
# ── CI coverage: computed CONDITIONALLY on selection.
#    A 95% CI is computed via confint(lm()) only for variables that appear in
#    the final selected model.  Variables not selected contribute NA and are
#    excluded from the coverage average.  Justification: a confidence interval
#    only exists for a variable in the analyst's final model; coverage for a
#    variable that was never estimated is not a meaningful quantity.
#    The denominator (number of signal variables selected) is also returned so
#    the reader can assess how thinly the coverage average is spread.
# =============================================================================
evaluate_replication <- function(sel_vars, sel_coefs, df) {
  
  # TPR: proportion of true signal vars selected
  TPR <- sum(TRUE_VARS %in% sel_vars) / P_TRUE
  
  # FPR: proportion of null vars selected
  FPR <- sum(NULL_VARS %in% sel_vars) / P_NULL
  
  # Bias (unconditional: non-selected = 0 estimate)
  coef_est <- setNames(rep(0, P_TRUE), TRUE_VARS)
  for (v in intersect(TRUE_VARS, sel_vars)) coef_est[v] <- sel_coefs[v]
  bias <- coef_est - BETA_TRUE[TRUE_VARS]
  
  # CI coverage (conditional on selection via OLS refit)
  covered <- setNames(rep(NA_real_, P_TRUE), TRUE_VARS)
  if (length(sel_vars) > 0) {
    tryCatch({
      refit <- lm(as.formula(paste("Y ~", paste(sel_vars, collapse = " + "))),
                  data = df)
      ci <- confint(refit, level = 0.95)
      for (v in TRUE_VARS)
        if (v %in% rownames(ci))
          covered[v] <- as.numeric(ci[v,1] <= BETA_TRUE[v] &
                                     BETA_TRUE[v] <= ci[v,2])
    }, error = function(e) NULL)
  }
  
  # Denominator for coverage average: how many signal vars were selected?
  n_sel_signal <- sum(TRUE_VARS %in% sel_vars)
  
  list(TPR = TPR, FPR = FPR, bias = bias,
       covered = covered, n_sel_signal = n_sel_signal)
}

# =============================================================================
# DATA GENERATION
# Uses MASS::mvrnorm() with compound-symmetry Sigma (Sigma_jk = rho, j != k).
# =============================================================================
generate_data <- function(N, rho) {
  Sigma       <- matrix(rho, P, P); diag(Sigma) <- 1
  X           <- mvrnorm(N, mu = rep(0, P), Sigma = Sigma)
  colnames(X) <- ALL_VARS
  Y           <- as.numeric(X %*% BETA_TRUE[ALL_VARS] + rnorm(N, 0, SIGMA))
  df          <- as.data.frame(X); df$Y <- Y; df
}

# =============================================================================
# MAIN SIMULATION LOOP
# =============================================================================
cat("Simulation started  |  B =", B, "| seed = 2026\n")

sim_rows <- list()

for (N in N_VALUES) {
  for (rho in RHO_VALUES) {
    cat(sprintf("\n  Condition: N = %d, rho = %.2f\n", N, rho))
    
    tpr_mat     <- matrix(NA_real_, B, 7)
    fpr_mat     <- matrix(NA_real_, B, 7)
    bias_mat    <- matrix(NA_real_, B, 7)   # mean unconditional bias over X1-X5
    cov_mat     <- matrix(NA_real_, B, 7)   # mean conditional coverage over X1-X5
    nsel_mat    <- matrix(NA_real_, B, 7)   # mean # signal vars selected (coverage denom)
    sel_arr     <- array(0L,        dim = c(B, 7, P_TRUE))  # 1 = selected
    bias_pv_arr <- array(NA_real_,  dim = c(B, 7, P_TRUE))  # per-variable bias
    
    for (b in 1:B) {
      if (b %% 200 == 0) cat(sprintf("    rep %d / %d\n", b, B))
      
      df <- generate_data(N, rho)
      Xm <- as.matrix(df[, ALL_VARS])
      Yv <- df$Y
      
      # ── Stepwise methods ──────────────────────────────────────────────────
      fits_sw <- list(fit_backward(df), fit_aic(df), fit_bic(df))
      
      # ── Penalized: ONE cv.glmnet call per penalty type ─────────────────
      # Using a single CV fit for both lambda choices ensures the two
      # thresholds share the same fold assignments (no randomness mismatch).
      cv_lasso <- cv.glmnet(Xm, Yv, alpha = 1,         nfolds = CV_FOLDS,
                            standardize = TRUE)
      cv_enet  <- cv.glmnet(Xm, Yv, alpha = ALPHA_ENET, nfolds = CV_FOLDS,
                            standardize = TRUE)
      
      pen_fits <- list(
        list(fit = glmnet(Xm, Yv, alpha = 1,         lambda = cv_lasso$lambda.min, standardize = TRUE)),
        list(fit = glmnet(Xm, Yv, alpha = 1,         lambda = cv_lasso$lambda.1se, standardize = TRUE)),
        list(fit = glmnet(Xm, Yv, alpha = ALPHA_ENET, lambda = cv_enet$lambda.min,  standardize = TRUE)),
        list(fit = glmnet(Xm, Yv, alpha = ALPHA_ENET, lambda = cv_enet$lambda.1se,  standardize = TRUE))
      )
      
      # ── Extract ──────────────────────────────────────────────────────────
      extracted <- c(
        lapply(fits_sw,   extract_lm),
        lapply(pen_fits,  extract_glmnet)
      )
      
      # ── Evaluate ─────────────────────────────────────────────────────────
      for (m in seq_len(7)) {
        ev <- evaluate_replication(extracted[[m]]$vars, extracted[[m]]$coefs, df)
        tpr_mat[b, m]      <- ev$TPR
        fpr_mat[b, m]      <- ev$FPR
        bias_mat[b, m]     <- mean(ev$bias)
        cov_mat[b, m]      <- mean(ev$covered, na.rm = TRUE)
        nsel_mat[b, m]     <- ev$n_sel_signal
        sel_arr[b, m, ]    <- as.integer(TRUE_VARS %in% extracted[[m]]$vars)
        bias_pv_arr[b, m, ]<- ev$bias
      }
    } # end replications
    
    # ── Condition-level summary ─────────────────────────────────────────────
    for (m in seq_len(7)) {
      sim_rows[[length(sim_rows) + 1]] <- data.frame(
        N             = N,
        rho           = rho,
        Method        = METHOD_LABELS[m],
        Mean_TPR      = mean(tpr_mat[, m]),
        SD_TPR        = sd(tpr_mat[, m]),
        Mean_FPR      = mean(fpr_mat[, m]),
        SD_FPR        = sd(fpr_mat[, m]),
        Mean_Bias     = mean(bias_mat[, m]),
        Mean_Cov      = mean(cov_mat[, m], na.rm = TRUE),
        # Average number of signal vars contributing to coverage each replication
        Mean_N_Sel_Signal = mean(nsel_mat[, m]),
        SelRate_X1    = mean(sel_arr[, m, 1]),
        SelRate_X2    = mean(sel_arr[, m, 2]),
        SelRate_X3    = mean(sel_arr[, m, 3]),
        SelRate_X4    = mean(sel_arr[, m, 4]),
        SelRate_X5    = mean(sel_arr[, m, 5]),
        Bias_X1       = mean(bias_pv_arr[, m, 1]),
        Bias_X2       = mean(bias_pv_arr[, m, 2]),
        Bias_X3       = mean(bias_pv_arr[, m, 3]),
        Bias_X4       = mean(bias_pv_arr[, m, 4]),
        Bias_X5       = mean(bias_pv_arr[, m, 5]),
        stringsAsFactors = FALSE
      )
    }
  } # end rho
} # end N

results        <- do.call(rbind, sim_rows)
results$Method <- factor(results$Method, levels = METHOD_LABELS)
cat("\nSimulation complete.\n")

# =============================================================================
# SAVE RESULTS
# =============================================================================
# Full results as CSV (human-readable, shareable)
write.csv(results, file = "Results/simulation_results.csv", row.names = FALSE)
cat("Saved: Results/simulation_results.csv\n")

# Full results as R object (preserves factor levels etc.)
saveRDS(results, file = "Results/simulation_results.rds")
cat("Saved: Results/simulation_results.rds\n")

# Session information for exact reproducibility (package versions, R version)
writeLines(capture.output(sessionInfo()), "Results/sessionInfo.txt")
cat("Saved: Results/sessionInfo.txt\n")

# =============================================================================
# HELPER: render a data frame as a PNG table + .tex file
# =============================================================================
save_table_png <- function(df, filename, caption, w = 9, h = 4) {
  png(paste0(tab_dir, filename, ".png"),
      width = w, height = h, units = "in", res = 200, bg = "white")
  grid.newpage()
  grid.text(caption, x = 0.5, y = 0.97,
            gp = gpar(fontsize = 11, fontface = "bold"))
  tt <- ttheme_minimal(
    core    = list(fg_params = list(fontsize = 8), padding = unit(c(3, 2.5), "mm")),
    colhead = list(fg_params = list(fontsize = 8, fontface = "bold"),
                   padding   = unit(c(3, 2.5), "mm")))
  tg <- tableGrob(df, rows = NULL, theme = tt)
  vp <- viewport(x = 0.5, y = 0.44, width = 0.98,
                 height = 0.86, just = c("center", "center"))
  pushViewport(vp); grid.draw(tg); popViewport()
  dev.off()
  cat("Saved:", paste0(tab_dir, filename, ".png\n"))
  
  print(xtable(df, caption = caption),
        file = paste0(tab_dir, filename, ".tex"),
        include.rownames = FALSE, booktabs = TRUE,
        caption.placement = "top", sanitize.text.function = identity)
}

# =============================================================================
# TABLE 1: Simulation Parameters
# =============================================================================
t1 <- data.frame(
  Parameter = c(
    "Replications (B)", "Total predictors (p)",
    "Signal predictors", "Null predictors",
    "True coefficients (X1-X5)", "Null coefficients (X6-X20)",
    "Error SD (sigma)", "Sample sizes (N)",
    "Correlation structure", "Correlation values (rho) [class decision]",
    "Backward selection threshold [class decision]",
    "LASSO lambda choices [class decision]",
    "Elastic net mixing alpha [class decision]",
    "Elastic net lambda choices [class decision]",
    "Cross-validation folds [class decision]"
  ),
  Value = c(
    "1,000", "20", "X1-X5", "X6-X20",
    "1/6, 2/6, 3/6, 4/6, 5/6", "0", "1",
    "250, 500", "Exchangeable (compound symmetry)", "0, 0.35, 0.70",
    "p-value < 0.10", "lambda.min, lambda.1se", "0.50",
    "lambda.min, lambda.1se", "10"
  ),
  stringsAsFactors = FALSE
)
save_table_png(t1, "Table1_SimulationParameters",
               "Table 1. Simulation Parameter Settings", w = 9, h = 5.5)

# =============================================================================
# TABLE 2: Mean TPR and FPR
# =============================================================================
t2 <- results %>%
  mutate(`Mean TPR (SD)` = sprintf("%.3f (%.3f)", Mean_TPR, SD_TPR),
         `Mean FPR (SD)` = sprintf("%.3f (%.3f)", Mean_FPR, SD_FPR)) %>%
  select(N, rho, Method, `Mean TPR (SD)`, `Mean FPR (SD)`) %>%
  arrange(N, rho, Method)
save_table_png(t2, "Table2_TPR_FPR",
               "Table 2. Mean TPR and FPR by Method and Condition (SD in parentheses; B = 1000)",
               w = 9, h = 14)

# =============================================================================
# TABLE 3: Bias and Coverage (with coverage denominator)
# =============================================================================
t3 <- results %>%
  mutate(
    `Mean Bias`         = round(Mean_Bias,         4),
    `95% CI Cov`        = round(Mean_Cov,           3),
    `Avg Sigs Selected` = round(Mean_N_Sel_Signal,  2)
  ) %>%
  select(N, rho, Method, `Mean Bias`, `95% CI Cov`, `Avg Sigs Selected`) %>%
  arrange(N, rho, Method)
save_table_png(t3, "Table3_Bias_Coverage",
               "Table 3. Mean Bias (unconditional), 95% CI Coverage (conditional), and Average Signal Variables Selected",
               w = 10, h = 14)

# =============================================================================
# TABLE 4: Confusion Matrix (N=250, rho=0.35)
# =============================================================================
t4 <- results %>%
  filter(N == 250, rho == 0.35) %>%
  mutate(`TP Rate` = round(Mean_TPR,       3),
         `FN Rate` = round(1 - Mean_TPR,   3),
         `FP Rate` = round(Mean_FPR,       3),
         `TN Rate` = round(1 - Mean_FPR,   3)) %>%
  select(Method, `TP Rate`, `FN Rate`, `FP Rate`, `TN Rate`)
save_table_png(t4, "Table4_ConfusionMatrix_N250_rho035",
               "Table 4. Confusion Matrix Rates: N=250, rho=0.35 (TP/FN over X1-X5; FP/TN over X6-X20)",
               w = 8, h = 3)

# =============================================================================
# TABLE 5: Per-Variable Selection Rates (N=250, rho=0 and 0.70)
# =============================================================================
t5 <- results %>%
  filter(N == 250, rho %in% c(0, 0.70)) %>%
  select(N, rho, Method,
         `X1 (b=0.167)` = SelRate_X1,
         `X2 (b=0.333)` = SelRate_X2,
         `X3 (b=0.500)` = SelRate_X3,
         `X4 (b=0.667)` = SelRate_X4,
         `X5 (b=0.833)` = SelRate_X5) %>%
  arrange(rho, Method)
save_table_png(t5, "Table5_PerVariable_SelectionRates",
               "Table 5. Per-Variable Selection Rate for N=250, rho=0 and rho=0.70",
               w = 10, h = 4.5)

# =============================================================================
# FIGURE AESTHETICS
# =============================================================================
method_colors <- c(
  "Backward (p-val)"         = "#1f77b4",
  "AIC"                      = "#ff7f0e",
  "BIC"                      = "#2ca02c",
  "LASSO (lambda.min)"       = "#d62728",
  "LASSO (lambda.1se)"       = "#aec7e8",
  "Elastic Net (lambda.min)" = "#ff69b4",
  "Elastic Net (lambda.1se)" = "#7f7f7f"
)
method_shapes <- c(
  "Backward (p-val)"         = 16,
  "AIC"                      = 17,
  "BIC"                      = 15,
  "LASSO (lambda.min)"       = 18,
  "LASSO (lambda.1se)"       =  3,
  "Elastic Net (lambda.min)" =  8,
  "Elastic Net (lambda.1se)" =  4
)
results$N_label   <- paste0("N = ", results$N)
results$rho_label <- paste0("rho = ", results$rho)
base_theme <- theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text  = element_text(face = "bold"),
        plot.title  = element_text(face = "bold"))

# =============================================================================
# FIGURE 1: Mean TPR
# =============================================================================
p1 <- ggplot(results, aes(x = Method, y = Mean_TPR, fill = Method)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "orange", linewidth = 0.8) +
  facet_grid(N_label ~ rho_label) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(breaks = seq(0.70, 1.0, 0.05)) +
  coord_cartesian(ylim = c(0.70, 1.01)) +
  labs(title    = "Figure 1. Mean True Positive Rate",
       subtitle = "Proportion of replications X1-X5 was selected (orange dashed = perfect TPR = 1.0)",
       x = "Variable Selection Method", y = "Mean True Positive Rate") +
  base_theme
ggsave(paste0(fig_dir, "Figure1_TruePositiveRate.png"),
       p1, width = 12, height = 7, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 2: Mean FPR
# (subtitle fixed: FPR is avg proportion of individual null vars selected,
#  not whether *any* null var was selected)
# =============================================================================
p2 <- ggplot(results, aes(x = Method, y = Mean_FPR, fill = Method)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "orange", linewidth = 0.8) +
  facet_grid(N_label ~ rho_label) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 0.65), breaks = seq(0, 0.6, 0.1)) +
  labs(title    = "Figure 2. Mean False Positive Rate",
       subtitle = "Average proportion of null variables (X6-X20) selected per replication; orange dashed = 0.10",
       x = "Variable Selection Method", y = "Mean False Positive Rate") +
  base_theme
ggsave(paste0(fig_dir, "Figure2_FalsePositiveRate.png"),
       p2, width = 12, height = 7, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 3: Mean Coefficient Bias (unconditional)
# =============================================================================
p3 <- ggplot(results, aes(x = Method, y = Mean_Bias, fill = Method)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  facet_grid(N_label ~ rho_label) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(title    = "Figure 3. Mean Coefficient Bias for Signal Variables X1-X5",
       subtitle = "Unconditional average of (beta_hat - beta_true); non-selected variables contribute 0 - beta_true",
       x = "Variable Selection Method", y = "Mean Bias") +
  base_theme
ggsave(paste0(fig_dir, "Figure3_CoefficientBias.png"),
       p3, width = 12, height = 7, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 4: CI Coverage (conditional)
# =============================================================================
p4 <- ggplot(results, aes(x = Method, y = Mean_Cov, fill = Method)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "steelblue", linewidth = 0.8) +
  facet_grid(N_label ~ rho_label) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(breaks = seq(0.88, 1.0, 0.02)) +
  coord_cartesian(ylim = c(0.88, 1.01)) +
  labs(title    = "Figure 4. Empirical 95% CI Coverage for Signal Variables X1-X5",
       subtitle = "Conditional on selection; post-selection OLS refit; blue dashed = nominal 95%",
       x = "Variable Selection Method", y = "Coverage Probability") +
  base_theme
ggsave(paste0(fig_dir, "Figure4_CICoverage.png"),
       p4, width = 12, height = 7, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 5: TPR vs FPR scatter (ROC-style)
# =============================================================================
p5 <- ggplot(results, aes(x = Mean_FPR, y = Mean_TPR, color = Method, shape = Method)) +
  geom_point(size = 3) +
  facet_grid(N_label ~ rho_label) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  labs(title    = "Figure 5. TPR vs FPR Scatter by Method (ROC-Style)",
       subtitle = "Upper-left corner = ideal (high sensitivity, low false discovery rate)",
       x = "Mean False Positive Rate", y = "Mean True Positive Rate",
       color = "Method", shape = "Method") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
ggsave(paste0(fig_dir, "Figure5_TPR_vs_FPR_scatter.png"),
       p5, width = 12, height = 8, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 6: Per-Variable Selection Rate Heatmap (N=250, rho=0.35)
# =============================================================================
hm_df <- results %>%
  filter(N == 250, rho == 0.35) %>%
  select(Method, SelRate_X1:SelRate_X5) %>%
  setNames(c("Method", "X1\n(b=0.167)", "X2\n(b=0.333)",
             "X3\n(b=0.500)", "X4\n(b=0.667)", "X5\n(b=0.833)")) %>%
  melt(id.vars = "Method", variable.name = "Variable", value.name = "SelRate")

p6 <- ggplot(hm_df, aes(x = Method, y = Variable, fill = SelRate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", SelRate)),
            color = ifelse(hm_df$SelRate > 0.5, "white", "black"),
            fontface = "bold", size = 3.5) +
  scale_fill_gradient2(low = "white", mid = "#56B4E9", high = "#00008B",
                       midpoint = 0.5, limits = c(0, 1), name = "Sel. Rate") +
  scale_x_discrete(limits = METHOD_LABELS) +
  labs(title    = "Figure 6. Per-Variable Selection Rate Heatmap (N=250, rho=0.35)",
       subtitle = "Fraction of 1000 replications each signal variable was retained",
       x = "Method", y = "Signal Variable (True beta)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold"))
ggsave(paste0(fig_dir, "Figure6_SelectionRate_Heatmap_N250_rho035.png"),
       p6, width = 9, height = 5, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 7: Per-Variable Coefficient Bias (N=500, rho=0; unconditional)
# =============================================================================
pv_long <- results %>%
  filter(N == 500, rho == 0) %>%
  select(Method, Bias_X1:Bias_X5) %>%
  setNames(c("Method", "X1 (b=0.167)", "X2 (b=0.333)",
             "X3 (b=0.500)", "X4 (b=0.667)", "X5 (b=0.833)")) %>%
  melt(id.vars = "Method", variable.name = "Variable", value.name = "MeanBias")

p7 <- ggplot(pv_long, aes(x = Variable, y = MeanBias,
                          color = Method, group = Method, shape = Method)) +
  geom_line() + geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  labs(title    = "Figure 7. Per-Variable Coefficient Bias (N=500, rho=0)",
       subtitle = "Unconditional bias (non-selected = 0); stepwise near-zero; penalized methods show shrinkage",
       x = "Signal Variable (True beta)", y = "Mean Bias",
       color = "Method", shape = "Method") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))
ggsave(paste0(fig_dir, "Figure7_PerVariable_Bias_N500_rho0.png"),
       p7, width = 9, height = 5, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 8: Effect of Correlation on TPR
# =============================================================================
p8 <- ggplot(results, aes(x = rho, y = Mean_TPR,
                          color = Method, group = Method, shape = Method)) +
  geom_line() + geom_point(size = 2.5) +
  facet_wrap(~ N_label) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_x_continuous(breaks = c(0, 0.35, 0.70)) +
  scale_y_continuous(limits = c(0.75, 1.01), breaks = seq(0.75, 1.0, 0.05)) +
  labs(title    = "Figure 8. Effect of Predictor Correlation on True Positive Rate",
       subtitle = "Elastic Net degrades least; BIC most sensitive to correlation",
       x = "Correlation (rho)", y = "Mean True Positive Rate",
       color = "Method", shape = "Method") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
ggsave(paste0(fig_dir, "Figure8_Correlation_Effect_on_TPR.png"),
       p8, width = 10, height = 5, dpi = 200, bg = "white")

cat("\nAll figures and tables saved.\n")