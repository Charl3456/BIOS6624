# =============================================================================
# Project 4: Investigation of Variable Selection Algorithms
# Course:    BIOS 6624 — Advanced Statistics
# File:      Code/Project4_VarSelection_Simulation.R
#
# Description:
#   Compares 7 variable selection approaches (Backward/p-val, AIC, BIC,
#   LASSO λ.min, LASSO λ.1se, Elastic Net λ.min, Elastic Net λ.1se)
#   across 6 simulation conditions (N=250/500 × ρ=0/0.35/0.70).
#
#   Data generation: 20 predictors, 5 truly important (β = 0.5/3 … 2.5/3),
#   15 null (β = 0). Exchangeable correlation among all Xs (class decision).
#
#   Metrics: True Positive Rate, False Positive Rate, coefficient bias,
#            empirical 95% CI coverage (naive post-selection refit).
#
# Output structure (relative to Code/ folder):
#   ../Results/Figures/   ← all .png figures
#   ../Results/Tables/    ← all .tex and .png tables
#   ../Results/           ← simulation_results.rds
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 0.  PACKAGE MANAGEMENT
# ─────────────────────────────────────────────────────────────────────────────

required_pkgs <- c("MASS", "glmnet", "xtable", "ggplot2", "gridExtra",
                   "grid", "dplyr", "reshape2", "RColorBrewer")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# Attempt to install hdrm from GitHub (instructor-specified data-generation pkg)
# If available, hdrm::genXY() can replace generate_data() below.
if (!requireNamespace("hdrm", quietly = TRUE)) {
  message("hdrm not found — attempting install from GitHub (requires devtools).")
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools", repos = "https://cloud.r-project.org")
  tryCatch(
    devtools::install_github("pbreheny/hdrm"),
    error = function(e)
      message("hdrm unavailable; using MASS::mvrnorm for data generation.")
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# 1.  DIRECTORY SETUP
# ─────────────────────────────────────────────────────────────────────────────

fig_dir <- "../Results/Figures/"
tab_dir <- "../Results/Tables/"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)


# ─────────────────────────────────────────────────────────────────────────────
# 2.  SIMULATION PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────

set.seed(2026)

N_SIM    <- 1000          # replications per condition
P        <- 20            # total predictors
P_TRUE   <- 5             # truly important predictors (X1–X5)
P_NULL   <- P - P_TRUE    # null predictors (X6–X20)

N_VALS   <- c(250, 500)         # sample sizes (Cases 1 and 2)
RHO_VALS <- c(0, 0.35, 0.70)   # exchangeable correlation (class decision)
#   rho=0    → Cases 1a / 2a (independent predictors)
#   rho=0.35 → correlated sub-case (moderate)
#   rho=0.70 → correlated sub-case (strong)

# True regression coefficients:
#   X1–X5: 0.5/3, 1/3, 1.5/3, 2.0/3, 2.5/3
#   X6–X20: 0 (null)
BETA_TRUE <- c(0.5/3, 1/3, 1.5/3, 2.0/3, 2.5/3, rep(0, P_NULL))

SIGMA_E      <- 1.0    # error SD (σ)
ALPHA_BACK   <- 0.10   # backward selection p-value threshold
ALPHA_ENET   <- 0.50   # elastic net mixing parameter (0=ridge, 1=lasso)

TRUE_VARS <- paste0("X", seq_len(P_TRUE))
NULL_VARS <- paste0("X", (P_TRUE + 1):P)

METHOD_NAMES <- c(
  "Backward_pval",
  "AIC",
  "BIC",
  "LASSO_min",
  "LASSO_1se",
  "ENet_min",
  "ENet_1se"
)

METHOD_LABELS <- c(
  "Backward (p-val)",
  "AIC",
  "BIC",
  "LASSO (λ.min)",
  "LASSO (λ.1se)",
  "Elastic Net (λ.min)",
  "Elastic Net (λ.1se)"
)


# ─────────────────────────────────────────────────────────────────────────────
# 3.  DATA GENERATION
# ─────────────────────────────────────────────────────────────────────────────

# Exchangeable (compound symmetry) covariance matrix for predictors
make_exch_cov <- function(p, rho) {
  S        <- matrix(rho, nrow = p, ncol = p)
  diag(S)  <- 1
  S
}

# Generate one simulated dataset
# If hdrm is installed, you can replace the body with:
#   hdrm::genXY(n, p, beta, rho, family="gaussian", ...)
generate_data <- function(n, p, rho, beta, sigma_e = 1) {
  Sigma <- make_exch_cov(p, rho)
  X     <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", seq_len(p))
  Y  <- drop(X %*% beta) + rnorm(n, 0, sigma_e)
  data.frame(Y = Y, X)
}


# ─────────────────────────────────────────────────────────────────────────────
# 4.  MODEL FITTING FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

full_formula <- as.formula(
  paste("Y ~", paste(paste0("X", seq_len(P)), collapse = " + "))
)

# ── 4a.  Backward selection by p-value ──────────────────────────────────────
fit_backward_pval <- function(df, alpha = ALPHA_BACK) {
  mod <- lm(full_formula, data = df)
  repeat {
    cf    <- summary(mod)$coefficients
    pvals <- cf[-1, "Pr(>|t|)", drop = FALSE]   # exclude intercept
    if (nrow(pvals) == 0 || max(pvals) <= alpha) break
    drop_var <- rownames(pvals)[which.max(pvals)]
    current  <- setdiff(attr(terms(mod), "term.labels"), drop_var)
    if (length(current) == 0) { mod <- lm(Y ~ 1, data = df); break }
    mod <- lm(as.formula(paste("Y ~", paste(current, collapse = " + "))),
              data = df)
  }
  mod
}

# ── 4b.  AIC (backward step) ────────────────────────────────────────────────
fit_aic <- function(df) {
  step(lm(full_formula, data = df), direction = "backward", trace = 0, k = 2)
}

# ── 4c.  BIC (backward step) ────────────────────────────────────────────────
fit_bic <- function(df) {
  step(lm(full_formula, data = df), direction = "backward", trace = 0,
       k = log(nrow(df)))
}

# ── 4d–4e.  LASSO (lambda.min and lambda.1se via 10-fold CV) ────────────────
fit_lasso <- function(df, lambda_choice = "lambda.min") {
  Xm  <- as.matrix(df[, paste0("X", seq_len(P))])
  Yv  <- df$Y
  cv  <- cv.glmnet(Xm, Yv, alpha = 1, nfolds = 10, standardize = TRUE)
  lam <- cv[[lambda_choice]]
  fit <- glmnet(Xm, Yv, alpha = 1, lambda = lam, standardize = TRUE)
  list(cv = cv, fit = fit, lambda = lam, X = Xm, Y = Yv)
}

# ── 4f–4g.  Elastic Net (alpha=0.5; lambda.min and lambda.1se) ──────────────
fit_enet <- function(df, alpha_en = ALPHA_ENET, lambda_choice = "lambda.min") {
  Xm  <- as.matrix(df[, paste0("X", seq_len(P))])
  Yv  <- df$Y
  cv  <- cv.glmnet(Xm, Yv, alpha = alpha_en, nfolds = 10, standardize = TRUE)
  lam <- cv[[lambda_choice]]
  fit <- glmnet(Xm, Yv, alpha = alpha_en, lambda = lam, standardize = TRUE)
  list(cv = cv, fit = fit, lambda = lam, X = Xm, Y = Yv)
}


# ─────────────────────────────────────────────────────────────────────────────
# 5.  EXTRACTING SELECTED VARIABLES AND COEFFICIENTS
# ─────────────────────────────────────────────────────────────────────────────

# lm-based methods → named coefficient vector (excluding intercept)
extract_lm <- function(mod) {
  cf   <- coef(mod)
  vars <- setdiff(names(cf), "(Intercept)")
  list(vars = vars, coefs = cf[vars])
}

# glmnet-based methods → selected vars + shrinkage coefs
extract_glmnet <- function(obj) {
  cf_mat   <- as.matrix(coef(obj$fit))           # (p+1) × 1
  all_vars <- rownames(cf_mat)[-1]               # drop intercept row
  cf_vals  <- cf_mat[-1, 1]
  sel      <- all_vars[cf_vals != 0]
  list(vars = sel, coefs = setNames(cf_vals[cf_vals != 0], sel))
}


# ─────────────────────────────────────────────────────────────────────────────
# 6.  METRICS FOR ONE REPLICATION
# ─────────────────────────────────────────────────────────────────────────────

#' @param sel_vars  character vector of selected variable names
#' @param sel_coefs named numeric vector of coefficients for selected vars
#' @param df        data.frame used to fit the model (for CI refit)
#' @param beta_true full numeric vector of true coefficients (length P)
#'
#' Returns a list: TPR, FPR, bias (length P_TRUE), covered (length P_TRUE)

evaluate_replication <- function(sel_vars, sel_coefs, df, beta_true) {

  # --- TP / FP rates ---
  TP  <- sum(TRUE_VARS  %in% sel_vars)
  FP  <- sum(NULL_VARS  %in% sel_vars)
  TPR <- TP / P_TRUE
  FPR <- FP / P_NULL

  # --- Per-variable coefficient estimate (0 if not selected) ---
  coef_est <- setNames(rep(0, P_TRUE), TRUE_VARS)
  for (v in intersect(TRUE_VARS, sel_vars)) coef_est[v] <- sel_coefs[v]

  bias <- coef_est - beta_true[seq_len(P_TRUE)]

  # --- Naive post-selection 95% CI coverage (refit with OLS on selected set) ---
  covered <- setNames(rep(NA_real_, P_TRUE), TRUE_VARS)

  if (length(sel_vars) > 0) {
    refit_formula <- as.formula(paste("Y ~", paste(sel_vars, collapse = " + ")))
    tryCatch({
      refit <- lm(refit_formula, data = df)
      ci    <- suppressMessages(confint(refit, level = 0.95))
      for (v in TRUE_VARS) {
        bt <- beta_true[which(TRUE_VARS == v)]
        if (v %in% rownames(ci)) {
          covered[v] <- as.numeric(ci[v, 1] <= bt && bt <= ci[v, 2])
        }
        # NA if variable was not in the final model (no CI available)
      }
    }, error = function(e) NULL)
  }

  list(TPR = TPR, FPR = FPR,
       bias = bias, coef_est = coef_est, covered = covered)
}


# ─────────────────────────────────────────────────────────────────────────────
# 7.  ONE-CONDITION SIMULATION RUNNER
# ─────────────────────────────────────────────────────────────────────────────

run_condition <- function(n, rho, B = N_SIM) {

  # Storage matrices: rows = replications, cols = methods
  tpr_mat  <- matrix(NA_real_, B, length(METHOD_NAMES),
                     dimnames = list(NULL, METHOD_NAMES))
  fpr_mat  <- matrix(NA_real_, B, length(METHOD_NAMES),
                     dimnames = list(NULL, METHOD_NAMES))

  # Bias and coverage: B × P_TRUE × n_methods
  bias_arr <- array(NA_real_,
                    dim      = c(B, P_TRUE, length(METHOD_NAMES)),
                    dimnames = list(NULL, TRUE_VARS, METHOD_NAMES))
  cov_arr  <- array(NA_real_,
                    dim      = c(B, P_TRUE, length(METHOD_NAMES)),
                    dimnames = list(NULL, TRUE_VARS, METHOD_NAMES))

  for (b in seq_len(B)) {

    df <- generate_data(n, P, rho, BETA_TRUE, SIGMA_E)

    # A small helper to run one method safely and store results
    store <- function(method_name, sel_vars, sel_coefs) {
      ev <- evaluate_replication(sel_vars, sel_coefs, df, BETA_TRUE)
      tpr_mat[b, method_name]    <<- ev$TPR
      fpr_mat[b, method_name]    <<- ev$FPR
      bias_arr[b, , method_name] <<- ev$bias
      cov_arr[b,  , method_name] <<- ev$covered
    }

    # ── Backward p-value ──
    tryCatch({
      m1 <- fit_backward_pval(df)
      r1 <- extract_lm(m1)
      store("Backward_pval", r1$vars, r1$coefs)
    }, error = function(e) NULL)

    # ── AIC ──
    tryCatch({
      m2 <- fit_aic(df)
      r2 <- extract_lm(m2)
      store("AIC", r2$vars, r2$coefs)
    }, error = function(e) NULL)

    # ── BIC ──
    tryCatch({
      m3 <- fit_bic(df)
      r3 <- extract_lm(m3)
      store("BIC", r3$vars, r3$coefs)
    }, error = function(e) NULL)

    # ── LASSO λ.min ──
    tryCatch({
      m4 <- fit_lasso(df, "lambda.min")
      r4 <- extract_glmnet(m4)
      store("LASSO_min", r4$vars, r4$coefs)
    }, error = function(e) NULL)

    # ── LASSO λ.1se ──
    tryCatch({
      m5 <- fit_lasso(df, "lambda.1se")
      r5 <- extract_glmnet(m5)
      store("LASSO_1se", r5$vars, r5$coefs)
    }, error = function(e) NULL)

    # ── Elastic Net λ.min ──
    tryCatch({
      m6 <- fit_enet(df, ALPHA_ENET, "lambda.min")
      r6 <- extract_glmnet(m6)
      store("ENet_min", r6$vars, r6$coefs)
    }, error = function(e) NULL)

    # ── Elastic Net λ.1se ──
    tryCatch({
      m7 <- fit_enet(df, ALPHA_ENET, "lambda.1se")
      r7 <- extract_glmnet(m7)
      store("ENet_1se", r7$vars, r7$coefs)
    }, error = function(e) NULL)

  } # end replication loop

  list(n       = n,
       rho     = rho,
       tpr     = tpr_mat,
       fpr     = fpr_mat,
       bias    = bias_arr,
       covered = cov_arr)
}


# ─────────────────────────────────────────────────────────────────────────────
# 8.  RUN ALL CONDITIONS
# ─────────────────────────────────────────────────────────────────────────────

conditions   <- expand.grid(n = N_VALS, rho = RHO_VALS, stringsAsFactors = FALSE)
n_conditions <- nrow(conditions)
all_results  <- vector("list", n_conditions)

cat("=======================================================\n")
cat("Project 4 Simulation — BIOS 6624\n")
cat(sprintf("B = %d replications × %d conditions\n", N_SIM, n_conditions))
cat("=======================================================\n\n")

for (ci in seq_len(n_conditions)) {
  n_ci  <- conditions$n[ci]
  rho_ci <- conditions$rho[ci]
  cat(sprintf("[%d/%d]  N = %d,  rho = %.2f  ... ",
              ci, n_conditions, n_ci, rho_ci))
  t0 <- proc.time()
  all_results[[ci]] <- run_condition(n_ci, rho_ci, B = N_SIM)
  elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
  cat(sprintf("done (%.1f s)\n", elapsed))
}

cat("\nSimulation complete.\n\n")


# ─────────────────────────────────────────────────────────────────────────────
# 9.  BUILD SUMMARY DATA FRAME
# ─────────────────────────────────────────────────────────────────────────────

# Aggregate across replications: mean TPR, FPR, bias, coverage
summary_rows <- lapply(all_results, function(res) {
  lapply(METHOD_NAMES, function(m) {
    # Mean bias across X1-X5 and replications
    mean_bias_per_sim <- apply(res$bias[, , m], 1, mean, na.rm = TRUE)
    # Mean coverage across X1-X5 and replications
    mean_cov_per_sim  <- apply(res$covered[, , m], 1, mean, na.rm = TRUE)

    data.frame(
      n         = res$n,
      rho       = res$rho,
      method    = m,
      mean_tpr  = mean(res$tpr[, m],  na.rm = TRUE),
      mean_fpr  = mean(res$fpr[, m],  na.rm = TRUE),
      mean_bias = mean(mean_bias_per_sim, na.rm = TRUE),
      mean_cov  = mean(mean_cov_per_sim,  na.rm = TRUE),
      sd_tpr    = sd(res$tpr[, m],  na.rm = TRUE),
      sd_fpr    = sd(res$fpr[, m],  na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
})

summary_df <- do.call(rbind, do.call(c, summary_rows))

# Add readable factor labels
summary_df$Method_label <- factor(summary_df$method,
  levels = METHOD_NAMES, labels = METHOD_LABELS)
summary_df$n_label   <- paste0("N = ",   summary_df$n)
summary_df$rho_label <- paste0("\u03c1 = ", summary_df$rho)   # ρ symbol

# Ordered factor for ρ labels (for correct facet ordering)
summary_df$rho_label <- factor(summary_df$rho_label,
  levels = paste0("\u03c1 = ", RHO_VALS))

# Ordered factor for N labels
summary_df$n_label <- factor(summary_df$n_label,
  levels = paste0("N = ", N_VALS))


# ─────────────────────────────────────────────────────────────────────────────
# 10. SAVE HELPERS
# ─────────────────────────────────────────────────────────────────────────────

# Save a data.frame as both .tex and .png
save_table <- function(df, filename, caption = "", digits = 3) {

  # ── .tex ──────────────────────────────────────────────────────────────────
  tex_path <- paste0(tab_dir, filename, ".tex")
  xt <- xtable::xtable(df, caption = caption, digits = digits)
  print(xt,
        file                  = tex_path,
        include.rownames      = FALSE,
        booktabs              = TRUE,
        caption.placement     = "top",
        sanitize.text.function = identity,
        comment               = FALSE)
  cat("  Saved:", tex_path, "\n")

  # ── .png ──────────────────────────────────────────────────────────────────
  png_path <- paste0(tab_dir, filename, ".png")

  tt <- gridExtra::ttheme_minimal(
    core    = list(fg_params = list(fontsize = 9),
                   padding   = grid::unit(c(4, 4), "mm")),
    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"),
                   padding   = grid::unit(c(4, 4), "mm"))
  )
  tbl_grob <- gridExtra::tableGrob(df, rows = NULL, theme = tt)

  w_tbl <- grid::convertWidth(sum(tbl_grob$widths),   "in", valueOnly = TRUE)
  h_tbl <- grid::convertHeight(sum(tbl_grob$heights), "in", valueOnly = TRUE)
  h_cap  <- if (nchar(caption) > 0) 0.5 else 0
  w_tot  <- max(w_tbl + 0.8, 8)
  h_tot  <- max(h_tbl + h_cap + 0.3, 2)

  png(png_path, width = w_tot, height = h_tot,
      units = "in", res = 200, bg = "white")
  grid::grid.newpage()

  if (nchar(caption) > 0) {
    grid::grid.text(caption,
                    x    = grid::unit(0.4, "in"),
                    y    = grid::unit(h_tot - 0.28, "in"),
                    gp   = grid::gpar(fontsize = 11, fontface = "bold"),
                    just = "left", default.units = "in")
  }

  x_off <- (w_tot - w_tbl) / 2
  grid::pushViewport(grid::viewport(
    x      = grid::unit(x_off, "in"),
    y      = grid::unit(0.15, "in"),
    width  = grid::unit(w_tbl, "in"),
    height = grid::unit(h_tbl, "in"),
    just   = c("left", "bottom")
  ))
  grid::grid.draw(tbl_grob)
  grid::popViewport()
  dev.off()
  cat("  Saved:", png_path, "\n")
}

# Save a ggplot as .png
save_fig <- function(plot_obj, filename, w = 14, h = 9) {
  path <- paste0(fig_dir, filename, ".png")
  ggplot2::ggsave(path, plot_obj, width = w, height = h, dpi = 200, bg = "white")
  cat("  Saved:", path, "\n")
}

# Common ggplot theme
proj_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    axis.text.x      = ggplot2::element_text(angle = 40, hjust = 1, size = 8.5),
    legend.position  = "none",
    strip.background = ggplot2::element_rect(fill = "grey20"),
    strip.text       = ggplot2::element_text(color = "white", face = "bold", size = 10),
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle    = ggplot2::element_text(size = 10, color = "grey30")
  )


# ─────────────────────────────────────────────────────────────────────────────
# 11. TABLE 1 — SIMULATION PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────

param_table <- data.frame(
  Parameter = c(
    "Replications (B)",
    "Total predictors (p)",
    "Signal predictors",
    "Null predictors",
    "True coefficients (X1–X5)",
    "Null coefficients (X6–X20)",
    "Error standard deviation (σ)",
    "Sample sizes (N)",
    "Correlation structure",
    "Correlation values (ρ)",
    "Backward selection threshold",
    "LASSO lambda choices",
    "Elastic net mixing (α)",
    "Elastic net lambda choices",
    "CV folds (LASSO / E-Net)"
  ),
  Value = c(
    format(N_SIM, big.mark = ","),
    as.character(P),
    paste0("X1–X", P_TRUE),
    paste0("X", P_TRUE + 1, "–X", P),
    paste(round(BETA_TRUE[1:P_TRUE], 4), collapse = ", "),
    "0",
    as.character(SIGMA_E),
    paste(N_VALS, collapse = ", "),
    "Exchangeable (compound symmetry)",
    paste(RHO_VALS, collapse = ", "),
    paste0("p-value < ", ALPHA_BACK),
    "lambda.min,  lambda.1se",
    as.character(ALPHA_ENET),
    "lambda.min,  lambda.1se",
    "10"
  ),
  stringsAsFactors = FALSE
)

cat("\nSaving Table 1 — Simulation Parameters\n")
save_table(param_table,
           "Table1_SimulationParameters",
           caption = "Table 1. Simulation Parameter Settings")


# ─────────────────────────────────────────────────────────────────────────────
# 12. TABLE 2 — MEAN TPR AND FPR (ACROSS ALL CONDITIONS)
# ─────────────────────────────────────────────────────────────────────────────

tpr_fpr_out <- summary_df %>%
  dplyr::transmute(
    N      = as.character(n),
    rho    = as.character(rho),
    Method = as.character(Method_label),
    TPR    = sprintf("%.3f", mean_tpr),
    FPR    = sprintf("%.3f", mean_fpr)
  ) %>%
  dplyr::arrange(N, rho, Method)

cat("\nSaving Table 2 — TPR and FPR\n")
save_table(tpr_fpr_out,
           "Table2_TPR_FPR",
           caption = paste0(
             "Table 2. Mean True Positive Rate (TPR) and False Positive Rate (FPR) ",
             "for Each Method and Simulation Condition (", N_SIM, " replications)"))


# ─────────────────────────────────────────────────────────────────────────────
# 13. TABLE 3 — MEAN BIAS AND 95% CI COVERAGE
# ─────────────────────────────────────────────────────────────────────────────

bias_cov_out <- summary_df %>%
  dplyr::transmute(
    N        = as.character(n),
    rho      = as.character(rho),
    Method   = as.character(Method_label),
    Bias     = sprintf("%.4f", mean_bias),
    Coverage = sprintf("%.3f", mean_cov)
  ) %>%
  dplyr::arrange(N, rho, Method)
names(bias_cov_out)[5] <- "95% CI Coverage"

cat("\nSaving Table 3 — Bias and Coverage\n")
save_table(bias_cov_out,
           "Table3_Bias_Coverage",
           caption = paste0(
             "Table 3. Mean Coefficient Bias and Empirical 95% CI Coverage ",
             "for True Signal Variables X1–X5 (post-selection OLS refit)"))


# ─────────────────────────────────────────────────────────────────────────────
# 14. TABLE 4 — CONFUSION MATRIX: N=250, rho=0.35 (representative condition)
# ─────────────────────────────────────────────────────────────────────────────

cm_idx <- which(
  vapply(all_results, function(r) r$n == 250 && r$rho == 0.35, logical(1))
)

if (length(cm_idx) > 0) {
  res_cm <- all_results[[cm_idx[1]]]

  cm_df <- do.call(rbind, lapply(seq_along(METHOD_NAMES), function(mi) {
    m       <- METHOD_NAMES[mi]
    tpr_val <- mean(res_cm$tpr[, m], na.rm = TRUE)
    fpr_val <- mean(res_cm$fpr[, m], na.rm = TRUE)
    data.frame(
      Method          = METHOD_LABELS[mi],
      `TP Rate`       = sprintf("%.3f", tpr_val),
      `FN Rate`       = sprintf("%.3f", 1 - tpr_val),
      `FP Rate`       = sprintf("%.3f", fpr_val),
      `TN Rate`       = sprintf("%.3f", 1 - fpr_val),
      check.names     = FALSE,
      stringsAsFactors = FALSE
    )
  }))

  cat("\nSaving Table 4 — Confusion Matrix (N=250, rho=0.35)\n")
  save_table(cm_df,
             "Table4_ConfusionMatrix_N250_rho035",
             caption = paste0(
               "Table 4. Confusion Matrix Rates: N = 250, \\u03c1 = 0.35. ",
               "TP/FN computed over signal variables X1-X5; ",
               "FP/TN computed over null variables X6-X20."))
}


# ─────────────────────────────────────────────────────────────────────────────
# 15. TABLE 5 — PER-VARIABLE SELECTION RATES  (N=250, rho=0 vs rho=0.70)
# ─────────────────────────────────────────────────────────────────────────────

per_var_rows <- lapply(all_results, function(res) {
  lapply(seq_along(METHOD_NAMES), function(mi) {
    m <- METHOD_NAMES[mi]
    sel_rates <- sapply(TRUE_VARS, function(v) {
      # variable selected in rep b ↔ coverage is NOT NA (we only set it when in model)
      mean(!is.na(res$covered[, v, m]), na.rm = TRUE)
    })
    row <- data.frame(
      N      = res$n,
      rho    = res$rho,
      Method = METHOD_LABELS[mi],
      t(round(sel_rates, 3)),
      stringsAsFactors = FALSE
    )
    names(row)[4:8] <- paste0(TRUE_VARS,
      " (beta=", round(BETA_TRUE[1:P_TRUE], 3), ")")
    row
  })
})
per_var_df <- do.call(rbind, do.call(c, per_var_rows))

# Print for two spotlight conditions
pv_subset <- per_var_df[per_var_df$rho %in% c(0, 0.70) &
                         per_var_df$N == 250, ]
pv_subset$N   <- as.character(pv_subset$N)
pv_subset$rho <- as.character(pv_subset$rho)

cat("\nSaving Table 5 — Per-Variable Selection Rates (N=250)\n")
save_table(pv_subset,
           "Table5_PerVariable_SelectionRates",
           caption = paste0(
             "Table 5. Per-Variable Selection Rate (proportion of replications) ",
             "for N = 250 under independence (\\u03c1=0) and high correlation (\\u03c1=0.70)"))


# ─────────────────────────────────────────────────────────────────────────────
# 16. FIGURE 1 — TRUE POSITIVE RATE
# ─────────────────────────────────────────────────────────────────────────────

p_tpr <- ggplot2::ggplot(summary_df,
    ggplot2::aes(x = Method_label, y = mean_tpr, fill = Method_label)) +
  ggplot2::geom_col(color = "black", linewidth = 0.25) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                      color = "red2", linewidth = 0.7) +
  ggplot2::facet_grid(n_label ~ rho_label) +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  ggplot2::labs(
    title    = "Figure 1. Mean True Positive Rate by Method, Sample Size, and Correlation",
    subtitle = "Proportion of replications in which each of X1\u2013X5 was included in the final model (red dashed = perfect)",
    x        = "Variable Selection Method",
    y        = "Mean True Positive Rate"
  ) +
  proj_theme

save_fig(p_tpr, "Figure1_TruePositiveRate")


# ─────────────────────────────────────────────────────────────────────────────
# 17. FIGURE 2 — FALSE POSITIVE RATE
# ─────────────────────────────────────────────────────────────────────────────

p_fpr <- ggplot2::ggplot(summary_df,
    ggplot2::aes(x = Method_label, y = mean_fpr, fill = Method_label)) +
  ggplot2::geom_col(color = "black", linewidth = 0.25) +
  ggplot2::geom_hline(yintercept = ALPHA_BACK, linetype = "dashed",
                      color = "red2", linewidth = 0.7) +
  ggplot2::annotate("text", x = 0.65, y = ALPHA_BACK + 0.008,
                    label = paste0("\u03b1 = ", ALPHA_BACK),
                    color = "red2", size = 3.2) +
  ggplot2::facet_grid(n_label ~ rho_label) +
  ggplot2::scale_fill_brewer(palette = "Set1") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  ggplot2::labs(
    title    = "Figure 2. Mean False Positive Rate by Method, Sample Size, and Correlation",
    subtitle = "Proportion of replications in which any of X6\u2013X20 was included in the final model",
    x        = "Variable Selection Method",
    y        = "Mean False Positive Rate"
  ) +
  proj_theme

save_fig(p_fpr, "Figure2_FalsePositiveRate")


# ─────────────────────────────────────────────────────────────────────────────
# 18. FIGURE 3 — MEAN COEFFICIENT BIAS
# ─────────────────────────────────────────────────────────────────────────────

p_bias <- ggplot2::ggplot(summary_df,
    ggplot2::aes(x = Method_label, y = mean_bias, fill = Method_label)) +
  ggplot2::geom_col(color = "black", linewidth = 0.25) +
  ggplot2::geom_hline(yintercept = 0, linetype = "solid",
                      color = "black", linewidth = 0.5) +
  ggplot2::facet_grid(n_label ~ rho_label) +
  ggplot2::scale_fill_brewer(palette = "Set3") +
  ggplot2::labs(
    title    = "Figure 3. Mean Coefficient Bias for True Signal Variables (X1\u2013X5)",
    subtitle = "Average of (\u03b2\u0302 \u2212 \u03b2) across replications and signal variables; negative = shrinkage bias",
    x        = "Variable Selection Method",
    y        = "Mean Bias"
  ) +
  proj_theme

save_fig(p_bias, "Figure3_CoefficientBias")


# ─────────────────────────────────────────────────────────────────────────────
# 19. FIGURE 4 — 95% CI COVERAGE
# ─────────────────────────────────────────────────────────────────────────────

p_cov <- ggplot2::ggplot(summary_df,
    ggplot2::aes(x = Method_label, y = mean_cov, fill = Method_label)) +
  ggplot2::geom_col(color = "black", linewidth = 0.25) +
  ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed",
                      color = "navy", linewidth = 0.8) +
  ggplot2::annotate("text", x = 0.65, y = 0.955,
                    label = "Nominal 95%", color = "navy", size = 3.2) +
  ggplot2::facet_grid(n_label ~ rho_label) +
  ggplot2::scale_fill_brewer(palette = "Paired") +
  ggplot2::scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.1)) +
  ggplot2::labs(
    title    = "Figure 4. Empirical 95% CI Coverage for True Signal Variables (X1\u2013X5)",
    subtitle = "Post-selection OLS refit; blue dashed = nominal 95% level",
    x        = "Variable Selection Method",
    y        = "Empirical Coverage Probability"
  ) +
  proj_theme

save_fig(p_cov, "Figure4_CICoverage")


# ─────────────────────────────────────────────────────────────────────────────
# 20. FIGURE 5 — TPR vs FPR SCATTER (ROC-STYLE)
# ─────────────────────────────────────────────────────────────────────────────

p_roc <- ggplot2::ggplot(summary_df,
    ggplot2::aes(x = mean_fpr, y = mean_tpr,
                 color = Method_label, shape = Method_label)) +
  ggplot2::geom_point(size = 3.5, alpha = 0.9) +
  ggplot2::geom_abline(slope = 1, intercept = 0,
                       linetype = "dashed", color = "grey55") +
  ggplot2::facet_grid(n_label ~ rho_label) +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::scale_shape_manual(values = c(15, 16, 17, 18, 8, 11, 13)) +
  ggplot2::labs(
    title    = "Figure 5. TPR vs. FPR Scatter by Method (ROC-Style)",
    subtitle = "Upper-left corner = ideal (high sensitivity, low false discovery); diagonal = random",
    x        = "Mean False Positive Rate",
    y        = "Mean True Positive Rate",
    color    = "Method",
    shape    = "Method"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position  = "bottom",
    legend.text      = ggplot2::element_text(size = 8.5),
    strip.background = ggplot2::element_rect(fill = "grey20"),
    strip.text       = ggplot2::element_text(color = "white", face = "bold"),
    plot.title       = ggplot2::element_text(face = "bold", size = 13)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
    shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
  )

save_fig(p_roc, "Figure5_TPR_vs_FPR_scatter", w = 14, h = 10)


# ─────────────────────────────────────────────────────────────────────────────
# 21. FIGURE 6 — PER-VARIABLE SELECTION RATE HEATMAP (N=250, rho=0.35)
# ─────────────────────────────────────────────────────────────────────────────

pv_heat_idx <- which(
  vapply(all_results, function(r) r$n == 250 && r$rho == 0.35, logical(1))
)

if (length(pv_heat_idx) > 0) {
  res_hm <- all_results[[pv_heat_idx[1]]]

  heat_rows <- lapply(seq_along(METHOD_NAMES), function(mi) {
    m  <- METHOD_NAMES[mi]
    sr <- sapply(TRUE_VARS, function(v) {
      mean(!is.na(res_hm$covered[, v, m]), na.rm = TRUE)
    })
    data.frame(
      Method      = METHOD_LABELS[mi],
      Variable    = TRUE_VARS,
      True_Beta   = BETA_TRUE[1:P_TRUE],
      Sel_Rate    = round(sr, 3),
      stringsAsFactors = FALSE
    )
  })
  heat_df <- do.call(rbind, heat_rows)

  heat_df$Method <- factor(heat_df$Method, levels = METHOD_LABELS)
  heat_df$var_label <- paste0(heat_df$Variable,
    "\n(\u03b2 = ", sprintf("%.3f", heat_df$True_Beta), ")")
  heat_df$var_label <- factor(heat_df$var_label,
    levels = unique(heat_df$var_label[order(heat_df$True_Beta)]))

  p_heatmap <- ggplot2::ggplot(heat_df,
      ggplot2::aes(x = Method, y = var_label, fill = Sel_Rate)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Sel_Rate)),
                       size = 4.2, fontface = "bold",
                       color = ifelse(heat_df$Sel_Rate > 0.6, "white", "black")) +
    ggplot2::scale_fill_gradient2(
      low = "white", mid = "#7fcdbb", high = "#08519c",
      midpoint = 0.5, limits = c(0, 1), name = "Selection\nRate"
    ) +
    ggplot2::labs(
      title    = "Figure 6. Per-Variable Selection Rate Heatmap (N = 250, \u03c1 = 0.35)",
      subtitle = paste0("Fraction of ", N_SIM, " replications in which each signal variable was included"),
      x        = "Variable Selection Method",
      y        = "Signal Variable (True \u03b2)"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 35, hjust = 1, size = 9),
      panel.grid       = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      legend.position  = "right"
    )

  save_fig(p_heatmap, "Figure6_SelectionRate_Heatmap_N250_rho035",
           w = 11, h = 6)
}


# ─────────────────────────────────────────────────────────────────────────────
# 22. FIGURE 7 — BIAS DECOMPOSITION BY INDIVIDUAL VARIABLE
#     (shows how shrinkage bias varies with coefficient magnitude)
# ─────────────────────────────────────────────────────────────────────────────

# Use N=500, rho=0 for a clean signal
bias_var_idx <- which(
  vapply(all_results, function(r) r$n == 500 && r$rho == 0, logical(1))
)

if (length(bias_var_idx) > 0) {
  res_bv <- all_results[[bias_var_idx[1]]]

  bv_rows <- lapply(seq_along(METHOD_NAMES), function(mi) {
    m <- METHOD_NAMES[mi]
    # Mean bias per variable (collapse across reps)
    per_var_bias <- apply(res_bv$bias[, , m], 2, mean, na.rm = TRUE)
    data.frame(
      Method    = METHOD_LABELS[mi],
      Variable  = TRUE_VARS,
      True_Beta = BETA_TRUE[1:P_TRUE],
      Mean_Bias = per_var_bias,
      stringsAsFactors = FALSE
    )
  })
  bv_df <- do.call(rbind, bv_rows)
  bv_df$Method <- factor(bv_df$Method, levels = METHOD_LABELS)
  bv_df$var_label <- paste0(bv_df$Variable,
    " (\u03b2=", sprintf("%.3f", bv_df$True_Beta), ")")

  p_bias_var <- ggplot2::ggplot(bv_df,
      ggplot2::aes(x = var_label, y = Mean_Bias,
                   group = Method, color = Method, shape = Method)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.75) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::scale_shape_manual(values = c(15, 16, 17, 18, 8, 11, 13)) +
    ggplot2::labs(
      title    = "Figure 7. Per-Variable Coefficient Bias (N = 500, \u03c1 = 0)",
      subtitle = "Shrinkage methods show downward bias; stepwise methods show upward bias for weaker signals",
      x        = "Signal Variable (True \u03b2)",
      y        = "Mean Bias (\u03b2\u0302 \u2212 \u03b2)",
      color    = "Method", shape = "Method"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      plot.title      = ggplot2::element_text(face = "bold", size = 13)
    )

  save_fig(p_bias_var, "Figure7_PerVariable_Bias_N500_rho0", w = 11, h = 6)
}


# ─────────────────────────────────────────────────────────────────────────────
# 23. FIGURE 8 — EFFECT OF CORRELATION ON TPR (line plot)
# ─────────────────────────────────────────────────────────────────────────────

p_rho_tpr <- ggplot2::ggplot(summary_df,
    ggplot2::aes(x = factor(rho), y = mean_tpr,
                 group = Method_label, color = Method_label, shape = Method_label)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 3) +
  ggplot2::facet_wrap(~ n_label) +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::scale_shape_manual(values = c(15, 16, 17, 18, 8, 11, 13)) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ggplot2::labs(
    title    = "Figure 8. Effect of Predictor Correlation on True Positive Rate",
    subtitle = "Each line = one method; x-axis = exchangeable correlation \u03c1",
    x        = "Correlation (\u03c1)",
    y        = "Mean True Positive Rate",
    color    = "Method", shape = "Method"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position  = "bottom",
    legend.text      = ggplot2::element_text(size = 8.5),
    strip.background = ggplot2::element_rect(fill = "grey20"),
    strip.text       = ggplot2::element_text(color = "white", face = "bold"),
    plot.title       = ggplot2::element_text(face = "bold", size = 13)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 2),
    shape = ggplot2::guide_legend(nrow = 2)
  )

save_fig(p_rho_tpr, "Figure8_Correlation_Effect_on_TPR", w = 12, h = 7)


# ─────────────────────────────────────────────────────────────────────────────
# 24. SAVE FULL RESULTS OBJECT
# ─────────────────────────────────────────────────────────────────────────────

saveRDS(
  list(
    simulation_params = list(
      N_SIM    = N_SIM,   P        = P,     P_TRUE   = P_TRUE,
      N_VALS   = N_VALS,  RHO_VALS = RHO_VALS,
      BETA     = BETA_TRUE, SIGMA_E  = SIGMA_E,
      ALPHA_BACK = ALPHA_BACK, ALPHA_ENET = ALPHA_ENET
    ),
    all_results = all_results,
    summary_df  = summary_df
  ),
  file = "../Results/simulation_results.rds"
)

cat("\n=======================================================\n")
cat("All outputs saved.\n")
cat("  Tables  →  ../Results/Tables/   (*.tex and *.png)\n")
cat("  Figures →  ../Results/Figures/  (*.png)\n")
cat("  RDS     →  ../Results/simulation_results.rds\n")
cat("=======================================================\n")

# Print a console summary
cat("\n--- Simulation Summary ---\n")
print(
  summary_df[, c("n", "rho", "method", "mean_tpr", "mean_fpr",
                 "mean_bias", "mean_cov")],
  digits = 3, row.names = FALSE
)
