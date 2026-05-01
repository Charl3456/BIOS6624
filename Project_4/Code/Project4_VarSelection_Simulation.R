# Project 4: Variable Selection Algorithms — BIOS 6624
# Place this file in: Project_4/Code/
# Outputs go to:      Project_4/Results/Figures/ and Project_4/Results/Tables/
#
# Run via: RStudio > Session > Set WD > To Source File Location > Source
#      or: cd Project_4/Code && Rscript Project4_VarSelection_Simulation.R


# ── Path setup ────────────────────────────────────────────────────────────────

get_script_dir <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0)
    return(normalizePath(dirname(sub("^--file=", "", file_arg)), mustWork = FALSE))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    if (!is.null(ctx) && nchar(ctx$path) > 0)
      return(normalizePath(dirname(ctx$path), mustWork = FALSE))
  }
  for (fr in rev(sys.frames())) {
    ofile <- tryCatch(get("ofile", envir = fr), error = function(e) NULL)
    if (!is.null(ofile) && nchar(ofile) > 0)
      return(normalizePath(dirname(ofile), mustWork = FALSE))
  }
  warning("Could not detect script location; using working directory: ", getwd())
  getwd()
}

SCRIPT_DIR  <- get_script_dir()
PROJECT_DIR <- dirname(SCRIPT_DIR)
RESULTS_DIR <- file.path(PROJECT_DIR, "Results")
FIG_DIR     <- file.path(RESULTS_DIR, "Figures")
TAB_DIR     <- file.path(RESULTS_DIR, "Tables")

for (d in c(RESULTS_DIR, FIG_DIR, TAB_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d))
    stop("Could not create directory: ", d)
}


# ── Packages ──────────────────────────────────────────────────────────────────

required_pkgs <- c("MASS", "glmnet", "xtable", "ggplot2",
                   "gridExtra", "grid", "dplyr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

if (!requireNamespace("hdrm", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools", repos = "https://cloud.r-project.org")
  tryCatch(devtools::install_github("pbreheny/hdrm"), error = function(e) NULL)
}


# ── Parameters ────────────────────────────────────────────────────────────────

set.seed(2026)

N_SIM      <- 1000
P          <- 20
P_TRUE     <- 5
P_NULL     <- P - P_TRUE
N_VALS     <- c(250, 500)
RHO_VALS   <- c(0, 0.35, 0.70)   # exchangeable correlation (class decision)
BETA_TRUE  <- c(0.5/3, 1/3, 1.5/3, 2.0/3, 2.5/3, rep(0, P_NULL))
SIGMA_E    <- 1.0
ALPHA_BACK <- 0.10
ALPHA_ENET <- 0.50

TRUE_VARS <- paste0("X", seq_len(P_TRUE))
NULL_VARS <- paste0("X", (P_TRUE + 1):P)

METHOD_NAMES <- c("Backward_pval", "AIC", "BIC",
                  "LASSO_min", "LASSO_1se", "ENet_min", "ENet_1se")

METHOD_LABELS <- c("Backward (p-val)", "AIC", "BIC",
                   "LASSO (lambda.min)", "LASSO (lambda.1se)",
                   "Elastic Net (lambda.min)", "Elastic Net (lambda.1se)")


# ── Data generation ───────────────────────────────────────────────────────────

make_exch_cov <- function(p, rho) {
  S <- matrix(rho, p, p); diag(S) <- 1; S
}

generate_data <- function(n, p, rho, beta, sigma_e = 1) {
  X           <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = make_exch_cov(p, rho))
  colnames(X) <- paste0("X", seq_len(p))
  data.frame(Y = drop(X %*% beta) + rnorm(n, 0, sigma_e), X)
}


# ── Model fitting ─────────────────────────────────────────────────────────────

full_formula <- as.formula(
  paste("Y ~", paste(paste0("X", seq_len(P)), collapse = " + "))
)

fit_backward_pval <- function(df, alpha = ALPHA_BACK) {
  mod <- lm(full_formula, data = df)
  repeat {
    pv <- summary(mod)$coefficients[-1, "Pr(>|t|)", drop = FALSE]
    if (nrow(pv) == 0 || max(pv) <= alpha) break
    current <- setdiff(attr(terms(mod), "term.labels"), rownames(pv)[which.max(pv)])
    if (length(current) == 0) { mod <- lm(Y ~ 1, data = df); break }
    mod <- lm(as.formula(paste("Y ~", paste(current, collapse = " + "))), data = df)
  }
  mod
}

fit_aic <- function(df)
  step(lm(full_formula, data = df), direction = "backward", trace = 0, k = 2)

fit_bic <- function(df)
  step(lm(full_formula, data = df), direction = "backward", trace = 0, k = log(nrow(df)))

fit_lasso <- function(df, lambda_choice = "lambda.min") {
  Xm <- as.matrix(df[, paste0("X", seq_len(P))]); Yv <- df$Y
  cv  <- cv.glmnet(Xm, Yv, alpha = 1, nfolds = 10, standardize = TRUE)
  list(fit = glmnet(Xm, Yv, alpha = 1, lambda = cv[[lambda_choice]], standardize = TRUE))
}

fit_enet <- function(df, alpha_en = ALPHA_ENET, lambda_choice = "lambda.min") {
  Xm <- as.matrix(df[, paste0("X", seq_len(P))]); Yv <- df$Y
  cv  <- cv.glmnet(Xm, Yv, alpha = alpha_en, nfolds = 10, standardize = TRUE)
  list(fit = glmnet(Xm, Yv, alpha = alpha_en, lambda = cv[[lambda_choice]], standardize = TRUE))
}


# ── Extraction ────────────────────────────────────────────────────────────────

extract_lm <- function(mod) {
  cf <- coef(mod); vars <- setdiff(names(cf), "(Intercept)")
  list(vars = vars, coefs = cf[vars])
}

extract_glmnet <- function(obj) {
  cm  <- as.matrix(coef(obj$fit))[-1, 1]
  sel <- names(cm)[cm != 0]
  list(vars = sel, coefs = cm[sel])
}


# ── Per-replication metrics ───────────────────────────────────────────────────

evaluate_replication <- function(sel_vars, sel_coefs, df, beta_true) {
  TPR <- sum(TRUE_VARS %in% sel_vars) / P_TRUE
  FPR <- sum(NULL_VARS %in% sel_vars) / P_NULL
  
  coef_est <- setNames(rep(0, P_TRUE), TRUE_VARS)
  for (v in intersect(TRUE_VARS, sel_vars)) coef_est[v] <- sel_coefs[v]
  bias <- coef_est - beta_true[seq_len(P_TRUE)]
  
  covered <- setNames(rep(NA_real_, P_TRUE), TRUE_VARS)
  if (length(sel_vars) > 0) {
    tryCatch({
      ci <- suppressMessages(confint(
        lm(as.formula(paste("Y ~", paste(sel_vars, collapse = " + "))), data = df),
        level = 0.95))
      for (v in TRUE_VARS)
        if (v %in% rownames(ci))
          covered[v] <- as.numeric(ci[v, 1] <= beta_true[which(TRUE_VARS == v)] &&
                                     beta_true[which(TRUE_VARS == v)] <= ci[v, 2])
    }, error = function(e) NULL)
  }
  
  list(TPR = TPR, FPR = FPR, bias = bias, covered = covered)
}


# ── Simulation loop ───────────────────────────────────────────────────────────

run_condition <- function(n, rho, B = N_SIM) {
  tpr_mat  <- matrix(NA_real_, B, length(METHOD_NAMES), dimnames = list(NULL, METHOD_NAMES))
  fpr_mat  <- matrix(NA_real_, B, length(METHOD_NAMES), dimnames = list(NULL, METHOD_NAMES))
  bias_arr <- array(NA_real_, c(B, P_TRUE, length(METHOD_NAMES)),
                    dimnames = list(NULL, TRUE_VARS, METHOD_NAMES))
  cov_arr  <- array(NA_real_, c(B, P_TRUE, length(METHOD_NAMES)),
                    dimnames = list(NULL, TRUE_VARS, METHOD_NAMES))
  
  for (b in seq_len(B)) {
    df <- generate_data(n, P, rho, BETA_TRUE, SIGMA_E)
    
    store <- function(mname, vars, coefs) {
      ev                   <- evaluate_replication(vars, coefs, df, BETA_TRUE)
      tpr_mat[b, mname]    <<- ev$TPR
      fpr_mat[b, mname]    <<- ev$FPR
      bias_arr[b, , mname] <<- ev$bias
      cov_arr[b,  , mname] <<- ev$covered
    }
    
    tryCatch({ r <- extract_lm(fit_backward_pval(df));           store("Backward_pval", r$vars, r$coefs) }, error = function(e) NULL)
    tryCatch({ r <- extract_lm(fit_aic(df));                     store("AIC",           r$vars, r$coefs) }, error = function(e) NULL)
    tryCatch({ r <- extract_lm(fit_bic(df));                     store("BIC",           r$vars, r$coefs) }, error = function(e) NULL)
    tryCatch({ r <- extract_glmnet(fit_lasso(df, "lambda.min")); store("LASSO_min",     r$vars, r$coefs) }, error = function(e) NULL)
    tryCatch({ r <- extract_glmnet(fit_lasso(df, "lambda.1se")); store("LASSO_1se",     r$vars, r$coefs) }, error = function(e) NULL)
    tryCatch({ r <- extract_glmnet(fit_enet(df, ALPHA_ENET, "lambda.min")); store("ENet_min", r$vars, r$coefs) }, error = function(e) NULL)
    tryCatch({ r <- extract_glmnet(fit_enet(df, ALPHA_ENET, "lambda.1se")); store("ENet_1se", r$vars, r$coefs) }, error = function(e) NULL)
  }
  
  list(n = n, rho = rho, tpr = tpr_mat, fpr = fpr_mat, bias = bias_arr, covered = cov_arr)
}

conditions  <- expand.grid(n = N_VALS, rho = RHO_VALS, stringsAsFactors = FALSE)
all_results <- vector("list", nrow(conditions))

for (ci in seq_len(nrow(conditions)))
  all_results[[ci]] <- run_condition(conditions$n[ci], conditions$rho[ci], B = N_SIM)


# ── Summary data frame ────────────────────────────────────────────────────────

summary_df <- do.call(rbind, do.call(c, lapply(all_results, function(res) {
  lapply(METHOD_NAMES, function(m) {
    mb <- apply(res$bias[, , m],    1, mean, na.rm = TRUE)
    mc <- apply(res$covered[, , m], 1, mean, na.rm = TRUE)
    data.frame(n = res$n, rho = res$rho, method = m,
               mean_tpr  = mean(res$tpr[, m], na.rm = TRUE),
               mean_fpr  = mean(res$fpr[, m], na.rm = TRUE),
               mean_bias = mean(mb, na.rm = TRUE),
               mean_cov  = mean(mc, na.rm = TRUE),
               sd_tpr    = sd(res$tpr[, m], na.rm = TRUE),
               sd_fpr    = sd(res$fpr[, m], na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
})))

summary_df$Method_label <- factor(summary_df$method, levels = METHOD_NAMES, labels = METHOD_LABELS)
summary_df$n_label      <- factor(paste0("N = ",    summary_df$n),   levels = paste0("N = ",    N_VALS))
summary_df$rho_label    <- factor(paste0("rho = ",  summary_df$rho), levels = paste0("rho = ",  RHO_VALS))


# ── Output helpers ────────────────────────────────────────────────────────────

save_table <- function(df, filename, caption = "") {
  tex_path <- file.path(TAB_DIR, paste0(filename, ".tex"))
  print(xtable::xtable(df, caption = caption), file = tex_path,
        include.rownames = FALSE, booktabs = TRUE,
        caption.placement = "top", sanitize.text.function = identity, comment = FALSE)
  
  png_path <- file.path(TAB_DIR, paste0(filename, ".png"))
  tt  <- gridExtra::ttheme_minimal(
    core    = list(fg_params = list(fontsize = 9), padding = grid::unit(c(4, 4), "mm")),
    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"), padding = grid::unit(c(4, 4), "mm")))
  tg    <- gridExtra::tableGrob(df, rows = NULL, theme = tt)
  w_tbl <- grid::convertWidth( sum(tg$widths),  "in", valueOnly = TRUE)
  h_tbl <- grid::convertHeight(sum(tg$heights), "in", valueOnly = TRUE)
  h_cap <- if (nchar(caption) > 0) 0.45 else 0.1
  w_tot <- max(w_tbl + 0.8, 7)
  h_tot <- max(h_tbl + h_cap + 0.2, 2)
  
  grDevices::png(png_path, width = w_tot, height = h_tot, units = "in", res = 200, bg = "white")
  grid::grid.newpage()
  if (nchar(caption) > 0)
    grid::grid.text(caption, x = grid::unit(0.35, "in"), y = grid::unit(h_tot - 0.28, "in"),
                    gp = grid::gpar(fontsize = 11, fontface = "bold"),
                    just = "left", default.units = "in")
  grid::pushViewport(grid::viewport(
    x = grid::unit((w_tot - w_tbl) / 2, "in"), y = grid::unit(0.12, "in"),
    width = grid::unit(w_tbl, "in"), height = grid::unit(h_tbl, "in"), just = c("left", "bottom")))
  grid::grid.draw(tg)
  grid::popViewport()
  grDevices::dev.off()
}

save_fig <- function(p, filename, w = 14, h = 9)
  ggplot2::ggsave(file.path(FIG_DIR, paste0(filename, ".png")),
                  p, width = w, height = h, dpi = 200, bg = "white")

proj_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    axis.text.x      = ggplot2::element_text(angle = 38, hjust = 1, size = 8.5),
    legend.position  = "none",
    strip.background = ggplot2::element_rect(fill = "grey20"),
    strip.text       = ggplot2::element_text(color = "white", face = "bold", size = 10),
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle    = ggplot2::element_text(size = 9.5, color = "grey35"))


# ── Tables ────────────────────────────────────────────────────────────────────

save_table(
  data.frame(
    Parameter = c("Replications (B)", "Total predictors (p)", "Signal predictors",
                  "Null predictors", "True coefficients (X1-X5)", "Null coefficients (X6-X20)",
                  "Error SD (sigma)", "Sample sizes (N)", "Correlation structure",
                  "Correlation values (rho)", "Backward selection threshold",
                  "LASSO lambda choices", "Elastic net mixing (alpha)",
                  "Elastic net lambda choices", "Cross-validation folds"),
    Value = c(format(N_SIM, big.mark = ","), as.character(P),
              paste0("X1-X", P_TRUE), paste0("X", P_TRUE + 1, "-X", P),
              paste(round(BETA_TRUE[1:P_TRUE], 4), collapse = ", "), "0",
              as.character(SIGMA_E), paste(N_VALS, collapse = ", "),
              "Exchangeable (compound symmetry)", paste(RHO_VALS, collapse = ", "),
              paste0("p-value < ", ALPHA_BACK), "lambda.min, lambda.1se",
              as.character(ALPHA_ENET), "lambda.min, lambda.1se", "10"),
    stringsAsFactors = FALSE),
  "Table1_SimulationParameters",
  "Table 1. Simulation Parameter Settings")

save_table(
  summary_df %>%
    dplyr::transmute(N = as.character(n), rho = as.character(rho),
                     Method = as.character(Method_label),
                     `Mean TPR (SD)` = sprintf("%.3f (%.3f)", mean_tpr, sd_tpr),
                     `Mean FPR (SD)` = sprintf("%.3f (%.3f)", mean_fpr, sd_fpr)) %>%
    dplyr::arrange(N, rho, Method),
  "Table2_TPR_FPR",
  paste0("Table 2. Mean TPR and FPR by Method and Condition (SD in parentheses; B = ", N_SIM, ")"))

save_table(
  summary_df %>%
    dplyr::transmute(N = as.character(n), rho = as.character(rho),
                     Method = as.character(Method_label),
                     `Mean Bias`  = sprintf("%.4f", mean_bias),
                     `95% CI Cov` = sprintf("%.3f", mean_cov)) %>%
    dplyr::arrange(N, rho, Method),
  "Table3_Bias_Coverage",
  "Table 3. Mean Coefficient Bias and 95% CI Coverage for X1-X5 (post-selection OLS refit)")

cm_idx <- which(vapply(all_results, function(r) r$n == 250 && r$rho == 0.35, logical(1)))
if (length(cm_idx) > 0) {
  res_cm <- all_results[[cm_idx[1]]]
  save_table(
    do.call(rbind, lapply(seq_along(METHOD_NAMES), function(mi) {
      m <- METHOD_NAMES[mi]; tpr <- mean(res_cm$tpr[, m], na.rm = TRUE); fpr <- mean(res_cm$fpr[, m], na.rm = TRUE)
      data.frame(Method = METHOD_LABELS[mi], `TP Rate` = sprintf("%.3f", tpr),
                 `FN Rate` = sprintf("%.3f", 1 - tpr), `FP Rate` = sprintf("%.3f", fpr),
                 `TN Rate` = sprintf("%.3f", 1 - fpr), check.names = FALSE, stringsAsFactors = FALSE)
    })),
    "Table4_ConfusionMatrix_N250_rho035",
    "Table 4. Confusion Matrix Rates: N=250, rho=0.35 (TP/FN over X1-X5; FP/TN over X6-X20)")
}

pv_df <- do.call(rbind, do.call(c, lapply(all_results, function(res) {
  lapply(seq_along(METHOD_NAMES), function(mi) {
    m   <- METHOD_NAMES[mi]
    sel <- sapply(TRUE_VARS, function(v) mean(!is.na(res$covered[, v, m]), na.rm = TRUE))
    row <- data.frame(N = res$n, rho = res$rho, Method = METHOD_LABELS[mi],
                      t(round(sel, 3)), stringsAsFactors = FALSE)
    colnames(row)[4:8] <- paste0(TRUE_VARS, " (b=", round(BETA_TRUE[1:P_TRUE], 3), ")")
    row
  })
})))
pv_sub <- pv_df[pv_df$rho %in% c(0, 0.70) & pv_df$N == 250, ]
pv_sub$N <- as.character(pv_sub$N); pv_sub$rho <- as.character(pv_sub$rho)
save_table(pv_sub, "Table5_PerVariable_SelectionRates",
           "Table 5. Per-Variable Selection Rate for N=250, rho=0 and rho=0.70")


# ── Figures ───────────────────────────────────────────────────────────────────

# Okabe-Ito colorblind-safe palette (7 methods)
CB_COLORS <- c(
  "Backward (p-val)"         = "#0072B2",   # blue
  "AIC"                      = "#E69F00",   # orange
  "BIC"                      = "#009E73",   # bluish green
  "LASSO (lambda.min)"       = "#D55E00",   # vermillion
  "LASSO (lambda.1se)"       = "#56B4E9",   # sky blue
  "Elastic Net (lambda.min)" = "#CC79A7",   # reddish purple
  "Elastic Net (lambda.1se)" = "#999999"    # grey
)

# All solid shapes — no fill aesthetic needed, distinguished by color + shape
CB_SHAPES <- c(
  "Backward (p-val)"         = 16,   # circle
  "AIC"                      = 17,   # triangle up
  "BIC"                      = 15,   # square
  "LASSO (lambda.min)"       = 18,   # diamond
  "LASSO (lambda.1se)"       = 7,    # square-X  (visible, not weird)
  "Elastic Net (lambda.min)" = 10,   # circle-plus
  "Elastic Net (lambda.1se)" = 13    # circle-X
)

# ── Figure 1: True Positive Rate ──────────────────────────────────────────────
# coord_cartesian zooms WITHOUT clipping bars (bars still start from 0)
save_fig(
  ggplot2::ggplot(summary_df, ggplot2::aes(x = Method_label, y = mean_tpr, fill = Method_label)) +
    ggplot2::geom_col(color = "grey20", linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
    ggplot2::facet_grid(n_label ~ rho_label) +
    ggplot2::scale_fill_manual(values = CB_COLORS) +
    ggplot2::scale_y_continuous(breaks = seq(0.70, 1.00, 0.05)) +
    ggplot2::coord_cartesian(ylim = c(0.70, 1.02)) +
    ggplot2::labs(title = "Figure 1. Mean True Positive Rate",
                  subtitle = "Proportion of replications X1-X5 was selected (orange dashed = perfect TPR = 1.0)",
                  x = "Variable Selection Method", y = "Mean True Positive Rate") + proj_theme,
  "Figure1_TruePositiveRate")

# ── Figure 2: False Positive Rate ─────────────────────────────────────────────
save_fig(
  ggplot2::ggplot(summary_df, ggplot2::aes(x = Method_label, y = mean_fpr, fill = Method_label)) +
    ggplot2::geom_col(color = "grey20", linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = ALPHA_BACK, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
    ggplot2::facet_grid(n_label ~ rho_label) +
    ggplot2::scale_fill_manual(values = CB_COLORS) +
    ggplot2::scale_y_continuous(limits = c(0, 0.60), breaks = seq(0, 0.60, 0.10)) +
    ggplot2::labs(title = "Figure 2. Mean False Positive Rate",
                  subtitle = paste0("Proportion of replications any null variable X6-X20 was selected; orange dashed = alpha=", ALPHA_BACK),
                  x = "Variable Selection Method", y = "Mean False Positive Rate") + proj_theme,
  "Figure2_FalsePositiveRate")

# ── Figure 3: Coefficient Bias ────────────────────────────────────────────────
save_fig(
  ggplot2::ggplot(summary_df, ggplot2::aes(x = Method_label, y = mean_bias, fill = Method_label)) +
    ggplot2::geom_col(color = "grey20", linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.6, color = "grey20") +
    ggplot2::facet_grid(n_label ~ rho_label) +
    ggplot2::scale_fill_manual(values = CB_COLORS) +
    ggplot2::scale_y_continuous(limits = c(-0.20, 0.02), breaks = seq(-0.20, 0.00, 0.05)) +
    ggplot2::labs(title = "Figure 3. Mean Coefficient Bias for Signal Variables X1-X5",
                  subtitle = "Average of (beta_hat - beta_true); negative = shrinkage bias",
                  x = "Variable Selection Method", y = "Mean Bias") + proj_theme,
  "Figure3_CoefficientBias")

# ── Figure 4: 95% CI Coverage ─────────────────────────────────────────────────
# coord_cartesian zooms WITHOUT clipping bars (bars still start from 0)
save_fig(
  ggplot2::ggplot(summary_df, ggplot2::aes(x = Method_label, y = mean_cov, fill = Method_label)) +
    ggplot2::geom_col(color = "grey20", linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", color = "#0072B2", linewidth = 0.8) +
    ggplot2::facet_grid(n_label ~ rho_label) +
    ggplot2::scale_fill_manual(values = CB_COLORS) +
    ggplot2::scale_y_continuous(breaks = seq(0.88, 1.00, 0.02)) +
    ggplot2::coord_cartesian(ylim = c(0.88, 1.0)) +
    ggplot2::labs(title = "Figure 4. Empirical 95% CI Coverage for Signal Variables X1-X5",
                  subtitle = "Post-selection OLS refit; blue dashed = nominal 95%",
                  x = "Variable Selection Method", y = "Coverage Probability") + proj_theme,
  "Figure4_CICoverage")

# ── Figure 5: TPR vs FPR scatter (ROC-style) ──────────────────────────────────
save_fig(
  ggplot2::ggplot(summary_df, ggplot2::aes(x = mean_fpr, y = mean_tpr,
                                           color = Method_label, shape = Method_label)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
    ggplot2::facet_grid(n_label ~ rho_label) +
    ggplot2::scale_color_manual(values = CB_COLORS) +
    ggplot2::scale_shape_manual(values = CB_SHAPES) +
    ggplot2::scale_x_continuous(limits = c(0, 0.60), breaks = seq(0, 0.60, 0.20)) +
    ggplot2::scale_y_continuous(limits = c(0.78, 1.02), breaks = seq(0.80, 1.00, 0.05)) +
    ggplot2::labs(title = "Figure 5. TPR vs FPR Scatter by Method (ROC-Style)",
                  subtitle = "Upper-left corner = ideal (high sensitivity, low false discovery); diagonal = chance",
                  x = "Mean False Positive Rate", y = "Mean True Positive Rate",
                  color = "Method", shape = "Method") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8.5),
                   strip.background = ggplot2::element_rect(fill = "grey20"),
                   strip.text = ggplot2::element_text(color = "white", face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold", size = 13)) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2), shape = ggplot2::guide_legend(nrow = 2)),
  "Figure5_TPR_vs_FPR_scatter", w = 14, h = 10)

# ── Figure 6: Selection Rate Heatmap (N=250, rho=0.35) ───────────────────────
hm_idx <- which(vapply(all_results, function(r) r$n == 250 && r$rho == 0.35, logical(1)))
if (length(hm_idx) > 0) {
  res_hm <- all_results[[hm_idx[1]]]
  hm_df  <- do.call(rbind, lapply(seq_along(METHOD_NAMES), function(mi) {
    m  <- METHOD_NAMES[mi]
    sr <- sapply(TRUE_VARS, function(v) mean(!is.na(res_hm$covered[, v, m]), na.rm = TRUE))
    data.frame(Method = METHOD_LABELS[mi], Variable = TRUE_VARS,
               True_Beta = BETA_TRUE[1:P_TRUE], Sel_Rate = round(sr, 3), stringsAsFactors = FALSE)
  }))
  hm_df <- hm_df[order(hm_df$True_Beta), ]
  hm_df$var_label <- factor(paste0(hm_df$Variable, "\n(b=", sprintf("%.3f", hm_df$True_Beta), ")"),
                            levels = unique(paste0(hm_df$Variable, "\n(b=", sprintf("%.3f", hm_df$True_Beta), ")")))
  hm_df$Method <- factor(hm_df$Method, levels = METHOD_LABELS)
  
  save_fig(
    ggplot2::ggplot(hm_df, ggplot2::aes(x = Method, y = var_label, fill = Sel_Rate)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Sel_Rate)), size = 4, fontface = "bold",
                         color = ifelse(hm_df$Sel_Rate > 0.55, "white", "black")) +
      ggplot2::scale_fill_gradient2(low = "white", mid = "#7fcdbb", high = "#08519c",
                                    midpoint = 0.5, limits = c(0, 1), name = "Sel. Rate") +
      ggplot2::labs(title = "Figure 6. Per-Variable Selection Rate Heatmap (N=250, rho=0.35)",
                    subtitle = paste0("Fraction of ", N_SIM, " replications each signal variable was retained"),
                    x = "Method", y = "Signal Variable (True beta)") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 9),
                     panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold", size = 13),
                     legend.position = "right"),
    "Figure6_SelectionRate_Heatmap_N250_rho035", w = 11, h = 6)
}

# ── Figure 7: Per-variable bias line (N=500, rho=0) ──────────────────────────
bv_idx <- which(vapply(all_results, function(r) r$n == 500 && r$rho == 0, logical(1)))
if (length(bv_idx) > 0) {
  res_bv <- all_results[[bv_idx[1]]]
  bv_df  <- do.call(rbind, lapply(seq_along(METHOD_NAMES), function(mi) {
    m <- METHOD_NAMES[mi]
    data.frame(Method = METHOD_LABELS[mi], Variable = TRUE_VARS,
               True_Beta = BETA_TRUE[1:P_TRUE],
               Mean_Bias = apply(res_bv$bias[, , m], 2, mean, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  bv_df$Method    <- factor(bv_df$Method, levels = METHOD_LABELS)
  bv_df$var_label <- paste0(bv_df$Variable, " (b=", sprintf("%.3f", bv_df$True_Beta), ")")
  
  save_fig(
    ggplot2::ggplot(bv_df, ggplot2::aes(x = var_label, y = Mean_Bias,
                                        group = Method, color = Method, shape = Method)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 3.5) +
      ggplot2::scale_color_manual(values = CB_COLORS) +
      ggplot2::scale_shape_manual(values = CB_SHAPES) +
      ggplot2::scale_y_continuous(limits = c(-0.20, 0.02), breaks = seq(-0.20, 0.00, 0.05)) +
      ggplot2::labs(title = "Figure 7. Per-Variable Coefficient Bias (N=500, rho=0)",
                    subtitle = "Shrinkage = downward bias; stepwise near-zero bias",
                    x = "Signal Variable (True beta)", y = "Mean Bias", color = "Method", shape = "Method") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(face = "bold", size = 13)),
    "Figure7_PerVariable_Bias_N500_rho0", w = 11, h = 6)
}

# ── Figure 8: Correlation effect on TPR line plot ─────────────────────────────
save_fig(
  ggplot2::ggplot(summary_df, ggplot2::aes(x = factor(rho), y = mean_tpr,
                                           group = Method_label, color = Method_label, shape = Method_label)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::facet_wrap(~ n_label) +
    ggplot2::scale_color_manual(values = CB_COLORS) +
    ggplot2::scale_shape_manual(values = CB_SHAPES) +
    ggplot2::scale_y_continuous(breaks = seq(0.75, 1.00, 0.05)) +
    ggplot2::coord_cartesian(ylim = c(0.75, 1.02)) +
    ggplot2::labs(title = "Figure 8. Effect of Predictor Correlation on True Positive Rate",
                  subtitle = "Elastic net degrades less than LASSO under high correlation; BIC most sensitive",
                  x = "Correlation (rho)", y = "Mean True Positive Rate", color = "Method", shape = "Method") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8.5),
                   strip.background = ggplot2::element_rect(fill = "grey20"),
                   strip.text = ggplot2::element_text(color = "white", face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold", size = 13)) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2), shape = ggplot2::guide_legend(nrow = 2)),
  "Figure8_Correlation_Effect_on_TPR", w = 12, h = 7)


# ── Save RDS ──────────────────────────────────────────────────────────────────

saveRDS(list(params      = list(N_SIM = N_SIM, P = P, P_TRUE = P_TRUE,
                                N_VALS = N_VALS, RHO_VALS = RHO_VALS,
                                BETA = BETA_TRUE, SIGMA_E = SIGMA_E),
             all_results = all_results,
             summary_df  = summary_df),
        file = file.path(RESULTS_DIR, "simulation_results.rds"))