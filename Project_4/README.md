# Project 4 — Variable Selection Methods in Linear Regression: A Simulation Study

## Overview

This project benchmarks **seven variable selection procedures** for linear regression through a Monte Carlo simulation study. The goal is to characterize each method's performance across competing sample sizes and predictor correlation structures, evaluating the accuracy–parsimony tradeoff and the inferential consequences (bias and confidence interval coverage) of post-selection OLS refitting.

---

## Research Questions

1. Which selection method achieves the highest **True Positive Rate (TPR)** — correctly identifying signal variables?
2. Which methods have the lowest **False Positive Rate (FPR)** — avoiding null variable inclusion?
3. Do penalized methods (LASSO, elastic net) introduce greater **coefficient shrinkage bias** than stepwise approaches?
4. Does post-selection OLS refitting produce well-calibrated **95% CI coverage**, or does empirical coverage fall below the nominal 0.95 level?

---

## Simulation Design

| Factor | Levels |
|---|---|
| Sample size (N) | 250, 500 |
| Predictor correlation (ρ) | 0, 0.35, 0.70 |
| Replications (B) | 1,000 per condition |
| **Total datasets** | **6,000** |

### Data-Generating Model

- **Predictors:** p = 20, drawn from MVN(0, Σ_ρ) with exchangeable (compound symmetry) covariance
- **Signal variables (X1–X5):** β = (1/6, 2/6, 3/6, 4/6, 5/6) ≈ (0.167, 0.333, 0.500, 0.667, 0.833)
- **Null variables (X6–X20):** β = 0
- **Error:** ε ~ N(0, 1); intercept β₀ = 0

---

## Variable Selection Methods Compared

| # | Method | Key Setting |
|---|---|---|
| 1 | Backward elimination (p-value) | Removal threshold p ≥ 0.10 |
| 2 | AIC stepwise | Bidirectional, minimize AIC |
| 3 | BIC stepwise | Bidirectional, minimize BIC |
| 4 | LASSO — λ.min | 10-fold CV, error-minimizing λ |
| 5 | LASSO — λ.1se | 10-fold CV, 1-SE rule (sparser) |
| 6 | Elastic Net — λ.min | α = 0.50, 10-fold CV, λ.min |
| 7 | Elastic Net — λ.1se | α = 0.50, 10-fold CV, λ.1se |

> **Post-selection refitting:** After penalized selection (LASSO/EN), the model is refit by OLS on the selected subset to obtain unpenalized coefficient estimates and 95% CIs.

---

## Performance Metrics

| Metric | Definition |
|---|---|
| **TPR** | Proportion of signal variables (X1–X5) correctly selected, averaged across the 5 signals |
| **FPR** | Proportion of null variables (X6–X20) incorrectly selected |
| **Coefficient Bias** | Mean (β̂ − β) for X1–X5 across replications |
| **CI Coverage** | Empirical proportion of 95% CIs that contain the true β |

All metrics are averaged across B = 1,000 replications within each of the 6 simulation conditions.

---

## Repository Structure

```
Project4/
├── README_Project4.md               ← This file
│
├── RawData/
│   └── (no real data — fully simulated)
│
├── Code/
│   ├── Interim_analysis.R           ← Simulation engine + all 7 methods
│   └── Project4_RegenerateFigures.R ← Standalone figure regeneration script
│
├── Results/
│   ├── Figures/
│   │   ├── Figure1_TruePositiveRate.png
│   │   ├── Figure2_FalsePositiveRate.png
│   │   ├── Figure3_CoefficientBias.png
│   │   ├── Figure4_CICoverage.png
│   │   ├── Figure5_TPR_vs_FPR_scatter.png
│   │   ├── Figure6_SelectionRate_Heatmap_N250_rho035.png
│   │   ├── Figure7_PerVariable_Bias_N500_rho0.png
│   │   └── Figure8_Correlation_Effect_on_TPR.png
│   └── Tables/
│       ├── Table1_SimulationParameters.{tex,png}
│       ├── Table2_TPR_FPR.{tex,png}
│       ├── Table3_Bias_Coverage.{tex,png}
│       ├── Table4_ConfusionMatrix_N250_rho035.{tex,png}
│       └── Table5_PerVariable_SelectionRates.{tex,png}
│
└── Report/
    ├── Project4_FinalReport.Rmd     ← Final reproducible report
    └── Project4_FinalReport.pdf     ← Compiled output
```

---

## How to Reproduce

### 1. Run the Simulation

```r
# In R, from the Code/ directory:
source("Interim_analysis.R")
```

This runs all 6,000 simulations, fits all 7 methods, and writes tables and figures to `Results/`.

> **Seed:** `set.seed(2026)` — results are fully reproducible.

### 2. Regenerate Figures Only

```r
source("Project4_RegenerateFigures.R")
```

### 3. Compile the Final Report

```r
rmarkdown::render("Report/Project4_FinalReport.Rmd", output_format = "pdf_document")
```

---

## R Package Dependencies

| Package | Purpose |
|---|---|
| `glmnet` | LASSO and elastic net via 10-fold CV |
| `MASS` | `mvrnorm()` for multivariate normal simulation |
| `ggplot2` | All figures |
| `knitr` / `kableExtra` | Tables in the Rmd report |
| `dplyr`, `tidyr` | Data wrangling |

Install all at once:
```r
install.packages(c("glmnet", "MASS", "ggplot2", "knitr", "kableExtra", "dplyr", "tidyr"))
```

---

## Key Findings (Summary)

- **BIC** achieves the lowest FPR (most conservative), but at the cost of lower TPR for weak signals (X1, X2) especially at N = 250.
- **LASSO λ.1se** and **Elastic Net λ.1se** behave similarly to BIC — sparse but conservative.
- **Backward (p < 0.10)** and **AIC** achieve higher TPR but admit more false positives under high correlation (ρ = 0.70).
- **CI coverage** falls below 0.95 for all methods under post-selection OLS refitting, consistent with the known post-selection inference problem; the drop is most severe for BIC and λ.1se methods that select aggressively sparse models.
- **Predictor correlation** (ρ) has a stronger negative effect on TPR than sample size alone, particularly for weak signal variables.

---

## Notes & Limitations

- All inference is post-selection: p-values and CIs should be interpreted with caution, as Type I error is not controlled at the nominal level after data-driven selection.
- The simulation uses a fixed linear model; results may not generalize to non-linear settings or models with interaction terms.
- Elastic net mixing parameter α = 0.50 was fixed per class decision; varying α was not explored.
- Lambda selection via cross-validation introduces additional Monte Carlo variability, especially at N = 250.

---

## Course Context

**Course:** Biostatistics — Applied Regression / Statistical Computing  
**Assignment:** Project 4 — Simulation Study: Variable Selection  
**Semester:** Spring 2026  

---

*All simulation code, figures, and the final report were developed iteratively using R (v4.x) with reproducible R Markdown output.*
