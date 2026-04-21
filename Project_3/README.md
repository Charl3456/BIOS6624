# Project 3: Risk Factors and 10-Year Probability of Stroke

## Overview

Analysis of the Framingham Heart Study data to identify important risk factors of stroke and compute 10-year stroke probabilities based on different risk profiles, stratified by sex (M/F). Conducted for BIOS 6624, Spring 2026.

**Investigators:** Elizabeth Wynn  
**Analysts:** [Your Name]

## Research Question

Are candidate cardiovascular risk factors (prevalent CHD, BP medications, current smoking, total cholesterol, BMI) associated with 10-year stroke risk after adjusting for age, diabetes, and systolic blood pressure? Do these associations differ between men and women?

## Methods Summary

- **Study design:** Prospective cohort (Framingham Heart Study), baseline Period 1 participants free of prevalent stroke (n = 4,402)
- **Outcome:** Time to first stroke within 10 years (censored at 3,652 days)
- **Models:** Sex-stratified Cox proportional hazards models
- **Variable selection:** AIC-based stepwise (backward and forward) with age, diabetes, and SBP force-included; individual LRT as secondary check
- **Model assessment:** Schoenfeld residuals (PH assumption), concordance index (discrimination), predicted 10-year probabilities for 15 risk profiles

## Directory Structure

```
Project_3/
├── README.md                  <- This file
├── Code/
│   ├── Interim_analysis.R     <- Standalone R script (all analyses, table/figure generation)
│   └── Final_Report.Rmd       <- R Markdown report (knits to PDF)
├── DataRaw/
│   └── frmgham2.csv           <- Framingham Heart Study dataset (4,434 participants, 39 variables)
└── Results/
    ├── Figures/               <- KM curves, forest plot (PNG)
    └── Tables/                <- Model results, descriptive stats (PNG + LaTeX)
```

## How to Reproduce

1. Clone the repository and navigate to `Project_3/`.
2. Ensure R (>= 4.4) and the following packages are installed:
   - `survival`, `survminer`, `gtsummary`, `ggplot2`, `dplyr`, `tidyr`, `knitr`, `kableExtra`, `gridExtra`, `xtable`, `flextable`
3. To generate all tables and figures independently:
   ```r
   setwd("Code/")
   source("Interim_analysis.R")
   ```
4. To generate the final PDF report:
   - Open `Code/Final_Report.Rmd` in RStudio
   - Knit to PDF (requires TinyTeX or another LaTeX distribution: `tinytex::install_tinytex()`)

## Key Results

- **Men:** Final model includes age, diabetes (HR = 4.81), systolic BP, and current smoking (HR = 1.92). C-index = 0.797.
- **Women:** Final model includes age, diabetes, and systolic BP only. No candidate variable improved AIC. C-index = 0.796.
- All three selection methods (backward AIC, forward AIC, individual LRT) converged on the same variable sets.
- Substantial within-individual risk factor changes across examination periods support consideration of time-varying covariate analysis in future work.

## Data Dictionary

See `Framingham_Longitudinal_Data_Documentation.pdf` in the project knowledge base for full variable descriptions. Key variables:

| Variable | Description |
|----------|-------------|
| RANDID | Unique participant ID |
| SEX | 1 = Men, 2 = Women |
| AGE | Age at exam (years) |
| SYSBP | Systolic blood pressure (mmHg) |
| DIABETES | 0/1 diabetes status |
| CURSMOKE | 0/1 current smoking status |
| PREVCHD | 0/1 prevalent coronary heart disease |
| BPMEDS | 0/1 antihypertensive medication use |
| TOTCHOL | Total cholesterol (mg/dL) |
| BMI | Body mass index (kg/m2) |
| STROKE | 0/1 stroke event |
| TIMESTRK | Days to stroke or last follow-up |
| PERIOD | Examination period (1, 2, or 3) |
