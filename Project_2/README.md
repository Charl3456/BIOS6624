# Peripheral Inflammation, Cognition, and Brain Structure in aMCI: Analysis Plan and Power Justification

**Course:** BIOS 6624 — Advanced Statistical Methods  
**Author:** Yunjing (Charlotte) Hong  
**Institution:** University of Colorado Anschutz Medical Campus  
**Date:** Spring 2026

---

## Overview

This project provides the **data analysis plan and sample size justification** for Aims 1 and 2 of a longitudinal study examining innate immune mechanisms of cognitive decline in amnestic mild cognitive impairment (aMCI). Data are collected at baseline and one-year follow-up at the Rocky Mountain Alzheimer's Disease Center (RMADC).

The primary research questions are:

> **Aim 1:** Are peripheral inflammatory markers (IL-6, TNF-α, MCP-1, Eotaxin-1) longitudinally associated with cognitive decline and cortical thinning in aMCI?

> **Aim 2:** Does the interaction between peripheral inflammation and amyloid deposition (Florbetapir-PET SUVR) predict clinical progression beyond either factor alone?

---

## Outcomes

Six primary outcomes are analyzed separately as one-year change scores (follow-up − baseline):

| # | Outcome | Measure | Variable |
|---|---------|---------|----------|
| 1 | CVLT Delayed Recall | 20-min delay, total correct | Change in CVLT |
| 2 | CVLT Recognition | Recognition d-prime | Change in d-prime |
| 3 | Benson Figure Recall | 15-min delay, total correct | Change in Benson |
| 4 | MST Pattern Separation | Proportion correct | Change in MST |
| 5 | Story Recall | 1-week delay, total correct | Change in Story Recall |
| 6 | AD-Signature Cortical Thickness | FreeSurfer 5.1 mean ROI | Change in cortical thickness |

---

## Methods

- **Aim 1:** Multivariable linear regression for each marker × outcome pair (4 markers × 6 outcomes = 24 models). Outcome is one-year change score with baseline outcome as covariate. Adjusted for age, sex, APOE ε4, BMI, hypercholesterolemia, NSAID use, immune-related conditions, and diagnostic group (aMCI vs. HC).
- **Aim 2:** Extends Aim 1 by adding amyloid SUVR and its interaction with each cytokine. Primary test is the interaction coefficient. Secondary analysis dichotomizes amyloid at SUVR ≥ 1.10.
- **Multiple comparisons:** FDR (Benjamini-Hochberg, *q* < 0.05) across four markers within each outcome; Bonferroni (α/4 = 0.0125) as sensitivity analysis.
- **Power (Aim 1):** `pwr.r.test()` — bivariate correlation approximation. Minimum detectable |*r*| = 0.21 (α = 0.05) and 0.25 (α = 0.0125) at 80% power, N = 175.
- **Power (Aim 2):** `pwr.f2.test()` — F-test for interaction R² increment. Minimum detectable *f*² = 0.048 (α = 0.05) and 0.068 (α = 0.0125) at 80% power, N = 175.

---

## Repository Structure

```
Project_2/
├── Code/
│   ├── Project2_PowerAnalysis.R    # Standalone power analysis script
│   └── Project2_Report.Rmd         # Report source (knit to produce final PDF)
├── DataRaw/
│   └── PrelimData.csv              # Preliminary data (n = 30)
├── Results/
│   ├── Aim1_PowerTable.csv         # Aim 1 power table
│   ├── Aim2_PowerTable.csv         # Aim 2 power table
│   └── Project2_PowerCurves.pdf    # Power curve figure (Figure 1)
└── README.md
```

**Absolute paths (from repository root):**

| File | Path |
|------|------|
| Power analysis script | `BIOS6624/Project_2/Code/Project2_PowerAnalysis.R` |
| Report source | `BIOS6624/Project_2/Code/Project2_Report.Rmd` |
| Preliminary data | `BIOS6624/Project_2/DataRaw/PrelimData.csv` |
| Aim 1 power table | `BIOS6624/Project_2/Results/Aim1_PowerTable.csv` |
| Aim 2 power table | `BIOS6624/Project_2/Results/Aim2_PowerTable.csv` |
| Power curves figure | `BIOS6624/Project_2/Results/Project2_PowerCurves.pdf` |

---

## Preliminary Data

`PrelimData.csv` contains 30 observations with 4 variables from the investigator's pilot studies:

| Variable | Description |
|----------|-------------|
| `CVLT_CNG3` | Change in CVLT memory score (follow-up − baseline) |
| `CORT_CNG3` | Change in AD-signature cortical thickness (follow-up − baseline) |
| `IL_6` | Baseline IL-6 level (pg/mL) |
| `MCP_1` | Baseline MCP-1 level (pg/mL) |

Key preliminary correlations:

| Pair | *r* | *p* |
|------|-----|-----|
| IL-6 vs. Memory Change | −0.259 | 0.168 |
| IL-6 vs. Cortical Thickness Change | −0.599 | < 0.001 |
| MCP-1 vs. Memory Change | −0.318 | 0.086 |
| MCP-1 vs. Cortical Thickness Change | −0.685 | < 0.001 |
| IL-6 vs. MCP-1 (inter-marker) | 0.934 | < 0.001 |

---

## How to Reproduce

1. Clone this repository.
2. Ensure `PrelimData.csv` is in `Project_2/DataRaw/`.
3. Open R and set working directory to `Project_2/Code/`.
4. Run `Project2_PowerAnalysis.R` to generate power tables and figure in `Project_2/Results/`.
5. Knit `Project2_Report.Rmd` to produce the final report PDF.

---

## Software

| Tool | Version |
|------|---------|
| R | 4.4.1 |
| Key packages | `pwr` (`pwr.r.test()`, `pwr.f2.test()`), `knitr`, `kableExtra` |

---

## Data Source

Preliminary data provided by the investigator (Dr. Bettcher) from ongoing studies at the Rocky Mountain Alzheimer's Disease Center (RMADC), University of Colorado Anschutz Medical Campus. Accessed via BIOS 6624 course materials.
