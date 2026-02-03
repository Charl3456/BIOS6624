# Project 0: SPIT Booklet Saliva Collection Study

**Course:** BIOS 6624 - Advanced Biostatistics  
**Author:** Hong Yunjing  
**Semester:** Spring 2026

---

## Overview

Analysis of cortisol and DHEA diurnal patterns using a novel saliva collection device (SPIT booklet). The study evaluates agreement between self-reported and electronically monitored sampling times, protocol adherence, and hormone patterns over time in 31 healthy subjects.

---

## Directory Structure

```
Project 0/
├── Code/
│   └── DataAnalysis_FINAL.Rmd       # Main analysis file
├── DataRaw/
│   └── Project0_Clean_v2.csv        # Original dataset
├── DataProcessed/
│   └── Project0_Processed.csv       # Cleaned data
├── Outputs/                          # Detailed analysis results (CSVs)
└── Reports/
    ├── ReportPlots/                  # Figures (PNG)
    └── ReportTables/                 # Tables (CSV/HTML/LaTeX)
```

---

## Quick Start

### Requirements
```r
# Install required packages
install.packages(c("tidyverse", "lubridate", "knitr", "kableExtra", 
                   "gtsummary", "gt", "lme4", "lmerTest", "broom.mixed"))
```
### Results
- **Report materials**: `Reports/ReportPlots/` and `Reports/ReportTables/`
- **Detailed outputs**: `Outputs/`
- **Processed data**: `DataProcessed/`

---

## Report Deliverables

**3 Figures** (PNG format):

- Figure 1: Agreement analysis (scatter + Bland-Altman plots)

- Figure 2: Adherence distributions (histograms)

- Figure 3: Hormone patterns (4-panel plot)

**3 Tables** (CSV/HTML/LaTeX formats):

- Table 1: Descriptive statistics 

- Table 2: Protocol adherence rates

- Table 3: Hormone model results

All materials available in `Reports/` folder.

---

## Key Methods

- **Q1 (Agreement):** Pearson correlation, mixed-effects regression, Bland-Altman analysis
- **Q2 (Adherence):** Time deviations from scheduled sampling times (±7.5, ±15 min)
- **Q3 (Hormones):** Piecewise linear mixed-effects models (knot at 30 min post-waking)

**Note:** All analyses use sleep diary wake times per investigator instructions.

