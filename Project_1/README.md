# Treatment Response and Hard Drug Use in HIV-Infected Men on HAART

**Course:** BIOS 6624 — Advanced Statistical Methods  
**Author:** Yunjing (Charlotte) Hong  
**Institution:** University of Colorado Anschutz Medical Campus  
**Date:** Spring 2026

---

## Overview

This project is a secondary data analysis of the **Multicenter AIDS Cohort Study (MACS)**, a prospective cohort study of HIV-1 infection in homosexual and bisexual men across four major U.S. cities.

The primary research question is:

> **Does treatment response 2 years after initiating HAART differ between men who reported hard drug use (heroin, opiates, or injection drugs) at baseline compared to those who did not?**

A secondary question examines whether any observed differences in outcomes can be explained by differences in ART adherence between the two groups.

---

## Outcomes

Four treatment response outcomes are examined at Year 2:

| Outcome | Variable | Description |
|---|---|---|
| Viral Load | `log10VLOAD` | Log10-transformed HIV copies per mL |
| CD4+ T Cell Count | `LEU3N` | Immunologic health marker (cells/mL) |
| Physical Quality of Life | `AGG_PHYS` | SF-36 physical component score (0–100) |
| Mental Quality of Life | `AGG_MENT` | SF-36 mental component score (0–100) |

---

## Methods

- **Frequentist:** Multivariable linear regression for each outcome, adjusting for baseline outcome value, age, race (White vs. Other), BMI, education (college degree+ vs. not), and smoking status
- **Bayesian:** Same regression structure with vague priors ($\beta_j \sim N(0, 100^2)$) and a skeptical sensitivity prior ($\beta_1 \sim N(0, 10^2)$) fitted via Stan (NUTS sampler)
- **Adherence analysis:** Coefficient attenuation approach — comparing the hard drug use coefficient before and after adjusting for ART adherence
- **Model comparison:** LOO-IC (Bayesian models with vs. without hard drug use)
- **Prior sensitivity:** Vague vs. skeptical prior posteriors compared

---

## Repository Structure

```
Project_1/
├── Code/
│   └── hiv_haart_analysis.Rmd     # Main analysis file
├── Raw_data/                       # ⚠️ Not tracked (see .gitignore)
│   └── hiv_6624_final (3).csv
├── Figures/                        # Auto-generated on knit
│   ├── fig1_adherence_by_drug_use.png
│   ├── fig2_freq_bayes_comparison.png
│   ├── figA1_trajectories.png
│   ├── figA2_outcome_distributions.png
│   ├── figA3_freq_diagnostics.png
│   ├── figA4a_mcmc_trace.png
│   └── figA4b_mcmc_density.png
├── Tables/                         # Auto-generated on knit
│   ├── tableA1_missing_data.tex
│   ├── table_freq_results.tex
│   ├── table_adherence_mediation.tex
│   ├── table_bayesian_posterior.tex
│   ├── table_loo_comparison.tex
│   ├── table_prior_sensitivity.tex
│   └── tableA3_convergence.tex
├── .gitignore
└── README.md
```

> Raw data are excluded from version control per data sharing restrictions. Contact the author for access questions.

---

## Key Findings

| Outcome | Frequentist β | p-value | Direction |
|---|---|---|---|
| Log10 Viral Load | ~ −0.04 | 0.85 | Not significant |
| CD4+ T Cell Count | ~ −168 | < 0.001 | Hard drug users lower ↓ |
| Physical QoL | ~ −3.3 | 0.014 | Hard drug users lower ↓ |
| Mental QoL | ~ 0.03 | 0.99 | Not significant |

Bayesian results under vague priors were concordant with frequentist estimates. ART adherence differed significantly between groups and partially explained outcome differences for some outcomes.

---

## Software

| Tool | Version |
|---|---|
| R | 4.4.1 |
| Stan (via cmdstanr) | latest |
| Key packages | tidyverse, cmdstanr, posterior, bayesplot, loo, HDInterval, tableone, kableExtra |

---

## Data Source

Multicenter AIDS Cohort Study (MACS). Data accessed via BIOS 6624 course materials, University of Colorado Anschutz. Not publicly redistributable.
