# segInc

**segInc** is an experimental R package for estimating and comparing **cancer incidence trends** in **CI5plus** using **segmented (piecewise) regression** via the [`segmented`](https://cran.r-project.org/package=segmented) R package, and selecting breakpoint models using a **weighted BIC (wBIC)** approach.

This repository is also a **methods case study** intended to align (as closely as practical) with analyses using the **NCI Joinpoint Regression Program**, as used in **Sung et al. (2025, *The Lancet*)** for CI5plus-based trend analyses—while documenting where results may differ due to different optimisation and model-selection strategies.

> Status: research / prototyping. API and outputs may change.

---

## Case study scope

### Cancer sites included
This case study focuses on:
- **Colorectal**
- **Upper GI group**:
  - oesophagus
  - stomach
  - pancreas
  - gallbladder
  - liver
- **Breast**
- **Kidney**

> Practical note: CI5plus extracts may encode sites as ICD-10, ICD-O, or CI5plus site groupings depending on how you export. You’ll need a consistent site mapping layer in your preprocessing.

### Age bands (primary + comparator)
Primary age bands:
- **15–39** (“early”)
- **40–49** (“early”)
- **50+** (“late”)

Comparator age bands to mirror the colorectal comparison in **Sung et al. (2025)**:
- **15–49** (“early”)
- **50–74** (“late”)

---

## Why segmented regression (and how it relates to Joinpoint)

The **NCI Joinpoint** software fits piecewise log-linear trends and typically selects joinpoints via permutation tests and related criteria.

In contrast, this project explores:
- fitting log-linear **Poisson** models (with person-time/population offset),
- estimating breakpoints using `segmented::segmented()`,
- selecting among candidate numbers of breakpoints using a **wBIC** criterion.

These approaches often agree qualitatively, but can differ because of:
- different breakpoint search/initialisation strategies,
- different model-selection criteria (permutation tests vs information criteria),
- constraints (minimum segment length, max joinpoints),
- overdispersion handling.

**Goal:** make the R workflow transparent and reproducible, and quantify sensitivity of breakpoint placement and trend estimates.

---

## Data: CI5plus

CI5plus is produced by the **International Agency for Research on Cancer (IARC)**.

1. Obtain CI5plus data through the appropriate channels.
2. Do **not** commit restricted CI5plus data to this repository unless you have explicit permission.
3. Ensure you comply with CI5plus terms of use and citation guidance.

### Expected analysis-ready structure (minimum)
Your derived dataset should be aggregatable to something like:

- `year` (integer)
- `age` or `age_group` (needed to create the study age bands)
- `cases` (incident counts)
- `py` (person-years / population-at-risk measure)
- `site` (or site code)
- optional: `sex`, `registry`, `country`, `quality_flags`, etc.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("cmclean5/segInc")
