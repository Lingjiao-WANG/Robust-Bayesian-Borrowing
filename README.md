# Reproducibility Code

## Manuscript Title

**Robust Bayesian Borrowing with Local Control of Type I Error and Power**

## Authors

Lingjiao Wang, Astrid Jullion, Xuan Zhu, Marc Vandemeulebroecke, Laurens Sluijterman, David Jesse, Francois Mercier, Johannes Vilsmeier, Teresa Engelbrecht, Martin Posch, Sarah Zohar and Nigel Stallard

## Computing Environment

The numerical results were obtained using the following configuration:

* **R version:** 4.5.1 (2025-06-13 ucrt)
* **Platform:** x86_64-w64-mingw32/x64
* **Operating environment:** Windows PC
* **Processor:** Intel(R) Core(TM) i5-1145G7 CPU
* **RAM:** 8.00 GB
* **Number of cores used:** 8

## Execution

To reproduce the numerical results, tables, and figures presented in the manuscript and Supplementary Material, please follow the instructions below.

## Main Manuscript

### 1. Functions

* `case_study_funs.R`
    * Contains the functions required to reproduce the numerical results.

### 2. Calibrated Analysis Priors — Section 5.2

* **Table 1:** `Table_1_calibrated_priors.R`
* **Figure 1:** `Figure_1_analysis_priors.R`

### 3. Type I Error, Power, and Calibrated Power Gain — Section 5.3

* **Table 2:** `Table_2_TIE_power.R`
* **Figure 2:** `Figure_2_TIE_power.R`

### 4. Effective Historical Sample Size — Section 5.4

* **Figure 3:** `Figure_3_EHSS.R`

### 5. Sample Size Determination — Section 5.5

* **Figure 4:** `Figure_4_sample_size.R`

## Supplementary Material

### 1. Numerical Calibration of the Borrowing Parameters — Section B

* **Figure A1:** `Figure_A1_cal_MAP_rMAP.R`

### 2. Calibration of the Widest Feasible Discrepancy Region for the APP — Section C

* **Figure A2:** `Figure_A2_TIE_pow_APP.R`

### 3. Average Type I Error and Average Power — Section D

* **Tables A2–A3:** `Tables_A2_A3_ave_TIE_power.R`
* **Figure A3:** `Figure_A3_design_priors.R`

### 4. Sample Size Determination — Section E

* **Figures A4–A6:** `Figures_A4_A6_sample_sizes.R`

## Notes

1. The `Rdatas` folder contains the saved data and results required to reproduce the tables and figures.

2. The following R scripts can also be used to reproduce the underlying numerical results:

    * `case_study_cal_analysis_priors.R`
        * Computes the calibrated borrowing parameters and posterior probability thresholds for the analysis priors.

    * `case_study_TIE_power.R`
        * Computes the Type I error (TIE), power, and related operating characteristics under the calibrated analysis priors.

    * `case_study_EHSS.R`
        * Computes the effective historical sample size (EHSS) under the analysis priors.

    * `case_study_ave_TIE_power.R`
        * Computes the average TIE and power under the considered design priors.

3. The scripts for reproducing the tables and figures use the saved results in the `Rdatas` folder. Re-running the underlying numerical calculations is therefore not required to reproduce the reported tables and figures.




