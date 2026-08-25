# Long-term air pollution exposure and childhood ADHD prevalence

Code accompanying the manuscript:

**Bentley, R. A. and L. Ozeryansky (2026)**  *Long-term air pollution exposure and childhood ADHD prevalence* in *Scientific Reports*

This repository contains the R code used to construct the analysis datasets, perform the statistical analyses and robustness checks, and generate the figures reported in the manuscript.

## Repository structure

The analysis is organized into three scripts. Script 01 documents construction of the ACS covariates from the original Census source files; Scripts 02 and 03 reproduce the analyses using the processed ACS covariates included in this repository.

### `01_build_ACS_covariates.R`

Constructs county-level socioeconomic and demographic covariates from the U.S. Census Bureau American Community Survey (ACS).

The script creates annual county-level covariates for the longitudinal low-birth-weight and diabetes analyses and 5-year ACS estimates for the childhood ADHD analysis.

### `02_build_health_PM25_analysis_data.R`

Constructs the analysis datasets by processing and merging:

- county-level PM2.5 exposure
- low birth weight
- diabetes prevalence
- childhood ADHD prevalence
- county-level covariates

Intermediate analysis datasets are saved as RDS and CSV files for use in the subsequent analyses.

### `03_PM2.5_ADHD_figures_and_robustness.R`

Performs the statistical analyses and robustness checks reported in the revised manuscript, including:

- county and year fixed-effects models for low birth weight and diabetes
- county-specific linear temporal trends
- adjustment for contemporaneous annual ACS covariates
- lagged fixed-effects analyses
- childhood ADHD models using ACS 5-year covariates and state fixed effects
- alternative PM2.5 exposure windows
- outlier sensitivity analyses
- Moran's I tests for residual spatial autocorrelation
- spatial-error modeling using county-contiguity weights
- construction of the manuscript figures and numerical results for the tables

## Running the analysis

To reproduce the analyses using the data and processed ACS covariates included in this repository, run:

```r
source("02_build_health_PM25_analysis_data.R")
source("03_PM2.5_ADHD_figures_and_robustness.R")
```

Script 02 generates the intermediate health and PM2.5 analysis datasets used by Script 03.

Script `01_build_ACS_covariates.R` documents the upstream construction of the ACS covariate panels from the original U.S. Census Bureau source files. The processed outputs of Script 01 are included in the repository, so Script 01 does not need to be rerun to reproduce the manuscript analyses.

## Data files

The repository includes the input datasets required to reproduce the analyses.

- `PM_data.csv` — County-level annual PM2.5 estimates derived from the Atmospheric Composition Analysis Group (ACAG) satellite-based gridded PM2.5 data. Annual 0.01° × 0.01° gridded estimates were aggregated to U.S. counties using population-weighted averages. The underlying ACAG SatPM2.5 data are distributed under CC BY 4.0. See the manuscript Methods for details and citations.

- `Births_2007-2024.csv` — County-level natality data obtained from the CDC WONDER Natality database, using the Infant Birth Weight (12-category) table for 2007–2024. Low birth weight is defined as birth weight below 2,500 g. CDC WONDER suppresses small county-level counts in accordance with NCHS confidentiality requirements. Users of these data remain subject to the CDC WONDER data-use restrictions.

- `Diabetes_2004-2021.csv` — Annual county-level age-adjusted estimates of diagnosed diabetes prevalence among adults aged 20 years or older, obtained from the CDC United States Diabetes Surveillance System (USDSS). These estimates are based on Behavioral Risk Factor Surveillance System survey responses and U.S. Census population data using small-area estimation methods.

- `ADHD_2016–2018.csv` — County-level childhood ADHD prevalence estimates from Zgodic et al. (2023), derived from the 2016–2018 National Survey of Children's Health using small-area estimation. Data were obtained from the article's supplementary file `1-s2.0-S1047279723000066-mmc2.docx`.

- `full_covariates_paper.csv` — County-level covariates retained to reproduce the covariate-adjusted PM2.5, low-birth-weight, and diabetes residual maps in Figure 4. The revised longitudinal robustness analyses instead use contemporaneous annual ACS covariates, while the primary ADHD analysis uses 2018 ACS 5-year estimates.

## ACS covariates

Socioeconomic and demographic covariates were obtained from the U.S. Census Bureau American Community Survey (ACS).

`01_build_ACS_covariates.R` contains the code used to construct the ACS covariate panels from the original Census ACS downloads. Because the raw ACS downloads consist of numerous publicly available source files, they are not duplicated in this repository.

For reproducibility of the downstream analyses, the processed outputs generated by Script 01 are included directly:

- `ACS_county_covariates_annual.rds`
- `ACS_county_covariates_annual.csv`
- `ACS_county_covariates_5year.rds`
- `ACS_county_covariates_5year.csv`

The annual panel contains contemporaneous county-level measures used in the low-birth-weight and diabetes robustness analyses, including median household income, median age, educational attainment, poverty, unemployment, and health-insurance coverage.

The 5-year panel provides the ACS estimates used for the cross-sectional ADHD analysis. The primary ADHD model uses the 2018 ACS 5-year estimates, representing 2014–2018.

The `.rds` files are read directly by the downstream R scripts. The equivalent `.csv` files are provided for transparency and convenient inspection.


## Code availability

The archived version of this repository associated with the published manuscript will be deposited in Zenodo and assigned a DOI.
