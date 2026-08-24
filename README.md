# Long-term air pollution exposure and childhood ADHD prevalence

Code accompanying the manuscript:

**Bentley, R. A. and Ozeryansky, L.**  
*Long-term air pollution exposure and childhood ADHD prevalence*

This repository contains the R code used to construct the analysis datasets, perform the statistical analyses and robustness checks, and generate the figures reported in the manuscript.

## Repository structure

The analysis is organized into three scripts that should be run in order.

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
- outlier and economic-connectedness sensitivity analyses
- Moran's I tests for residual spatial autocorrelation
- spatial-error modeling using county-contiguity weights
- construction of the manuscript figures and numerical results for the tables

## Running the analysis

The scripts should be run from the project directory in the following order:

```r
source("01_build_ACS_covariates.R")
source("02_build_health_PM25_analysis_data.R")
source("03_PM2.5_ADHD_figures_and_robustness.R")
```

Scripts 01 and 02 generate intermediate RDS files that are subsequently read by Script 03.

## R packages

The analyses use the following R packages:

```r
dplyr
tidyr
purrr
readr
ggplot2
fixest
tigris
sf
spdep
spatialreg
viridis
patchwork
grid
```

Additional packages required for importing or processing the ACS source files are specified in Script 01.

## Data

The analysis combines publicly available or previously published county-level data from several sources, including:

- U.S. Census Bureau American Community Survey (ACS)
- county-level PM2.5 exposure estimates
- U.S. county-level birth data
- county-level diabetes prevalence estimates
- county-level childhood ADHD prevalence estimates

The principal input files used by the scripts are:

```text
PM_data.csv — County-level annual PM₂.₅ estimates derived from the Atmospheric Composition Analysis Group (ACAG) satellite-based gridded PM₂.₅ data. Annual 0.01° × 0.01° gridded estimates were aggregated to U.S. counties using population-weighted averages. The underlying ACAG SatPM₂.₅ data are distributed under CC BY 4.0. See the manuscript Methods for details and citations.

Births_2007-2024.csv — County-level natality data obtained from the CDC WONDER Natality database, using the Infant Birth Weight (12-category) table for 2007–2024. Low birth weight is defined as birth weight below 2,500 g. CDC WONDER suppresses small county-level counts in accordance with NCHS confidentiality requirements. Users of these data remain subject to the CDC WONDER data-use restrictions.

Diabetes_2004-2021.csv — Annual county-level age-adjusted estimates of diagnosed diabetes prevalence among adults aged 20 years or older, obtained from the CDC United States Diabetes Surveillance System (USDSS). County estimates are based on BRFSS and U.S. Census population data and were generated using small-area estimation methods. This file preserves the USDSS estimates used in the present analysis; CDC has subsequently revised its county-level estimation methodology and historical estimates.

ADHD_2016–2018.csv — County-level childhood ADHD prevalence estimates from Zgodic et al. (2023), County-Level Prevalence Estimates of ADHD in Children in the United States, Annals of Epidemiology 79:56–64. The estimates were generated using small-area estimation applied to the 2016–2018 National Survey of Children’s Health for children aged 5–17 years. The data were obtained from the article’s supplementary file 1-s2.0-S1047279723000066-mmc2.docx; the prevalence estimate and 95% confidence interval are reported for each county. The CSV used here contains the county-level estimates extracted from that supplementary table.

full_covariates_paper.csv — County-level covariate file retained only to reproduce the covariate-adjusted PM₂.₅, low-birth-weight, and diabetes residual maps in Figure 4. The revised longitudinal robustness analyses use contemporaneous annual ACS measures of household income, educational attainment, median age, poverty, unemployment, and health-insurance coverage constructed by Script 01; the primary ADHD analysis uses 2018 ACS 5-year estimates.
```



Script 01 additionally uses annual and 5-year ACS source files organized into the directories specified in that script.

Detailed provenance and access information for the input datasets will be provided in the manuscript and repository documentation. Data derived from external sources remain subject to the terms and conditions of their original providers.

## Reproducibility

The code is organized so that the analysis datasets can be reconstructed from the input data and the statistical analyses can then be reproduced by running Scripts 01–03 sequentially.

Script 03 creates the statistical results underlying the manuscript tables and generates the manuscript figures in R. Figure files themselves are not automatically exported by the script.

## Code availability

The archived version of this repository associated with the published manuscript will be deposited in Zenodo and assigned a DOI.

## License

License information will be added prior to archival release.
