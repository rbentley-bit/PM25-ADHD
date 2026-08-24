# ============================================================
# 02_build_health_PM25_analysis_data.R
#
# Construct the health and PM2.5 datasets used in the paper.
#
# Statistical models, robustness analyses, maps, and figures
# are handled in 03_PM2.5_ADHD_robustness.R.
#
# Inputs:
#   PM_data.csv
#   Births_2007-2024.csv
#   Diabetes_2004-2021.csv
#   ADHD_2016–2018.csv
#   full_covariates_paper.csv
#
# Outputs:
#   PM25_county_year.rds
#   LBW_county_year.rds
#   LBW_PM25_county_year.rds
#   Diabetes_county_year.rds
#   Diabetes_PM25_county_year.rds
#   ADHD_PM25_county.rds
#
# CSV versions are also written for inspection.
# ============================================================


# ============================================================
# 0. Setup
# ============================================================

library(dplyr)
library(tidyr)
library(readr)


# ============================================================
# 1. Helper functions
# ============================================================

# Standardize county FIPS as both:
#
#   FIPS     = numeric form used in several original datasets
#   FIPS_chr = five-character Census form, e.g. "01003"
#
# Keeping both makes later joins explicit and avoids repeatedly
# converting FIPS in downstream scripts.

add_fips_chr <- function(data) {
  
  data %>%
    mutate(
      FIPS = as.numeric(FIPS),
      FIPS_chr = sprintf(
        "%05d",
        FIPS
      )
    )
}


# Safe first non-missing value

first_nonmissing <- function(x) {
  
  x <- x[
    !is.na(x)
  ]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  x[1]
}


# ============================================================
# 2. Load source datasets
# ============================================================

pm <- read_csv(
  "PM_data.csv",
  show_col_types = FALSE
)

births <- read_csv(
  "Births_2007-2024.csv",
  show_col_types = FALSE
)

diabetes_raw <- read_csv(
  "Diabetes_2004-2021.csv",
  show_col_types = FALSE
)

adhd_raw <- read_csv(
  "ADHD_2016–2018.csv",
  show_col_types = FALSE
)

county_covariates_original <- read_csv(
  "full_covariates_paper.csv",
  show_col_types = FALSE
)


# Standardize FIPS in source datasets

pm <- add_fips_chr(
  pm
)

births <- add_fips_chr(
  births
)

diabetes_raw <- add_fips_chr(
  diabetes_raw
)

adhd_raw <- add_fips_chr(
  adhd_raw
)

county_covariates_original <- add_fips_chr(
  county_covariates_original
)


# ============================================================
# PART I
# PM2.5 COUNTY-YEAR DATA
# ============================================================


# ============================================================
# 3. Extract PM2.5 observations
# ============================================================

pm25 <- pm %>%
  filter(
    pollutant == "pm25"
  ) %>%
  rename(
    PM25 = pred_wght
  ) %>%
  arrange(
    FIPS,
    Year
  )


# Inspect county-year uniqueness

pm25 %>%
  count(
    FIPS,
    Year
  ) %>%
  filter(
    n > 1
  )


# ============================================================
# 4. Save PM2.5 county-year data
# ============================================================

saveRDS(
  pm25,
  "PM25_county_year.rds"
)

write_csv(
  pm25,
  "PM25_county_year.csv"
)



# ============================================================
# PART II
# LOW BIRTH WEIGHT
# ============================================================


# ============================================================
# 5. Construct county-year low-birth-weight outcomes
# ============================================================

lbw <- births %>%
  mutate(
    Ave_Birth_Weight =
      as.numeric(Ave_Birth_Weight),
    Ave_Age_Mother =
      as.numeric(Ave_Age_Mother),
    Birth_Weight_12_Code_num =
      as.numeric(Birth_Weight_12_Code),
    low_bw =
      Birth_Weight_12_Code_num <= 5
  ) %>%
  
  group_by(
    County,
    FIPS,
    FIPS_chr,
    Year
  ) %>%
  
  summarise(
    
    total_births =
      sum(
        Births,
        na.rm = TRUE
      ),
    
    low_bw_births =
      sum(
        Births[low_bw],
        na.rm = TRUE
      ),
    
    Ave_Age_Mother =
      first_nonmissing(
        Ave_Age_Mother
      ),
    
    Ave_Birth_Weight =
      sum(
        Births *
          Ave_Birth_Weight,
        na.rm = TRUE
      ) /
      sum(
        Births,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    low_bw_rate =
      low_bw_births /
      total_births,
    
    low_bw_pct =
      100 *
      low_bw_rate
    
  ) %>%
  
  filter(
    !is.na(FIPS)
  ) %>%
  
  arrange(
    FIPS,
    Year
  )


# ============================================================
# 6. Merge low birth weight with PM2.5
# ============================================================

lbw_pm25 <- pm25 %>%
  
  left_join(
    lbw,
    by = c(
      "FIPS",
      "FIPS_chr",
      "Year"
    )
  ) %>%
  
  filter(
    !is.na(
      low_bw_rate
    )
  ) %>%
  
  arrange(
    FIPS,
    Year
  )


# ============================================================
# 7. LBW quality checks
# ============================================================

dim(
  lbw
)

dim(
  lbw_pm25
)


# County-year uniqueness

lbw %>%
  count(
    FIPS,
    Year
  ) %>%
  filter(
    n > 1
  )


# Mean low-birth-weight prevalence by year

lbw %>%
  group_by(
    Year
  ) %>%
  summarise(
    n_counties =
      n(),
    
    mean_lbw_pct =
      mean(
        low_bw_pct,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# Plausible ranges

summary(
  lbw %>%
    select(
      total_births,
      low_bw_pct,
      Ave_Age_Mother,
      Ave_Birth_Weight
    )
)


# ============================================================
# 8. Save LBW datasets
# ============================================================

saveRDS(
  lbw,
  "LBW_county_year.rds"
)

write_csv(
  lbw,
  "LBW_county_year.csv"
)


saveRDS(
  lbw_pm25,
  "LBW_PM25_county_year.rds"
)

write_csv(
  lbw_pm25,
  "LBW_PM25_county_year.csv"
)



# ============================================================
# PART III
# DIABETES
# ============================================================


# ============================================================
# 9. Reshape diabetes prevalence from wide to long
# ============================================================

diabetes_long <- diabetes_raw %>%
  
  mutate(
    across(
      starts_with(
        "Diabetes_"
      ),
      ~ parse_number(
        as.character(.x)
      )
    )
  ) %>%
  
  pivot_longer(
    
    cols =
      starts_with(
        "Diabetes_"
      ),
    
    names_to =
      "Year",
    
    values_to =
      "Diabetes"
  ) %>%
  
  mutate(
    
    Year =
      as.numeric(
        sub(
          "Diabetes_",
          "",
          Year
        )
      )
    
  ) %>%
  
  arrange(
    FIPS,
    Year
  )


# ============================================================
# 10. Keep clean diabetes county-year dataset
# ============================================================

# Preserve geographic identifiers from the diabetes source
# where available.

diabetes <- diabetes_long %>%
  
  select(
    any_of(
      c(
        "FIPS",
        "FIPS_chr",
        "County",
        "State",
        "Year",
        "Diabetes"
      )
    )
  ) %>%
  
  filter(
    !is.na(
      Diabetes
    )
  ) %>%
  
  arrange(
    FIPS,
    Year
  )


# ============================================================
# 11. Merge diabetes with PM2.5
# ============================================================

diabetes_pm25 <- pm25 %>%
  
  select(
    FIPS,
    FIPS_chr,
    Year,
    PM25,
    everything()
  ) %>%
  
  left_join(
    diabetes %>%
      select(
        -any_of(
          c(
            "County",
            "State"
          )
        )
      ),
    by = c(
      "FIPS",
      "FIPS_chr",
      "Year"
    )
  ) %>%
  
  filter(
    !is.na(
      Diabetes
    )
  ) %>%
  
  arrange(
    FIPS,
    Year
  )


# ============================================================
# 12. Diabetes quality checks
# ============================================================

dim(
  diabetes
)

dim(
  diabetes_pm25
)


# County-year uniqueness

diabetes %>%
  count(
    FIPS,
    Year
  ) %>%
  filter(
    n > 1
  )


# Coverage by year

diabetes %>%
  count(
    Year
  )


summary(
  diabetes$Diabetes
)


# ============================================================
# 13. Save diabetes datasets
# ============================================================

saveRDS(
  diabetes,
  "Diabetes_county_year.rds"
)

write_csv(
  diabetes,
  "Diabetes_county_year.csv"
)


saveRDS(
  diabetes_pm25,
  "Diabetes_PM25_county_year.rds"
)

write_csv(
  diabetes_pm25,
  "Diabetes_PM25_county_year.csv"
)



# ============================================================
# PART IV
# ADHD AND LONG-TERM PM2.5 EXPOSURE WINDOWS
# ============================================================


# ============================================================
# 14. Helper for county-level PM2.5 averages
# ============================================================

make_pm25_window <- function(
    start_year,
    end_year,
    variable_name
) {
  
  out <- pm25 %>%
    
    filter(
      Year >= start_year,
      Year <= end_year
    ) %>%
    
    group_by(
      FIPS,
      FIPS_chr
    ) %>%
    
    summarise(
      
      value =
        mean(
          PM25,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  names(out)[
    names(out) == "value"
  ] <- variable_name
  
  out
}


# ============================================================
# 15. Construct all ADHD PM2.5 exposure windows
# ============================================================

# Primary long-term exposure measure

pm25_2005_2015 <- make_pm25_window(
  2005,
  2015,
  "pm25_2005_2015"
)


# Alternative exposure windows used for sensitivity analysis

pm25_2005_2006 <- make_pm25_window(
  2005,
  2006,
  "pm25_2005_2006"
)

pm25_2007_2009 <- make_pm25_window(
  2007,
  2009,
  "pm25_2007_2009"
)

pm25_2010_2012 <- make_pm25_window(
  2010,
  2012,
  "pm25_2010_2012"
)

pm25_2013_2015 <- make_pm25_window(
  2013,
  2015,
  "pm25_2013_2015"
)


# ============================================================
# 16. Construct one canonical ADHD analysis dataset
# ============================================================

adhd_data <- adhd_raw %>%
  
  left_join(
    pm25_2005_2015,
    by = c(
      "FIPS",
      "FIPS_chr"
    )
  ) %>%
  
  left_join(
    pm25_2005_2006,
    by = c(
      "FIPS",
      "FIPS_chr"
    )
  ) %>%
  
  left_join(
    pm25_2007_2009,
    by = c(
      "FIPS",
      "FIPS_chr"
    )
  ) %>%
  
  left_join(
    pm25_2010_2012,
    by = c(
      "FIPS",
      "FIPS_chr"
    )
  ) %>%
  
  left_join(
    pm25_2013_2015,
    by = c(
      "FIPS",
      "FIPS_chr"
    )
  ) %>%
  
  # Retain economic connectedness and state identifier from
  # the original county covariate dataset for sensitivity
  # analyses in Script 03.
  left_join(
    
    county_covariates_original %>%
      
      select(
        FIPS,
        FIPS_chr,
        any_of(
          c(
            "State",
            "econ_connect"
          )
        )
      ),
    
    by = c(
      "FIPS",
      "FIPS_chr"
    ),
    
    suffix = c(
      "",
      "_covariates"
    )
  )


# ------------------------------------------------------------
# Resolve State if both source files provide it
# ------------------------------------------------------------

# Prefer the State value already contained in the ADHD source.
# If it is absent there, use the value from the county
# covariates file.

if (
  "State_covariates" %in%
  names(
    adhd_data
  )
) {
  
  if (
    "State" %in%
    names(
      adhd_data
    )
  ) {
    
    adhd_data <- adhd_data %>%
      mutate(
        State =
          coalesce(
            State,
            State_covariates
          )
      ) %>%
      select(
        -State_covariates
      )
    
  } else {
    
    adhd_data <- adhd_data %>%
      rename(
        State =
          State_covariates
      )
  }
}


# Backward-compatible alias used in some earlier scripts

adhd_data <- adhd_data %>%
  mutate(
    pm25_mean =
      pm25_2005_2015
  )


# ============================================================
# 17. ADHD quality checks
# ============================================================

dim(
  adhd_data
)


# County uniqueness

adhd_data %>%
  count(
    FIPS
  ) %>%
  filter(
    n > 1
  )


# PM2.5 coverage for each exposure window

adhd_data %>%
  summarise(
    
    n =
      n(),
    
    missing_2005_2015 =
      sum(
        is.na(
          pm25_2005_2015
        )
      ),
    
    missing_2005_2006 =
      sum(
        is.na(
          pm25_2005_2006
        )
      ),
    
    missing_2007_2009 =
      sum(
        is.na(
          pm25_2007_2009
        )
      ),
    
    missing_2010_2012 =
      sum(
        is.na(
          pm25_2010_2012
        )
      ),
    
    missing_2013_2015 =
      sum(
        is.na(
          pm25_2013_2015
        )
      )
  )


summary(
  adhd_data %>%
    select(
      ADHD,
      pm25_2005_2015,
      pm25_2005_2006,
      pm25_2007_2009,
      pm25_2010_2012,
      pm25_2013_2015
    )
)


# ============================================================
# 18. Save canonical ADHD dataset
# ============================================================

saveRDS(
  adhd_data,
  "ADHD_PM25_county.rds"
)

write_csv(
  adhd_data,
  "ADHD_PM25_county.csv"
)



# ============================================================
# 19. Final summary
# ============================================================

message(
  "Health/PM2.5 data construction complete."
)

message(
  "PM2.5 panel: ",
  nrow(pm25),
  " county-year observations."
)

message(
  "LBW panel: ",
  nrow(lbw_pm25),
  " matched county-year observations."
)

message(
  "Diabetes panel: ",
  nrow(diabetes_pm25),
  " matched county-year observations."
)

message(
  "ADHD dataset: ",
  nrow(adhd_data),
  " counties."
)

message(
  "All five ADHD PM2.5 exposure measures were constructed before",
  " the final ADHD dataset was saved."
)

