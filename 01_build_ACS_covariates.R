# ============================================================
# 01_build_ACS_covariates.R
#
# Build county-level ACS covariate datasets used in the paper:
#
# 1. Annual ACS 1-year estimates
#    Used for low-birth-weight robustness analyses
#
# 2. ACS 5-year estimates
#    Used for childhood ADHD robustness analyses
#
# Outputs:
#   ACS_county_covariates_annual.csv
#   ACS_county_covariates_annual.rds
#   ACS_county_covariates_5year.csv
#   ACS_county_covariates_5year.rds
# ============================================================


# ============================================================
# 0. Setup
# ============================================================

library(dplyr)
library(purrr)
library(readr)
library(stringr)


# ============================================================
# Shared helper functions
# ============================================================

# Convert Census GEO_ID:
# 0500000US01003
# to five-digit county FIPS:
# 01003

get_fips <- function(x) {
  str_extract(x, "\\d{5}$")
}


# Safely convert ACS values to numeric

to_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}



# ============================================================
# PART I
# ACS 1-YEAR COUNTY COVARIATES
# ============================================================


# ============================================================
# 1. ACS 1-year folders
# ============================================================

# Only folders needed for the revised LBW robustness analyses

acs1_folders <- c(
  "ACS_median_income_1-year",
  "ACS_age&sex_1-year",
  "ACS_education_1-year",
  "ACS_poverty_1-year",
  "ACS_employment_1-year",
  "ACS_insurance_1-year"
)


# ============================================================
# 2. ACS 1-year helper functions
# ============================================================

# Extract year from filenames such as:
# ACSST1Y2010.S1501-Data.csv

get_year1 <- function(x) {
  as.integer(
    str_extract(
      basename(x),
      "(?<=1Y)\\d{4}"
    )
  )
}


# Find Data.csv file for a given folder and year

find_data_file1 <- function(folder, year) {
  
  files <- list.files(
    folder,
    pattern = "-Data\\.csv$",
    full.names = TRUE
  )
  
  hit <- files[
    get_year1(files) == year
  ]
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  hit[1]
}


# ============================================================
# 3. Determine ACS 1-year endpoint years available
# ============================================================

years_by_folder1 <- lapply(
  acs1_folders,
  function(folder) {
    
    files <- list.files(
      folder,
      pattern = "-Data\\.csv$",
      full.names = TRUE
    )
    
    sort(
      unique(
        get_year1(files)
      )
    )
  }
)


# Use only years represented in every required folder.
# This naturally excludes years for which standard ACS
# 1-year products are unavailable.

available_years1 <- Reduce(
  intersect,
  years_by_folder1
)

available_years1


# ============================================================
# 4. Extract one year of ACS 1-year data
# ============================================================

extract_acs1_year <- function(year) {
  
  message(
    "Processing ACS 1-year: ",
    year
  )
  
  
  # ----------------------------------------------------------
  # MEDIAN HOUSEHOLD INCOME
  # S1903
  # ----------------------------------------------------------
  
  f <- find_data_file1(
    "ACS_median_income_1-year",
    year
  )
  
  income <- read_csv(
    f,
    show_col_types = FALSE
  ) %>%
    slice(-1) %>%
    transmute(
      FIPS = get_fips(GEO_ID),
      County = NAME,
      
      median_household_income =
        to_numeric(
          S1903_C02_001E
        )
    )
  
  
  # ----------------------------------------------------------
  # AGE
  # S0101
  #
  # Table structure changed in 2017
  # ----------------------------------------------------------
  
  f <- find_data_file1(
    "ACS_age&sex_1-year",
    year
  )
  
  age_raw <- read_csv(
    f,
    show_col_types = FALSE
  ) %>%
    slice(-1)
  
  if (year <= 2016) {
    
    age <- age_raw %>%
      transmute(
        FIPS = get_fips(GEO_ID),
        
        median_age =
          to_numeric(
            S0101_C01_030E
          ),
        
        pct_age65 =
          to_numeric(
            S0101_C01_028E
          )
      )
    
  } else {
    
    age <- age_raw %>%
      transmute(
        FIPS = get_fips(GEO_ID),
        
        median_age =
          to_numeric(
            S0101_C01_032E
          ),
        
        pct_age65 =
          to_numeric(
            S0101_C02_030E
          )
      )
  }
  
  
  # ----------------------------------------------------------
  # EDUCATION
  # S1501
  #
  # Table structure changed in 2015
  # ----------------------------------------------------------
  
  f <- find_data_file1(
    "ACS_education_1-year",
    year
  )
  
  education_raw <- read_csv(
    f,
    show_col_types = FALSE
  ) %>%
    slice(-1)
  
  if (year <= 2014) {
    
    education <- education_raw %>%
      transmute(
        FIPS = get_fips(GEO_ID),
        
        pct_bachelors =
          to_numeric(
            S1501_C01_015E
          )
      )
    
  } else {
    
    education <- education_raw %>%
      transmute(
        FIPS = get_fips(GEO_ID),
        
        pct_bachelors =
          to_numeric(
            S1501_C02_015E
          )
      )
  }
  
  
  # ----------------------------------------------------------
  # POVERTY
  # S1701
  # ----------------------------------------------------------
  
  f <- find_data_file1(
    "ACS_poverty_1-year",
    year
  )
  
  poverty <- read_csv(
    f,
    show_col_types = FALSE
  ) %>%
    slice(-1) %>%
    transmute(
      FIPS = get_fips(GEO_ID),
      
      poverty_population =
        to_numeric(
          S1701_C01_001E
        ),
      
      below_poverty =
        to_numeric(
          S1701_C02_001E
        )
    )
  
  
  # ----------------------------------------------------------
  # EMPLOYMENT
  # S2301
  # ----------------------------------------------------------
  
  f <- find_data_file1(
    "ACS_employment_1-year",
    year
  )
  
  employment <- read_csv(
    f,
    show_col_types = FALSE
  ) %>%
    slice(-1) %>%
    transmute(
      FIPS = get_fips(GEO_ID),
      
      unemployment_rate =
        to_numeric(
          S2301_C04_001E
        )
    )
  
  
  # ----------------------------------------------------------
  # HEALTH INSURANCE
  # S2701
  #
  # Table definition changed after 2014
  # ----------------------------------------------------------
  
  f <- find_data_file1(
    "ACS_insurance_1-year",
    year
  )
  
  insurance_raw <- read_csv(
    f,
    show_col_types = FALSE
  ) %>%
    slice(-1)
  
  if (year <= 2014) {
    
    insurance <- insurance_raw %>%
      transmute(
        FIPS = get_fips(GEO_ID),
        
        pct_uninsured =
          to_numeric(
            S2701_C03_001E
          )
      )
    
  } else {
    
    insurance <- insurance_raw %>%
      transmute(
        FIPS = get_fips(GEO_ID),
        
        pct_uninsured =
          100 -
          to_numeric(
            S2701_C03_001E
          )
      )
  }
  
  
  # ----------------------------------------------------------
  # MERGE ONE-YEAR COVARIATES
  # ----------------------------------------------------------
  
  out <- income %>%
    
    left_join(
      age,
      by = "FIPS"
    ) %>%
    
    left_join(
      education,
      by = "FIPS"
    ) %>%
    
    left_join(
      poverty,
      by = "FIPS"
    ) %>%
    
    left_join(
      employment,
      by = "FIPS"
    ) %>%
    
    left_join(
      insurance,
      by = "FIPS"
    ) %>%
    
    mutate(
      
      Year = year,
      
      pct_poverty =
        100 *
        below_poverty /
        poverty_population
      
    ) %>%
    
    select(
      Year,
      FIPS,
      County,
      median_household_income,
      median_age,
      pct_age65,
      pct_bachelors,
      pct_poverty,
      unemployment_rate,
      pct_uninsured
    )
  
  return(out)
}


# ============================================================
# 5. Build full ACS 1-year panel
# ============================================================

acs_covariates_annual <- map_dfr(
  available_years1,
  extract_acs1_year
) %>%
  arrange(
    FIPS,
    Year
  )


# ============================================================
# 6. Quality checks: ACS 1-year
# ============================================================

dim(
  acs_covariates_annual
)

head(
  acs_covariates_annual
)


# Number of counties represented by year

acs_covariates_annual %>%
  count(Year)


# Missing values by year

acs_covariates_annual %>%
  group_by(Year) %>%
  summarise(
    
    n = n(),
    
    income_missing =
      sum(
        is.na(
          median_household_income
        )
      ),
    
    age_missing =
      sum(
        is.na(
          median_age
        )
      ),
    
    education_missing =
      sum(
        is.na(
          pct_bachelors
        )
      ),
    
    poverty_missing =
      sum(
        is.na(
          pct_poverty
        )
      ),
    
    unemployment_missing =
      sum(
        is.na(
          unemployment_rate
        )
      ),
    
    insurance_missing =
      sum(
        is.na(
          pct_uninsured
        )
      )
  )


# Plausible ranges

summary(
  acs_covariates_annual %>%
    select(
      median_household_income,
      median_age,
      pct_age65,
      pct_bachelors,
      pct_poverty,
      unemployment_rate,
      pct_uninsured
    )
)


# Impossible percentages

acs_covariates_annual %>%
  filter(
    pct_age65 < 0 |
      pct_age65 > 100 |
      pct_bachelors < 0 |
      pct_bachelors > 100 |
      pct_poverty < 0 |
      pct_poverty > 100 |
      unemployment_rate < 0 |
      unemployment_rate > 100 |
      pct_uninsured < 0 |
      pct_uninsured > 100
  )


# ============================================================
# 7. Save ACS 1-year panel
# ============================================================

write_csv(
  acs_covariates_annual,
  "ACS_county_covariates_annual.csv"
)

saveRDS(
  acs_covariates_annual,
  "ACS_county_covariates_annual.rds"
)



# ============================================================
# PART II
# ACS 5-YEAR COUNTY COVARIATES
# ============================================================


# ============================================================
# 8. ACS 5-year folders
# ============================================================

acs5_folders <- c(
  "ACS_age&sex_5-year",
  "ACS_median_income_5-year",
  "ACS_education_5-year"
)


# ============================================================
# 9. ACS 5-year helper functions
# ============================================================

# Extract endpoint year from filenames such as:
# ACSST5Y2018.S1501-Data.csv

get_year5 <- function(x) {
  
  as.integer(
    str_extract(
      basename(x),
      "(?<=5Y)\\d{4}"
    )
  )
}


# Find Data.csv file for a given folder/year

find_data_file5 <- function(
    folder,
    year
) {
  
  files <- list.files(
    folder,
    pattern = "-Data\\.csv$",
    full.names = TRUE
  )
  
  hit <- files[
    get_year5(files) == year
  ]
  
  if (length(hit) == 0) {
    return(
      NA_character_
    )
  }
  
  hit[1]
}


# Some folders contain more than one ACS table.
# Allow table ID to be specified explicitly.

find_data_file5_table <- function(
    folder,
    year,
    table_id
) {
  
  files <- list.files(
    folder,
    pattern = "-Data\\.csv$",
    full.names = TRUE
  )
  
  hit <- files[
    get_year5(files) == year &
      str_detect(
        basename(files),
        fixed(table_id)
      )
  ]
  
  if (length(hit) == 0) {
    return(
      NA_character_
    )
  }
  
  hit[1]
}


# ============================================================
# 10. Determine available ACS 5-year endpoint years
# ============================================================

years_by_folder5 <- lapply(
  acs5_folders,
  function(folder) {
    
    files <- list.files(
      folder,
      pattern = "-Data\\.csv$",
      full.names = TRUE
    )
    
    sort(
      unique(
        get_year5(files)
      )
    )
  }
)

available_years5 <- Reduce(
  intersect,
  years_by_folder5
)

available_years5


# ============================================================
# 11. Extract one ACS 5-year endpoint year
# ============================================================

extract_acs5_year <- function(year) {
  
  message(
    "Processing ACS 5-year: ",
    year
  )
  
  
  # ----------------------------------------------------------
  # AGE
  # S0101
  #
  # Table structure changed in 2017
  # ----------------------------------------------------------
  
  f <- find_data_file5(
    "ACS_age&sex_5-year",
    year
  )
  
  age_raw <- read_csv(
    f,
    show_col_types = FALSE
  ) %>%
    slice(-1)
  
  if (year <= 2016) {
    
    age <- age_raw %>%
      transmute(
        FIPS = get_fips(GEO_ID),
        County = NAME,
        
        median_age =
          to_numeric(
            S0101_C01_030E
          ),
        
        pct_age65 =
          to_numeric(
            S0101_C01_028E
          )
      )
    
  } else {
    
    age <- age_raw %>%
      transmute(
        FIPS = get_fips(GEO_ID),
        County = NAME,
        
        median_age =
          to_numeric(
            S0101_C01_032E
          ),
        
        pct_age65 =
          to_numeric(
            S0101_C02_030E
          )
      )
  }
  
  
  # ----------------------------------------------------------
  # MEDIAN HOUSEHOLD INCOME
  # S1901
  #
  # Explicit table selection is required because this folder
  # may also contain S1903 files.
  # ----------------------------------------------------------
  
  f <- find_data_file5_table(
    "ACS_median_income_5-year",
    year,
    "S1901"
  )
  
  income <- read_csv(
    f,
    show_col_types = FALSE
  ) %>%
    slice(-1) %>%
    transmute(
      FIPS = get_fips(GEO_ID),
      
      median_household_income =
        to_numeric(
          S1901_C01_012E
        )
    )
  
  
  # ----------------------------------------------------------
  # EDUCATION
  # S1501
  #
  # Table structure changed in 2015
  # ----------------------------------------------------------
  
  f <- find_data_file5(
    "ACS_education_5-year",
    year
  )
  
  education_raw <- read_csv(
    f,
    show_col_types = FALSE
  ) %>%
    slice(-1)
  
  if (year <= 2014) {
    
    education <- education_raw %>%
      transmute(
        FIPS = get_fips(GEO_ID),
        
        pct_bachelors =
          to_numeric(
            S1501_C01_015E
          )
      )
    
  } else {
    
    education <- education_raw %>%
      transmute(
        FIPS = get_fips(GEO_ID),
        
        pct_bachelors =
          to_numeric(
            S1501_C02_015E
          )
      )
  }
  
  
  # ----------------------------------------------------------
  # MERGE FIVE-YEAR COVARIATES
  # ----------------------------------------------------------
  
  out <- age %>%
    
    left_join(
      income,
      by = "FIPS"
    ) %>%
    
    left_join(
      education,
      by = "FIPS"
    ) %>%
    
    mutate(
      
      Year = year,
      
      ACS_start_year =
        year - 4,
      
      ACS_end_year =
        year
    ) %>%
    
    select(
      Year,
      ACS_start_year,
      ACS_end_year,
      FIPS,
      County,
      median_household_income,
      median_age,
      pct_age65,
      pct_bachelors
    )
  
  return(out)
}


# ============================================================
# 12. Build complete ACS 5-year panel
# ============================================================

acs_covariates_5year <- map_dfr(
  available_years5,
  extract_acs5_year
) %>%
  arrange(
    FIPS,
    Year
  )


# ============================================================
# 13. Quality checks: ACS 5-year
# ============================================================

dim(
  acs_covariates_5year
)

head(
  acs_covariates_5year
)


# County counts by endpoint year

acs_covariates_5year %>%
  count(Year)


# Missing values

acs_covariates_5year %>%
  group_by(Year) %>%
  summarise(
    
    n = n(),
    
    income_missing =
      sum(
        is.na(
          median_household_income
        )
      ),
    
    age_missing =
      sum(
        is.na(
          median_age
        )
      ),
    
    age65_missing =
      sum(
        is.na(
          pct_age65
        )
      ),
    
    education_missing =
      sum(
        is.na(
          pct_bachelors
        )
      )
  )


# Plausible ranges

summary(
  acs_covariates_5year %>%
    select(
      median_household_income,
      median_age,
      pct_age65,
      pct_bachelors
    )
)


# Impossible percentage values

acs_covariates_5year %>%
  filter(
    pct_age65 < 0 |
      pct_age65 > 100 |
      pct_bachelors < 0 |
      pct_bachelors > 100
  )


# Specifically inspect missing household-income observations

acs_covariates_5year %>%
  filter(
    is.na(
      median_household_income
    )
  ) %>%
  select(
    Year,
    FIPS,
    County
  )


# ============================================================
# 14. Save ACS 5-year panel
# ============================================================

write_csv(
  acs_covariates_5year,
  "ACS_county_covariates_5year.csv"
)

saveRDS(
  acs_covariates_5year,
  "ACS_county_covariates_5year.rds"
)


# ============================================================
# End
# ============================================================

message(
  "ACS covariate construction complete."
)

message(
  "Annual panel: ",
  nrow(acs_covariates_annual),
  " county-year observations."
)

message(
  "Five-year panel: ",
  nrow(acs_covariates_5year),
  " county-year observations."
)