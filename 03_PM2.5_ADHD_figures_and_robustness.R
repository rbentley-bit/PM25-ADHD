# ============================================================
# 03_PM2.5_ADHD_robustness.R
#
# Statistical analyses and manuscript figures for the revised manuscript.
#
# Assumes these scripts have already been run:
#   01_build_ACS_covariates.R
#   02_build_health_PM25_analysis_data.R
#
# Main analyses:
#   1. Low birth weight
#      - baseline county/year fixed effects
#      - county-specific linear temporal trends
#      - contemporaneous annual ACS covariates
#      - lagged models with and without county trends
#   2. Diabetes
#      - baseline county/year fixed effects
#      - county-specific linear temporal trends
#      - contemporaneous annual ACS covariates
#      - lagged models with and without county trends
#   3. Childhood ADHD
#      - common-sample unadjusted model
#      - ACS 5-year + state fixed-effects model
#      - alternative PM2.5 exposure windows
#      - outlier sensitivity
#      - Moran's I
#      - spatial-error model
#   4. Manuscript figures
#      - Figure 2 descriptive county maps
#      - Figure 3 ADHD partial regression and exposure windows
#      - Figure 4 covariate-adjusted residual maps
#
# This script does NOT reconstruct raw health or ACS datasets.
# ============================================================

# ============================================================
# 0. Setup
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(fixest)
library(ggplot2)
library(viridis)
library(tigris)
library(sf)
library(spdep)
library(spatialreg)
library(patchwork)
library(grid)

options(tigris_use_cache = TRUE)

# ============================================================
# 1. Load clean outputs from Scripts 01 and 02
# ============================================================

acs_annual <- readRDS("ACS_county_covariates_annual.rds") %>% mutate(FIPS = as.integer(FIPS))
acs_5year <- readRDS("ACS_county_covariates_5year.rds") %>% mutate(FIPS = as.integer(FIPS))
pm25 <- readRDS("PM25_county_year.rds") %>% mutate(FIPS = as.integer(FIPS))
lbw_pm25 <- readRDS("LBW_PM25_county_year.rds") %>% mutate(FIPS = as.integer(FIPS))
diabetes_pm25 <- readRDS("Diabetes_PM25_county_year.rds") %>% mutate(FIPS = as.integer(FIPS))
adhd_pm25 <- readRDS("ADHD_PM25_county.rds") %>% mutate(FIPS = as.integer(FIPS))

# Original county covariates are retained only to reproduce
# the covariate-adjusted residual maps in Figure 4.
county_covariates_original <- read_csv("full_covariates_paper.csv", show_col_types = FALSE) %>%
  mutate(FIPS = as.integer(FIPS), FIPS_chr = sprintf("%05d", FIPS))

# ============================================================
# PART I. LOW BIRTH WEIGHT
# ============================================================

# ============================================================
# 2. Baseline LBW fixed-effects model
# ============================================================

model_lbw_fe <- feols(
  low_bw_pct ~ PM25 | FIPS + Year,
  data = lbw_pm25, weights = ~ total_births, cluster = ~ FIPS
)
summary(model_lbw_fe)

# ============================================================
# 3. LBW model with county-specific linear temporal trends
# ============================================================

lbw_trend_data <- lbw_pm25 %>% mutate(Year_c = Year - min(Year, na.rm = TRUE))

model_lbw_fe_trend <- feols(
  low_bw_pct ~ PM25 | FIPS + Year + FIPS[Year_c],
  data = lbw_trend_data, weights = ~ total_births, cluster = ~ FIPS
)
summary(model_lbw_fe_trend)

coef(model_lbw_fe)["PM25"]
coef(model_lbw_fe_trend)["PM25"]

# ============================================================
# 4. LBW + contemporaneous annual ACS covariates
# ============================================================

lbw_acs <- lbw_pm25 %>%
  left_join(
    acs_annual %>% select(
      FIPS, Year, median_household_income, median_age, pct_bachelors,
      pct_poverty, unemployment_rate, pct_uninsured
    ),
    by = c("FIPS", "Year")
  )

lbw_acs_complete <- lbw_acs %>%
  filter(Year >= 2010, Year <= 2019) %>%
  drop_na(
    low_bw_pct, PM25, total_births, median_household_income, median_age,
    pct_bachelors, pct_poverty, unemployment_rate, pct_uninsured
  ) %>%
  mutate(median_household_income_10K = median_household_income / 10000)

nrow(lbw_acs_complete)  # Expected: 5523

model_lbw_fe_same_sample <- feols(
  low_bw_pct ~ PM25 | FIPS + Year,
  data = lbw_acs_complete, weights = ~ total_births, cluster = ~ FIPS
)

model_lbw_fe_annual <- feols(
  low_bw_pct ~ PM25 + median_household_income_10K + pct_bachelors +
    median_age + pct_poverty | FIPS + Year,
  data = lbw_acs_complete, weights = ~ total_births, cluster = ~ FIPS
)

model_lbw_fe_annual_full <- feols(
  low_bw_pct ~ PM25 + median_household_income_10K + pct_bachelors +
    median_age + pct_poverty + unemployment_rate + pct_uninsured | FIPS + Year,
  data = lbw_acs_complete, weights = ~ total_births, cluster = ~ FIPS
)


# ============================================================
# 5. Lagged LBW fixed-effects models
# ============================================================

run_lbw_lag_model <- function(lag_x, data, county_trends = FALSE) {
  pm_lagged <- data %>%
    select(FIPS, Year, PM25) %>%
    mutate(Year = Year + lag_x) %>%
    rename(pm25_lagged = PM25)

  d <- data %>%
    select(FIPS, Year, low_bw_pct, total_births) %>%
    left_join(pm_lagged, by = c("FIPS", "Year")) %>%
    drop_na(low_bw_pct, pm25_lagged, total_births)

  if (nrow(d) < 50 || n_distinct(d$Year) < 2) {
    return(tibble(lag = lag_x, n = nrow(d), beta = NA_real_, se = NA_real_, p_value = NA_real_))
  }

  d <- d %>% mutate(Year_c = Year - min(Year))
  if (county_trends) {
    fit <- tryCatch(
      feols(
        low_bw_pct ~ pm25_lagged | FIPS + Year + FIPS[Year_c],
        data = d, weights = ~ total_births, cluster = ~ FIPS
      ),
      error = function(e) NULL
    )
  } else {
    fit <- feols(
      low_bw_pct ~ pm25_lagged | FIPS + Year,
      data = d, weights = ~ total_births, cluster = ~ FIPS
    )
  }

  if (is.null(fit)) {
    return(tibble(lag = lag_x, n = nrow(d), beta = NA_real_, se = NA_real_, p_value = NA_real_))
  }

  out <- summary(fit)$coeftable["pm25_lagged", ]
  tibble(
    lag = lag_x, n = nobs(fit),
    beta = unname(out["Estimate"]),
    se = unname(out["Std. Error"]),
    p_value = unname(out["Pr(>|t|)"])
  )
}

lbw_lag_results <- map_dfr(0:10, ~ run_lbw_lag_model(.x, lbw_pm25, county_trends = FALSE)) %>%
  mutate(ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se)

lbw_lag_results_trend <- map_dfr(0:10, ~ run_lbw_lag_model(.x, lbw_pm25, county_trends = TRUE)) %>%
  mutate(ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se)

lbw_lag_results
lbw_lag_results_trend

# ============================================================
# PART II. DIABETES
# ============================================================

# ============================================================
# 6. Baseline diabetes fixed-effects model
# ============================================================

model_diabetes_fe <- feols(
  Diabetes ~ PM25 | FIPS + Year,
  data = diabetes_pm25, cluster = ~ FIPS
)
summary(model_diabetes_fe)

# ============================================================
# 7. Diabetes with county-specific linear trends
# ============================================================

diabetes_trend_data <- diabetes_pm25 %>% mutate(Year_c = Year - min(Year, na.rm = TRUE))

model_diabetes_fe_trend <- feols(
  Diabetes ~ PM25 | FIPS + Year + FIPS[Year_c],
  data = diabetes_trend_data, cluster = ~ FIPS
)
summary(model_diabetes_fe_trend)

# ============================================================
# 8. Diabetes + contemporaneous annual ACS covariates
# ============================================================

diabetes_acs <- diabetes_pm25 %>%
  left_join(
    acs_annual %>% select(
      FIPS, Year, median_household_income, median_age, pct_bachelors,
      pct_poverty, unemployment_rate, pct_uninsured
    ),
    by = c("FIPS", "Year")
  )

diabetes_acs_complete <- diabetes_acs %>%
  filter(Year >= 2010, Year <= 2019) %>%
  drop_na(
    Diabetes, PM25, median_household_income, median_age, pct_bachelors,
    pct_poverty, unemployment_rate, pct_uninsured
  ) %>%
  mutate(median_household_income_10K = median_household_income / 10000)

nrow(diabetes_acs_complete)  # Expected: 8015

model_diab_fe_same_sample <- feols(
  Diabetes ~ PM25 | FIPS + Year,
  data = diabetes_acs_complete, cluster = ~ FIPS
)

model_diab_fe_annual <- feols(
  Diabetes ~ PM25 + median_household_income_10K + pct_bachelors +
    median_age + pct_poverty | FIPS + Year,
  data = diabetes_acs_complete, cluster = ~ FIPS
)

model_diab_fe_annual_full <- feols(
  Diabetes ~ PM25 + median_household_income_10K + pct_bachelors +
    median_age + pct_poverty + unemployment_rate + pct_uninsured | FIPS + Year,
  data = diabetes_acs_complete, cluster = ~ FIPS
)


# ============================================================
# 9. Lagged diabetes fixed-effects models
# ============================================================

run_diabetes_lag_model <- function(lag_x, data, county_trends = FALSE) {
  pm_lagged <- data %>%
    select(FIPS, Year, PM25) %>%
    mutate(Year = Year + lag_x) %>%
    rename(pm25_lagged = PM25)

  d <- data %>%
    select(FIPS, Year, Diabetes) %>%
    left_join(pm_lagged, by = c("FIPS", "Year")) %>%
    drop_na(Diabetes, pm25_lagged)

  if (nrow(d) < 50 || n_distinct(d$Year) < 2) {
    return(tibble(
      lag = lag_x, n = nrow(d), n_years = n_distinct(d$Year),
      beta = NA_real_, se = NA_real_, p_value = NA_real_
    ))
  }

  d <- d %>% mutate(Year_c = Year - min(Year))
  if (county_trends) {
    fit <- tryCatch(
      feols(
        Diabetes ~ pm25_lagged | FIPS + Year + FIPS[Year_c],
        data = d, cluster = ~ FIPS
      ),
      error = function(e) NULL
    )
  } else {
    fit <- feols(
      Diabetes ~ pm25_lagged | FIPS + Year,
      data = d, cluster = ~ FIPS
    )
  }

  if (is.null(fit)) {
    return(tibble(
      lag = lag_x, n = nrow(d), n_years = n_distinct(d$Year),
      beta = NA_real_, se = NA_real_, p_value = NA_real_
    ))
  }

  out <- summary(fit)$coeftable["pm25_lagged", ]
  tibble(
    lag = lag_x, n = nobs(fit), n_years = n_distinct(d$Year),
    beta = unname(out["Estimate"]),
    se = unname(out["Std. Error"]),
    p_value = unname(out["Pr(>|t|)"])
  )
}

diabetes_lag_results <- map_dfr(0:15, ~ run_diabetes_lag_model(.x, diabetes_pm25, county_trends = FALSE)) %>%
  mutate(ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se)

diabetes_lag_results_trend <- map_dfr(0:15, ~ run_diabetes_lag_model(.x, diabetes_pm25, county_trends = TRUE)) %>%
  mutate(ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se)

diabetes_lag_results
diabetes_lag_results_trend

# ============================================================
# PART III. CHILDHOOD ADHD
# ============================================================

# ============================================================
# 10. Extract 2018 ACS 5-year covariates
# ============================================================

# 2018 ACS 5-year estimate represents 2014–2018.
acs5_2018 <- acs_5year %>%
  filter(Year == 2018) %>%
  select(FIPS, median_household_income, median_age, pct_age65, pct_bachelors) %>%
  mutate(median_household_income_10K = median_household_income / 10000)

# ============================================================
# 11. Merge ACS 5-year covariates with canonical ADHD dataset
# ============================================================

adhd_acs5 <- adhd_pm25 %>% left_join(acs5_2018, by = "FIPS")

required_adhd_variables <- c(
  "ADHD", "pm25_2005_2015", "pm25_2005_2006", "pm25_2007_2009",
  "pm25_2010_2012", "pm25_2013_2015", "median_household_income_10K",
  "median_age", "pct_bachelors", "State"
)

missing_adhd_variables <- setdiff(required_adhd_variables, names(adhd_acs5))
if (length(missing_adhd_variables) > 0) {
  stop("Missing required ADHD variables: ", paste(missing_adhd_variables, collapse = ", "))
}

# ============================================================
# 12. Create one common ADHD analysis sample
# ============================================================

adhd_common <- adhd_acs5 %>%
  drop_na(
    ADHD, pm25_2005_2015, pm25_2005_2006, pm25_2007_2009,
    pm25_2010_2012, pm25_2013_2015, median_household_income_10K,
    median_age, pct_bachelors, State
  )

nrow(adhd_common)  # Expected: 3106

adhd_common %>% count(FIPS) %>% filter(n > 1)

# ============================================================
# 13. Primary ADHD models on identical sample
# ============================================================

model_adhd_same_sample <- lm(ADHD ~ pm25_2005_2015, data = adhd_common)

model_adhd_main <- lm(
  ADHD ~ pm25_2005_2015 + median_household_income_10K +
    median_age + pct_bachelors + factor(State),
  data = adhd_common
)

summary(model_adhd_same_sample)
summary(model_adhd_main)


# ============================================================
# 14. Alternative PM2.5 exposure windows
# ============================================================

run_adhd_window_model <- function(pm_var, data) {
  formula_window <- as.formula(
    paste0(
      "ADHD ~ ", pm_var,
      " + median_household_income_10K + median_age + pct_bachelors + factor(State)"
    )
  )

  fit <- lm(formula_window, data = data)
  out <- summary(fit)$coefficients[pm_var, ]

  tibble(
    exposure_window = pm_var, n = nobs(fit),
    beta = unname(out["Estimate"]),
    se = unname(out["Std. Error"]),
    p_value = unname(out["Pr(>|t|)"]),
    r_squared = summary(fit)$r.squared
  )
}

adhd_window_variables <- c(
  "pm25_2005_2006", "pm25_2007_2009", "pm25_2010_2012", "pm25_2013_2015"
)

adhd_window_results <- map_dfr(
  adhd_window_variables,
  ~ run_adhd_window_model(.x, adhd_common)
) %>%
  mutate(
    ci_low = beta - 1.96 * se,
    ci_high = beta + 1.96 * se,
    exposure_window = factor(
      exposure_window,
      levels = c("pm25_2005_2006", "pm25_2007_2009", "pm25_2010_2012", "pm25_2013_2015"),
      labels = c("2005–2006", "2007–2009", "2010–2012", "2013–2015")
    )
  )

adhd_window_results
write_csv(adhd_window_results, "ADHD_exposure_window_results.csv")
saveRDS(adhd_window_results, "ADHD_exposure_window_results.rds")

# ============================================================
# 15. ADHD sensitivity: PM2.5 distribution outliers
# ============================================================

pm25_limits <- quantile(adhd_common$pm25_2005_2015, probs = c(0.01, 0.99), na.rm = TRUE)

adhd_trimmed <- adhd_common %>%
  filter(pm25_2005_2015 >= pm25_limits[1], pm25_2005_2015 <= pm25_limits[2])

model_adhd_trimmed <- lm(
  ADHD ~ pm25_2005_2015 + median_household_income_10K +
    median_age + pct_bachelors + factor(State),
  data = adhd_trimmed
)

coef(summary(model_adhd_trimmed))["pm25_2005_2015", ]


# ============================================================
# PART IV. ADHD SPATIAL ROBUSTNESS
# ============================================================

# ============================================================
# 17. County and state geometry
# ============================================================

counties_sf <- counties(cb = TRUE, year = 2020, class = "sf") %>%
  mutate(FIPS = as.integer(GEOID), FIPS_chr = GEOID)

states_sf <- states(cb = TRUE, year = 2020, class = "sf") %>%
  filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78"))

counties_sf_contig <- counties_sf %>%
  filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78"))

# ============================================================
# 18. Exact spatial ADHD sample
# ============================================================

adhd_sf <- counties_sf_contig %>%
  inner_join(
    adhd_common %>% select(
      FIPS, ADHD, pm25_2005_2015, median_household_income_10K,
      median_age, pct_bachelors, State
    ),
    by = "FIPS"
  )

nrow(adhd_sf)

# ============================================================
# 19. Queen-contiguity weights
# ============================================================

adhd_nb <- poly2nb(adhd_sf, queen = TRUE)
adhd_lw <- nb2listw(adhd_nb, style = "W", zero.policy = TRUE)
which(card(adhd_nb) == 0)

# ============================================================
# 20. Ordinary adjusted model and Moran's I
# ============================================================

model_adhd_ols_spatial_sample <- lm(
  ADHD ~ pm25_2005_2015 + median_household_income_10K +
    median_age + pct_bachelors + factor(State),
  data = adhd_sf
)

adhd_sf$ols_resid <- residuals(model_adhd_ols_spatial_sample)

moran_adhd_ols <- moran.test(adhd_sf$ols_resid, adhd_lw, zero.policy = TRUE)
moran_adhd_ols

# ============================================================
# 21. Spatial-error model
# ============================================================

model_adhd_spatial_error <- errorsarlm(
  ADHD ~ pm25_2005_2015 + median_household_income_10K +
    median_age + pct_bachelors + factor(State),
  data = adhd_sf, listw = adhd_lw, zero.policy = TRUE
)

summary(model_adhd_spatial_error)
coef(summary(model_adhd_spatial_error))["pm25_2005_2015", ]
model_adhd_spatial_error$lambda

# ============================================================
# 22. Moran's I after spatial-error adjustment
# ============================================================

adhd_sf$spatial_error_resid <- residuals(model_adhd_spatial_error)
moran_adhd_spatial <- moran.test(adhd_sf$spatial_error_resid, adhd_lw, zero.policy = TRUE)
moran_adhd_spatial

# ============================================================
# PART V. SAVE KEY NUMERIC RESULTS
# ============================================================

# ============================================================
# 23. Coefficient-extraction helpers
# ============================================================

extract_lm_coef <- function(fit, term, model_name) {
  tab <- coef(summary(fit))
  tibble(
    model = model_name, term = term,
    beta = unname(tab[term, "Estimate"]),
    se = unname(tab[term, "Std. Error"]),
    p_value = unname(tab[term, grep("^Pr", colnames(tab))]),
    n = nobs(fit)
  )
}

extract_fe_coef <- function(fit, term, model_name, outcome) {
  tab <- summary(fit)$coeftable
  tibble(
    outcome = outcome, model = model_name,
    beta = unname(tab[term, "Estimate"]),
    se = unname(tab[term, "Std. Error"]),
    p_value = unname(tab[term, grep("^Pr", colnames(tab))]),
    n = nobs(fit)
  )
}

# ============================================================
# 24. Table 1: LBW and diabetes
# ============================================================

table1_results <- bind_rows(
  extract_fe_coef(model_lbw_fe_same_sample, "PM25", "County + year FE", "Low birth weight"),
  extract_fe_coef(model_lbw_fe_annual, "PM25", "+ annual ACS", "Low birth weight"),
  extract_fe_coef(model_lbw_fe_annual_full, "PM25", "+ full annual ACS", "Low birth weight"),
  extract_fe_coef(model_diab_fe_same_sample, "PM25", "County + year FE", "Diabetes prevalence"),
  extract_fe_coef(model_diab_fe_annual, "PM25", "+ annual ACS", "Diabetes prevalence"),
  extract_fe_coef(model_diab_fe_annual_full, "PM25", "+ full annual ACS", "Diabetes prevalence")
) %>%
  mutate(ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se)

table1_results
write_csv(table1_results, "Table1_LBW_diabetes_results.csv")

# ============================================================
# 25. Table 2: ADHD robustness
# ============================================================

spatial_tab <- coef(summary(model_adhd_spatial_error))

table2_results <- bind_rows(
  extract_lm_coef(model_adhd_same_sample, "pm25_2005_2015", "Unadjusted, same sample"),
  extract_lm_coef(model_adhd_main, "pm25_2005_2015", "ACS 5-year + State FE"),
  extract_lm_coef(model_adhd_trimmed, "pm25_2005_2015", "Trimmed 1% tails"),
  tibble(
    model = "Spatial error model", term = "pm25_2005_2015",
    beta = unname(spatial_tab["pm25_2005_2015", "Estimate"]),
    se = unname(spatial_tab["pm25_2005_2015", "Std. Error"]),
    p_value = unname(spatial_tab["pm25_2005_2015", grep("^Pr", colnames(spatial_tab))]),
    n = nobs(model_adhd_spatial_error)
  )
) %>%
  mutate(ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se)

table2_results
write_csv(table2_results, "Table2_ADHD_robustness_results.csv")

# ============================================================
# 26. Table 3: Alternative PM2.5 exposure windows
# ============================================================

table3_results <- adhd_window_results %>%
  select(exposure_window, beta, se, ci_low, ci_high, p_value, n)

table3_results
write_csv(table3_results, "Table3_ADHD_exposure_windows.csv")

# ============================================================
# 27. Save key model objects and lag-series results
# ============================================================

saveRDS(model_lbw_fe, "model_lbw_fe.rds")
saveRDS(model_lbw_fe_trend, "model_lbw_fe_trend.rds")
saveRDS(model_lbw_fe_same_sample, "model_lbw_fe_same_sample.rds")
saveRDS(model_lbw_fe_annual, "model_lbw_fe_annual.rds")
saveRDS(model_lbw_fe_annual_full, "model_lbw_fe_annual_full.rds")
saveRDS(model_diabetes_fe, "model_diabetes_fe.rds")
saveRDS(model_diabetes_fe_trend, "model_diabetes_fe_trend.rds")
saveRDS(model_adhd_same_sample, "model_adhd_same_sample.rds")
saveRDS(model_adhd_main, "model_adhd_main.rds")
saveRDS(model_adhd_spatial_error, "model_adhd_spatial_error.rds")

write_csv(lbw_lag_results, "LBW_lag_FE_results.csv")
write_csv(lbw_lag_results_trend, "LBW_lag_FE_county_trend_results.csv")
write_csv(diabetes_lag_results, "Diabetes_lag_FE_results.csv")
write_csv(diabetes_lag_results_trend, "Diabetes_lag_FE_county_trend_results.csv")

# ============================================================
# PART VI. MANUSCRIPT FIGURES
# ============================================================

dir.create("figures", showWarnings = FALSE)

small_legend <- theme(
  legend.position = c(0.87, 0.10),
  legend.justification = c(0, 0),
  legend.direction = "vertical",
  legend.key.height = unit(0.6, "cm"),
  legend.key.width = unit(0.25, "cm"),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 6),
  legend.background = element_rect(fill = "white", color = NA)
)

map_coord <- coord_sf(
  xlim = c(-125, -66), ylim = c(24, 50),
  expand = FALSE, datum = NA
)

# ============================================================
# 28. Figure 2: Descriptive county maps
# ============================================================

pm_map_data <- counties_sf_contig %>%
  left_join(adhd_pm25 %>% select(FIPS, pm25_2005_2015), by = "FIPS")

lbw_map_data <- counties_sf_contig %>%
  left_join(
    lbw_pm25 %>% filter(Year == 2016) %>% select(FIPS, low_bw_rate),
    by = "FIPS"
  )

diab_map_data <- counties_sf_contig %>%
  left_join(
    diabetes_pm25 %>% filter(Year == 2020) %>% select(FIPS, Diabetes),
    by = "FIPS"
  )

adhd_map_data <- counties_sf_contig %>%
  left_join(
    adhd_pm25 %>% select(FIPS, ADHD) %>% rename(ADHD_prev = ADHD),
    by = "FIPS"
  )

pm_map <- ggplot(pm_map_data) +
  geom_sf(aes(fill = pm25_2005_2015), color = NA) +
  geom_sf(data = states_sf, fill = NA, color = "gray85", linewidth = 0.15) +
  map_coord +
  scale_fill_viridis_c(
    option = "magma", na.value = "gray90",
    guide = guide_colorbar(barheight = unit(2, "cm"))
  ) +
  theme_void() + small_legend +
  labs(fill = "(µg/m³)", title = "(a) PM2.5 (µg/m³)")

lbw_map <- ggplot(lbw_map_data) +
  geom_sf(aes(fill = low_bw_rate * 100), color = NA) +
  geom_sf(data = states_sf, fill = NA, color = "gray85", linewidth = 0.15) +
  map_coord +
  scale_fill_viridis_c(
    option = "magma", na.value = "gray90",
    guide = guide_colorbar(barheight = unit(2, "cm"))
  ) +
  theme_void() + small_legend +
  labs(fill = "%", title = "(b) Low birth weight (%)")

diabetes_map <- ggplot(diab_map_data) +
  geom_sf(aes(fill = Diabetes), color = NA) +
  geom_sf(data = states_sf, fill = NA, color = "gray85", linewidth = 0.15) +
  map_coord +
  scale_fill_viridis_c(
    option = "magma", na.value = "gray90",
    guide = guide_colorbar(barheight = unit(2, "cm"))
  ) +
  theme_void() + small_legend +
  labs(fill = "%", title = "(c) Diabetes prevalence (%)")

adhd_map <- ggplot(adhd_map_data) +
  geom_sf(aes(fill = ADHD_prev), color = NA) +
  geom_sf(data = states_sf, fill = NA, color = "gray85", linewidth = 0.15) +
  map_coord +
  scale_fill_viridis_c(
    option = "magma", na.value = "gray90",
    guide = guide_colorbar(barheight = unit(2, "cm"))
  ) +
  theme_void() + small_legend +
  labs(fill = "%", title = "(d) ADHD prevalence (%)")

combined_map <- (pm_map + lbw_map) / (diabetes_map + adhd_map) +
  plot_layout(widths = c(1, 1.1))

combined_map

# ============================================================
# 29. Figure 3a: Partial regression for ADHD
# ============================================================

adhd_resid_new <- resid(
  lm(
    ADHD ~ median_household_income_10K + median_age +
      pct_bachelors + factor(State),
    data = adhd_common
  )
)

pm_resid_new <- resid(
  lm(
    pm25_2005_2015 ~ median_household_income_10K + median_age +
      pct_bachelors + factor(State),
    data = adhd_common
  )
)

partial_plot_data <- adhd_common %>%
  transmute(FIPS, PM_resid = pm_resid_new, ADHD_resid = adhd_resid_new)

binscatter_data_adhd <- partial_plot_data %>%
  mutate(bin = ntile(PM_resid, 20)) %>%
  group_by(bin) %>%
  summarise(
    PM_resid = mean(PM_resid),
    ADHD_resid = mean(ADHD_resid),
    .groups = "drop"
  )

partial_reg_new <- ggplot() +
  geom_point(
    data = partial_plot_data,
    aes(x = PM_resid, y = ADHD_resid),
    alpha = 0.20, size = 1
  ) +
  geom_point(
    data = binscatter_data_adhd,
    aes(x = PM_resid, y = ADHD_resid),
    size = 3, color = "#2171B5"
  ) +
  geom_smooth(
    data = partial_plot_data,
    aes(x = PM_resid, y = ADHD_resid),
    method = "lm", se = TRUE, color = "black",
    fill = "gray80", linewidth = 0.9
  ) +
  coord_cartesian(ylim = c(-5, 5), xlim = c(-3, 3)) +
  theme_classic() +
  theme(
    aspect.ratio = 0.9,
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  labs(
    x = "Residual PM2.5 (µg/m³)",
    y = "Residual ADHD prevalence (%)",
    title = "(a)"
  )

partial_reg_new


# Verify that the partial-regression slope matches the adjusted ADHD model.
coef(lm(ADHD_resid ~ PM_resid, data = partial_plot_data))["PM_resid"]

# ============================================================
# 30. Figure 3b: Alternative PM2.5 exposure windows
# ============================================================

adhd_window_plot_data <- adhd_window_results %>% mutate(window_num = 1:n())

adhd_plot_new <- ggplot(adhd_window_plot_data, aes(x = window_num, y = beta)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "#9ecae1", alpha = 0.35) +
  geom_line(linewidth = 1, color = "#2171B5") +
  geom_point(size = 3, color = "#2171B5") +
  scale_x_continuous(
    breaks = adhd_window_plot_data$window_num,
    labels = as.character(adhd_window_plot_data$exposure_window)
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 0.9,
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  labs(
    x = "PM2.5 exposure window",
    y = "Change in ADHD prevalence (% pts per µg/m³)",
    title = "(b)"
  )

adhd_plot_new

# ============================================================
# 31. Figure 4: Covariate-adjusted residual maps
# ============================================================

# The original county covariates are retained here to reproduce the
# cross-sectional PM2.5, LBW, and diabetes residual maps used in the manuscript.
# The ADHD residual map uses the final 3,106-county ACS-adjusted specification.

map_covariates <- county_covariates_original %>%
  transmute(
    FIPS,
    State_map = State,
    income_map = income,
    mean_age_map = mean_age,
    bachelors_map = bachelors_5y_25over_2022
  )

# PM2.5 residuals
pm_model_data <- adhd_pm25 %>%
  select(FIPS, pm25_2005_2015) %>%
  left_join(map_covariates, by = "FIPS") %>%
  drop_na()

pm_model <- lm(
  pm25_2005_2015 ~ income_map + mean_age_map + bachelors_map + factor(State_map),
  data = pm_model_data
)
pm_model_data <- pm_model_data %>% mutate(pm_resid = residuals(pm_model))

# LBW residuals, 2020
lbw_model_data <- lbw_pm25 %>%
  filter(Year == 2020) %>%
  select(FIPS, low_bw_rate) %>%
  left_join(map_covariates, by = "FIPS") %>%
  drop_na()

lbw_model <- lm(
  low_bw_rate ~ income_map + mean_age_map + bachelors_map + factor(State_map),
  data = lbw_model_data
)
lbw_model_data <- lbw_model_data %>% mutate(lbw_resid = 100 * residuals(lbw_model))

# Diabetes residuals, 2020
diab_model_data <- diabetes_pm25 %>%
  filter(Year == 2020) %>%
  select(FIPS, Diabetes) %>%
  left_join(map_covariates, by = "FIPS") %>%
  drop_na()

diab_model <- lm(
  Diabetes ~ income_map + mean_age_map + bachelors_map + factor(State_map),
  data = diab_model_data
)
diab_model_data <- diab_model_data %>% mutate(diab_resid = residuals(diab_model))

# Map data
pm_resid_map_data <- counties_sf_contig %>%
  left_join(pm_model_data %>% select(FIPS, pm_resid), by = "FIPS")

lbw_resid_map_data <- counties_sf_contig %>%
  left_join(lbw_model_data %>% select(FIPS, lbw_resid), by = "FIPS")

diab_resid_map_data <- counties_sf_contig %>%
  left_join(diab_model_data %>% select(FIPS, diab_resid), by = "FIPS")

adhd_resid_map_data <- counties_sf_contig %>%
  left_join(
    st_drop_geometry(adhd_sf) %>% select(FIPS, adhd_resid_pp = ols_resid),
    by = "FIPS"
  )

# Moran's I for the residual maps reported in the manuscript
pm_resid_sf <- counties_sf_contig %>%
  inner_join(pm_model_data %>% select(FIPS, pm_resid), by = "FIPS")
pm_resid_nb <- poly2nb(pm_resid_sf)
pm_resid_lw <- nb2listw(pm_resid_nb, style = "W", zero.policy = TRUE)
moran_pm_resid <- moran.test(pm_resid_sf$pm_resid, pm_resid_lw, zero.policy = TRUE)
moran_pm_resid

diab_resid_sf <- counties_sf_contig %>%
  inner_join(diab_model_data %>% select(FIPS, diab_resid), by = "FIPS")
diab_resid_nb <- poly2nb(diab_resid_sf)
diab_resid_lw <- nb2listw(diab_resid_nb, style = "W", zero.policy = TRUE)
moran_diab_resid <- moran.test(diab_resid_sf$diab_resid, diab_resid_lw, zero.policy = TRUE)
moran_diab_resid

# ADHD Moran's I is the value already calculated from the final adjusted model.
moran_adhd_ols

pm_resid_map <- ggplot(pm_resid_map_data) +
  geom_sf(aes(fill = pm_resid), color = NA) +
  geom_sf(data = states_sf, fill = NA, color = "gray85", linewidth = 0.15) +
  map_coord +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    na.value = "gray90", guide = guide_colorbar(barheight = unit(2, "cm"))
  ) +
  theme_void() + small_legend +
  labs(fill = "(µg/m³)", title = "(a) Residual PM2.5 exposure")

lbw_resid_map <- ggplot(lbw_resid_map_data) +
  geom_sf(aes(fill = lbw_resid), color = NA) +
  geom_sf(data = states_sf, fill = NA, color = "gray85", linewidth = 0.15) +
  map_coord +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    na.value = "gray90", guide = guide_colorbar(barheight = unit(2, "cm"))
  ) +
  theme_void() + small_legend +
  labs(fill = "(%)", title = "(b) Residual low birth weight")

diab_resid_map <- ggplot(diab_resid_map_data) +
  geom_sf(aes(fill = diab_resid), color = NA) +
  geom_sf(data = states_sf, fill = NA, color = "gray85", linewidth = 0.15) +
  map_coord +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    na.value = "gray90", guide = guide_colorbar(barheight = unit(2, "cm"))
  ) +
  theme_void() + small_legend +
  labs(fill = "(%)", title = "(c) Residual diabetes prevalence")

adhd_resid_map <- ggplot(adhd_resid_map_data) +
  geom_sf(aes(fill = adhd_resid_pp), color = NA) +
  geom_sf(data = states_sf, fill = NA, color = "gray85", linewidth = 0.15) +
  map_coord +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    limits = c(-8, 8), oob = scales::squish,
    na.value = "gray90", guide = guide_colorbar(barheight = unit(2, "cm"))
  ) +
  theme_void() + small_legend +
  labs(fill = "(%)", title = "(d) Residual ADHD prevalence")

four_part_resid_map <- (pm_resid_map + lbw_resid_map) /
  (diab_resid_map + adhd_resid_map) +
  plot_layout(widths = c(1, 1.1))

four_part_resid_map

# ============================================================
# 32. Finish
# ============================================================

message("")
message("Analyses and manuscript figures complete.")
message("LBW annual ACS sample: ", nrow(lbw_acs_complete), " county-year observations.")
message("Diabetes annual ACS sample: ", nrow(diabetes_acs_complete), " county-year observations.")
message("ADHD common sample: ", nrow(adhd_common), " counties.")
message("ADHD spatial sample: ", nrow(adhd_sf), " counties.")
message("")
message("Figure files written to ./figures/")
message("Check ADHD_exposure_window_results.csv for revised alternative-window coefficients.")
