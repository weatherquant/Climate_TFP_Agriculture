## AGRICULTURE TFP - CLIMATE ANALYSIS: REVISION VERSION

library(urca)
library(vars)
library(tseries)
library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)
library(tibble)

## Optional but recommended
library(uroot)       # ADF-GLS / ERS alternatives if available
library(strucchange) # breakpoints / Chow-type tests

setwd("C:/Users/jwest.INTERNAL/OneDrive - Bureau of Meteorology/Documents/Data/Economics/Agriculture/")

tfp    <- read_csv("TFP_indexes.csv")
output <- read_csv("Agricultural_Growth.csv")
inputs <- read_csv("Ag_Indexes.csv")
temp   <- read_csv("Global_Land_Temp.csv")
spei   <- read_csv("Main_Spei.csv")

years <- 1961:2021

tfp    <- tfp    |> filter(Year %in% years)
output <- output |> filter(Year %in% years)
inputs <- inputs |> filter(Year %in% years)
temp   <- temp   |> filter(Year %in% years)
spei   <- spei   |> filter(Year %in% years)

## ------------------------------------------------------------------
## 1. Standardise names
## ------------------------------------------------------------------

colnames(tfp) <- c("Year", "SSA", "CentralEurope", "NorthAmerica", "Australia")

inputs_raw <- read_csv("Ag_Indexes.csv", show_col_types = FALSE)

# Force expected names, matching your original working script
names(inputs_raw) <- c(
  "Country", "Region", "Subregion", "Inc",
  "Year", "Attribute", "Value"
)

inputs <- inputs_raw |>
  filter(Year %in% years) |>
  rename(
    Country_raw   = Country,
    Region_raw    = Region,
    Subregion_raw = Subregion,
    Inc_raw       = Inc,
    variable      = Attribute,
    value         = Value
  ) |>
  mutate(
    Year  = as.integer(Year),
    value = as.numeric(value),
    Country = case_when(
      Country_raw == "United States"   ~ "UnitedStates",
      Country_raw == "Canada"          ~ "Canada",
      Country_raw == "Australia"       ~ "Australia",
      Country_raw == "Europe, Central" ~ "CentralEurope",
      Country_raw == "SSA, Total"      ~ "SSA",
      TRUE                             ~ Country_raw
    )
  )


temp <- temp |>
  mutate(
    Location = case_when(
      Location == "US"          ~ "UnitedStates",
      Location == "Canada"      ~ "Canada",
      Location == "Australia"   ~ "Australia",
      Location == "Central_Eur" ~ "CentralEurope",
      Location == "CentralEurope" ~ "CentralEurope",
      Location == "SSA"         ~ "SSA",
      TRUE                      ~ Location
    )
  )

spei <- spei |>
  mutate(
    Country = case_when(
      Country == "US"          ~ "UnitedStates",
      Country == "Canada"      ~ "Canada",
      Country == "Australia"   ~ "Australia",
      Country == "Central_Eur" ~ "CentralEurope",
      Country == "Sudan"       ~ "Sudan",
      Country == "Nigeria"     ~ "Nigeria",
      Country == "SSA"         ~ "SSA",
      TRUE                     ~ Country
    )
  )

## ------------------------------------------------------------------
## 2. Build climate series and robustness proxies
## ------------------------------------------------------------------

## If no production weights are available, start with equal-weight composites.
## Replace weights later with agricultural output or cropland weights if available.

make_equal_composite <- function(df, id_col, value_col, members, new_name) {
  id_col <- rlang::ensym(id_col)
  value_col <- rlang::ensym(value_col)
  
  df |>
    filter(!!id_col %in% members) |>
    group_by(Year) |>
    summarise(value = mean(!!value_col, na.rm = TRUE), .groups = "drop") |>
    mutate(series = new_name)
}

na_temp_comp <- make_equal_composite(
  temp, Location, Value,
  members = c("UnitedStates", "Canada"),
  new_name = "NorthAmerica"
)

na_spei_comp <- make_equal_composite(
  spei, Country, SPEI,
  members = c("UnitedStates", "Canada"),
  new_name = "NorthAmerica"
)

## SSA robustness proxies:
## 1. Use regional SSA if available.
## 2. Use Sudan proxy.
## 3. Use Nigeria proxy.
## 4. Equal average of available SSA proxies if both Sudan/Nigeria exist.

ssa_spei_multi <- spei |>
  filter(Country %in% c("Sudan", "Nigeria", "SSA")) |>
  select(Year, Country, SPEI)

ssa_spei_composite <- ssa_spei_multi |>
  filter(Country %in% c("Sudan", "Nigeria")) |>
  group_by(Year) |>
  summarise(SPEI = mean(SPEI, na.rm = TRUE), .groups = "drop") |>
  mutate(Country = "SSA_SudanNigeriaComposite")

## Wide temp: replace or add NorthAmerica composite
temp_clean <- temp |>
  select(Year, series = Location, value = Value) |>
  bind_rows(na_temp_comp)

temp.w <- temp_clean |>
  group_by(Year, series) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = series, values_from = value)

## Wide SPEI: include NA composite and SSA alternatives
spei_clean <- spei |>
  select(Year, series = Country, value = SPEI) |>
  bind_rows(na_spei_comp) |>
  bind_rows(
    ssa_spei_composite |>
      transmute(Year, series = Country, value = SPEI)
  )

spei.w <- spei_clean |>
  group_by(Year, series) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = series, values_from = value)

## ------------------------------------------------------------------
## 3. Inputs: construct consistent plant/animal/input variables
## ------------------------------------------------------------------

inputs_base <- inputs |>
  filter(Country %in% c("Australia", "UnitedStates", "Canada", "CentralEurope", "SSA"))

## If Canada input data are available, construct NA composite.
## Otherwise keep NorthAmerica as UnitedStates-only and label robustness appropriately.

inputs_na <- inputs_base |>
  filter(Country %in% c("UnitedStates", "Canada")) |>
  group_by(Year, variable) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(Country = "NorthAmerica")

inputs_regional <- inputs_base |>
  mutate(
    Country = case_when(
      Country == "UnitedStates" ~ "NorthAmerica_USonly",
      TRUE ~ Country
    )
  ) |>
  bind_rows(inputs_na)

## Prefer composite NorthAmerica where possible
inputs_use <- inputs_regional |>
  filter(Country %in% c("Australia", "NorthAmerica", "CentralEurope", "SSA"))

## Construct ratios
inputs_ratios <- inputs_use |>
  pivot_wider(names_from = variable, values_from = value) |>
  mutate(
    CapLab_Index = Capital_Index / Labor_Index,
    CropPast_Index = Cropland_Q / Pasture_Q
  ) |>
  pivot_longer(
    cols = c(CapLab_Index, CropPast_Index),
    names_to = "variable",
    values_to = "value"
  )

inputs_long <- inputs_use |>
  bind_rows(inputs_ratios)

## ------------------------------------------------------------------
## 4. Correct unit-root diagnostics
## ------------------------------------------------------------------

adf_safe <- function(x) {
  x <- na.omit(as.numeric(x))
  out <- tryCatch(adf.test(x), error = function(e) NULL)
  if (is.null(out)) {
    return(tibble(statistic = NA_real_, p_value = NA_real_))
  }
  tibble(
    statistic = as.numeric(out$statistic),
    p_value = out$p.value
  )
}

kpss_safe <- function(x) {
  x <- na.omit(as.numeric(x))
  out <- tryCatch(kpss.test(x), error = function(e) NULL)
  if (is.null(out)) {
    return(tibble(statistic = NA_real_, p_value = NA_real_))
  }
  tibble(
    statistic = as.numeric(out$statistic),
    p_value = out$p.value
  )
}

unit_root_table <- function(series_df, series_name, x) {
  tibble(
    Series = series_name,
    ADF_level_stat = adf_safe(x)$statistic,
    ADF_level_p    = adf_safe(x)$p_value,
    KPSS_level_stat = kpss_safe(x)$statistic,
    KPSS_level_p    = kpss_safe(x)$p_value,
    ADF_diff_stat  = adf_safe(diff(x))$statistic,
    ADF_diff_p     = adf_safe(diff(x))$p_value,
    KPSS_diff_stat = kpss_safe(diff(x))$statistic,
    KPSS_diff_p    = kpss_safe(diff(x))$p_value
  )
}

regions <- c("Australia", "NorthAmerica", "CentralEurope", "SSA")

diagnostics_tfp <- map_dfr(regions, function(r) {
  bind_rows(
    unit_root_table(NULL, paste0("TFP_", r), log(tfp[[r]])),
    unit_root_table(NULL, paste0("TEMP_", r), temp.w[[r]]),
    unit_root_table(NULL, paste0("SPEI_", r), spei.w[[r]])
  )
})

print(diagnostics_tfp)

## ------------------------------------------------------------------
## 5. VECM helpers
## ------------------------------------------------------------------

select_rank_trace <- function(joh) {
  teststat <- joh@teststat
  cval <- joh@cval[, "5pct"]
  
  ## ca.jo reports r <= k tests.
  ## Conservative simple rule: selected rank = number of rejected nulls.
  rank <- sum(teststat > cval)
  rank <- max(1, min(rank, ncol(joh@x) - 1))
  return(rank)
}

estimate_vecm_system <- function(df, max_lag = 8, forced_rank = NULL) {
  
  df_ts <- ts(df, start = 1961, frequency = 1)
  lag_sel <- VARselect(df_ts, lag.max = max_lag, type = "const")
  p <- as.integer(lag_sel$selection["AIC(n)"])
  p <- max(2, p)
  
  joh <- ca.jo(
    df_ts,
    type = "trace",
    ecdet = "const",
    K = p
  )
  
  r_hat <- if (is.null(forced_rank)) select_rank_trace(joh) else forced_rank
  vecm <- cajorls(joh, r = r_hat)
  
  n_eff <- nrow(df) - p
  k <- ncol(df)
  approx_params <- k * (r_hat + k * (p - 1) + 1)
  
  list(
    johansen = joh,
    vecm = vecm,
    lag_AIC = p,
    rank = r_hat,
    n_eff = n_eff,
    k = k,
    approx_params = approx_params,
    param_obs_ratio = approx_params / n_eff
  )
}

extract_response_coeffs <- function(model, response_regex) {
  summ <- summary(model$vecm$rlm)
  eq_name <- names(summ)[grepl(response_regex, names(summ))][1]
  if (is.na(eq_name)) stop("No matching response equation.")
  as.data.frame(summ[[eq_name]]$coefficients) |>
    rownames_to_column("term")
}

extract_ect_only <- function(model, response_regex, region, model_name) {
  extract_response_coeffs(model, response_regex) |>
    filter(grepl("^ect", term)) |>
    transmute(
      Region = region,
      Model = model_name,
      ECT = term,
      Estimate = Estimate,
      StdError = `Std. Error`,
      p_value = `Pr(>|t|)`
    )
}

## ------------------------------------------------------------------
## 6. TFP models
## ------------------------------------------------------------------

estimate_vecm_tfp <- function(region, spei_series = region, temp_series = region) {
  
  df <- data.frame(
    TFP  = log(tfp[[region]]),
    TEMP = temp.w[[temp_series]],
    SPEI = spei.w[[spei_series]]
  )
  
  estimate_vecm_system(df, max_lag = 8)
}

####################################
## SIDETRIP
tfp_specs <- tibble::tribble(
  ~region,          ~temp_series,      ~spei_series,
  "Australia",      "Australia",       "Australia",
  "NorthAmerica",   "NorthAmerica",    "NorthAmerica",
  "CentralEurope",  "CentralEurope",   "CentralEurope",
  "SSA",            "SSA",             "SSA_SudanNigeriaComposite"
)

diagnostics_tfp <- purrr::pmap_dfr(
  list(tfp_specs$region, tfp_specs$temp_series, tfp_specs$spei_series),
  function(region, temp_series, spei_series) {
    dplyr::bind_rows(
      unit_root_table(NULL, paste0("TFP_", region), log(tfp[[region]])),
      unit_root_table(NULL, paste0("TEMP_", region), temp.w[[temp_series]]),
      unit_root_table(NULL, paste0("SPEI_", region, "_", spei_series), spei.w[[spei_series]])
    )
  }
)

print(diagnostics_tfp)

tfp_models <- purrr::pmap(
  list(tfp_specs$region, tfp_specs$spei_series, tfp_specs$temp_series),
  ~ estimate_vecm_tfp(region = ..1, spei_series = ..2, temp_series = ..3)
)

names(tfp_models) <- tfp_specs$region

tfp_ect_table <- map_dfr(regions, function(r) {
  extract_ect_only(tfp_models[[r]], "^Response\\s+TFP\\.d", r, "TFP")
})

tfp_model_diagnostics <- map_dfr(regions, function(r) {
  m <- tfp_models[[r]]
  tibble(
    Region = r,
    lag_AIC = m$lag_AIC,
    rank = m$rank,
    n_eff = m$n_eff,
    k = m$k,
    approx_params = m$approx_params,
    param_obs_ratio = m$param_obs_ratio
  )
})

print(tfp_ect_table)
print(tfp_model_diagnostics)

### PRODUCE TABLE 1
format_ect_cell <- function(est, se, p) {
  stars <- dplyr::case_when(
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE ~ ""
  )
  paste0(sprintf("%.3f", est), stars, " (", sprintf("%.3f", se), ")")
}

tfp_table1_new <- tfp_ect_table |>
  mutate(
    Cell = format_ect_cell(Estimate, StdError, p_value)
  ) |>
  select(Region, ECT, Cell) |>
  tidyr::pivot_wider(names_from = ECT, values_from = Cell) |>
  left_join(
    tfp_model_diagnostics |> select(Region, lag_AIC, rank, n_eff, param_obs_ratio),
    by = "Region"
  )

print(tfp_table1_new)

### DIECT COMPARISON FOR AUSTRALIA
estimate_vecm_tfp_forced <- function(region, spei_series = region, temp_series = region, forced_rank = 2) {
  df <- data.frame(
    TFP  = log(tfp[[region]]),
    TEMP = temp.w[[temp_series]],
    SPEI = spei.w[[spei_series]]
  ) |> tidyr::drop_na()
  
  estimate_vecm_system(df, max_lag = 8, forced_rank = forced_rank)
}

tfp_models_forced_r2 <- purrr::pmap(
  list(tfp_specs$region, tfp_specs$spei_series, tfp_specs$temp_series),
  ~ estimate_vecm_tfp_forced(region = ..1, spei_series = ..2, temp_series = ..3, forced_rank = 2)
)

names(tfp_models_forced_r2) <- tfp_specs$region

tfp_ect_forced_r2 <- purrr::map_dfr(tfp_specs$region, function(r) {
  extract_ect_only(tfp_models_forced_r2[[r]], "^Response\\s+TFP\\.d", r, "TFP_forced_r2")
})

print(tfp_ect_forced_r2)



## ------------------------------------------------------------------
## 7. Component models
## ------------------------------------------------------------------

estimate_vecm_component <- function(varname, region,
                                    temp_series = region,
                                    spei_series = region,
                                    start_year = 1961,
                                    end_year = 2021,
                                    max_lag = 6) {
  
  y_df <- inputs_long |>
    filter(
      Country == region,
      variable == varname,
      Year >= start_year,
      Year <= end_year
    ) |>
    arrange(Year)
  
  clim_df <- tibble(
    Year = years,
    TEMP = temp.w[[temp_series]],
    SPEI = spei.w[[spei_series]]
  ) |>
    filter(Year >= start_year, Year <= end_year)
  
  df <- y_df |>
    select(Year, value) |>
    inner_join(clim_df, by = "Year") |>
    arrange(Year) |>
    transmute(
      Y = log(value),
      TEMP = TEMP,
      SPEI = SPEI
    )
  
  estimate_vecm_system(df, max_lag = max_lag)
}

components <- c(
  CapLab_Index = "Capital-Labour",
  CropPast_Index = "Crop-Pasture",
  Livestock_Q = "Livestock",
  Fertilizer_Q = "Fertiliser"
)

component_models <- list()

# Component model specifications
component_specs <- tibble::tribble(
  ~region,          ~temp_series,      ~spei_series,
  "Australia",      "Australia",       "Australia",
  "NorthAmerica",   "NorthAmerica",    "NorthAmerica",
  "CentralEurope",  "CentralEurope",   "CentralEurope",
  "SSA",            "SSA",             "SSA_SudanNigeriaComposite"
)

# Estimate all component models
component_models <- list()

for (v in names(components)) {
  component_models[[v]] <- purrr::pmap(
    list(
      component_specs$region,
      component_specs$temp_series,
      component_specs$spei_series
    ),
    function(region, temp_series, spei_series) {
      estimate_vecm_component(
        varname     = v,
        region      = region,
        temp_series = temp_series,
        spei_series = spei_series
      )
    }
  )
  
  names(component_models[[v]]) <- component_specs$region
}

# Extract ECT coefficients
component_ect_table <- purrr::map_dfr(names(components), function(v) {
  purrr::map_dfr(component_specs$region, function(r) {
    extract_ect_only(
      component_models[[v]][[r]],
      "^Response\\s+Y\\.d",
      r,
      components[[v]]
    )
  })
})

# Component model diagnostics
component_model_diagnostics <- purrr::map_dfr(names(components), function(v) {
  purrr::map_dfr(component_specs$region, function(r) {
    m <- component_models[[v]][[r]]
    
    tibble::tibble(
      Component = components[[v]],
      Region = r,
      lag_AIC = m$lag_AIC,
      rank = m$rank,
      n_eff = m$n_eff,
      k = m$k,
      approx_params = m$approx_params,
      param_obs_ratio = m$param_obs_ratio
    )
  })
})

print(component_ect_table)
print(component_model_diagnostics)

# Dominant ECT for Table 2
table2_dom <- component_ect_table |>
  dplyr::group_by(Region, Model) |>
  dplyr::slice_max(abs(Estimate), n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Stars = dplyr::case_when(
      p_value < 0.01 ~ "***",
      p_value < 0.05 ~ "**",
      p_value < 0.10 ~ "*",
      TRUE ~ ""
    ),
    Cell = paste0(
      sprintf("%.3f", Estimate),
      Stars,
      " (",
      sprintf("%.3f", StdError),
      ")"
    )
  )

table2_final <- table2_dom |>
  dplyr::select(Component = Model, Region, Cell) |>
  tidyr::pivot_wider(
    names_from = Region,
    values_from = Cell
  )

print(table2_final)

## ADDED ANALYSIS
adjustment_summary <- component_ect_table |>
  mutate(
    convergent = Estimate < 0,
    significant_10 = p_value < 0.10,
    significant_5 = p_value < 0.05,
    convergent_sig_10 = Estimate < 0 & p_value < 0.10,
    convergent_sig_5 = Estimate < 0 & p_value < 0.05
  ) |>
  group_by(Model) |>
  summarise(
    n_terms = n(),
    n_negative = sum(convergent, na.rm = TRUE),
    n_positive = sum(Estimate > 0, na.rm = TRUE),
    n_sig_10 = sum(significant_10, na.rm = TRUE),
    n_negative_sig_10 = sum(convergent_sig_10, na.rm = TRUE),
    n_negative_sig_5 = sum(convergent_sig_5, na.rm = TRUE),
    mean_abs_estimate = mean(abs(Estimate), na.rm = TRUE),
    .groups = "drop"
  )

print(adjustment_summary)

###################################################
## Revised Figures

## ============================================================
## FIGURES 3–6: CORRECTED CONDITIONAL PROJECTIONS, 1961–2050
## ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)
library(vars)

## ------------------------------------------------------------
## 1. Region specifications
## ------------------------------------------------------------

figure_specs <- tibble::tribble(
  ~region,          ~region_label,           ~temp_series,      ~spei_series,
  "NorthAmerica",   "North America",         "NorthAmerica",    "NorthAmerica",
  "CentralEurope",  "Central Europe",        "CentralEurope",   "CentralEurope",
  "SSA",            "Sub-Saharan Africa",    "SSA",             "SSA_SudanNigeriaComposite",
  "Australia",      "Australia",             "Australia",       "Australia"
)

components <- c(
  CapLab_Index  = "Capital-labour ratio",
  CropPast_Index = "Crop-pasture ratio",
  Livestock_Q   = "Livestock",
  Fertilizer_Q  = "Fertiliser"
)

## ------------------------------------------------------------
## 2. Helper: fit component model if not already stored
## ------------------------------------------------------------
## Assumes you already have:
## - inputs_long
## - temp.w
## - spei.w
## - estimate_vecm_component()
## - component_models
##
## If component_models already exists, this will use it.
## Otherwise it will estimate models from scratch.

if (!exists("component_models")) {
  
  component_models <- list()
  
  for (v in names(components)) {
    component_models[[v]] <- purrr::pmap(
      list(
        figure_specs$region,
        figure_specs$temp_series,
        figure_specs$spei_series
      ),
      function(region, temp_series, spei_series) {
        estimate_vecm_component(
          varname     = v,
          region      = region,
          temp_series = temp_series,
          spei_series = spei_series
        )
      }
    )
    
    names(component_models[[v]]) <- figure_specs$region
  }
}

## ------------------------------------------------------------
## 3. Helper: convert VECM to VAR and forecast
## ------------------------------------------------------------

forecast_component <- function(varname,
                               component_label,
                               region,
                               region_label,
                               horizon = 29) {
  
  m <- component_models[[varname]][[region]]
  
  # Convert cajorls/ca.jo object to VAR representation
  var_obj <- vars::vec2var(m$johansen, r = m$rank)
  
  fc <- predict(var_obj, n.ahead = horizon)
  
  # Extract Y forecast
  y_fc <- as.data.frame(fc$fcst$Y)
  
  # Column names are usually fcst, lower, upper, CI
  y_fc <- y_fc |>
    mutate(
      Year = 2022:(2021 + horizon),
      value = fcst,
      lower = lower,
      upper = upper,
      source = "Projection"
    ) |>
    select(Year, value, lower, upper, source)
  
  # Historical observed component
  y_hist <- inputs_long |>
    filter(
      Country == region,
      variable == varname,
      Year >= 1961,
      Year <= 2021
    ) |>
    arrange(Year) |>
    transmute(
      Year,
      value = log(value),
      lower = NA_real_,
      upper = NA_real_,
      source = "Observed"
    )
  
  bind_rows(y_hist, y_fc) |>
    mutate(
      Region = region_label,
      Variable = component_label
    )
}

## ------------------------------------------------------------
## 4. Helper: climate historical + simple conditional projection
## ------------------------------------------------------------
## For climate variables, use ARIMA-free persistence/trend approach:
## - Temperature: linear trend extrapolation from historical annual series.
## - SPEI: mean-reverting projection to historical mean.
##
## These are not climate forecasts; they are conditional paths for
## visualising model-aligned adjustment.

project_temperature <- function(region, region_label, temp_series, horizon = 29) {
  
  hist_df <- tibble::tibble(
    Year = temp.w$Year,
    value = temp.w[[temp_series]]
  ) |>
    dplyr::filter(Year >= 1961, Year <= 2021) |>
    tidyr::drop_na()
  
  fit <- lm(value ~ Year, data = hist_df)
  
  fut_df <- tibble::tibble(
    Year = 2022:(2021 + horizon)
  )
  
  fut_df <- fut_df |>
    dplyr::mutate(
      value = as.numeric(predict(fit, newdata = fut_df)),
      source = "Projection"
    )
  
  hist_df <- hist_df |>
    dplyr::mutate(source = "Observed")
  
  dplyr::bind_rows(hist_df, fut_df) |>
    dplyr::mutate(
      Region = region_label,
      Variable = "Temperature anomaly"
    )
}

project_spei <- function(region, region_label, spei_series, horizon = 29) {
  
  hist_df <- tibble(
    Year = spei.w$Year,
    value = spei.w[[spei_series]]
  ) |>
    filter(Year >= 1961, Year <= 2021) |>
    drop_na()
  
  spei_mean <- mean(hist_df$value, na.rm = TRUE)
  last_val <- dplyr::last(hist_df$value)
  
  # Simple gradual reversion to historical mean
  fut_vals <- purrr::map_dbl(1:horizon, function(h) {
    spei_mean + (last_val - spei_mean) * (0.85 ^ h)
  })
  
  fut_df <- tibble(
    Year = 2022:(2021 + horizon),
    value = fut_vals,
    source = "Projection"
  )
  
  hist_df <- hist_df |>
    mutate(source = "Observed")
  
  bind_rows(hist_df, fut_df) |>
    mutate(
      Region = region_label,
      Variable = "SPEI"
    )
}

## ------------------------------------------------------------
## 5. Build plotting dataset for one region
## ------------------------------------------------------------

build_region_projection_data <- function(region,
                                         region_label,
                                         temp_series,
                                         spei_series,
                                         horizon = 29) {
  
  climate_df <- bind_rows(
    project_temperature(region, region_label, temp_series, horizon),
    project_spei(region, region_label, spei_series, horizon)
  ) |>
    mutate(
      lower = NA_real_,
      upper = NA_real_
    )
  
  component_df <- purrr::map_dfr(names(components), function(v) {
    forecast_component(
      varname = v,
      component_label = components[[v]],
      region = region,
      region_label = region_label,
      horizon = horizon
    )
  })
  
  bind_rows(climate_df, component_df) |>
    mutate(
      Variable = factor(
        Variable,
        levels = c(
          "Temperature anomaly",
          "SPEI",
          "Capital-labour ratio",
          "Crop-pasture ratio",
          "Livestock",
          "Fertiliser"
        )
      )
    )
}

## ------------------------------------------------------------
## 6. Plot function
## ------------------------------------------------------------

plot_region_projection <- function(plot_df, region_label) {
  
  ggplot(plot_df, aes(x = Year, y = value, colour = source)) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_vline(xintercept = 2021, linetype = "dashed", colour = "grey40") +
    facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
    scale_colour_manual(
      values = c(
        "Observed" = "black",
        "Projection" = "#0072B2"
      )
    ) +
    labs(
      title = paste0("Conditional adjustment projections: ", region_label),
      subtitle = "Observed series to 2021; conditional projections to 2050 based on corrected VECM specification",
      x = NULL,
      y = NULL,
      colour = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

## ------------------------------------------------------------
## 7. Generate and save Figures 3–6
## ------------------------------------------------------------

figures_out <- purrr::pmap(
  list(
    figure_specs$region,
    figure_specs$region_label,
    figure_specs$temp_series,
    figure_specs$spei_series
  ),
  function(region, region_label, temp_series, spei_series) {
    
    df <- build_region_projection_data(
      region = region,
      region_label = region_label,
      temp_series = temp_series,
      spei_series = spei_series,
      horizon = 29
    )
    
    p <- plot_region_projection(df, region_label)
    
    file_stub <- gsub("[^A-Za-z0-9]+", "_", region_label)
    
    ggsave(
      filename = paste0("Figure_", file_stub, "_conditional_projection.png"),
      plot = p,
      width = 8,
      height = 7,
      dpi = 300
    )
    
    list(data = df, plot = p)
  }
)

names(figures_out) <- figure_specs$region_label

## Print plots interactively
figures_out[["North America"]]$plot
figures_out[["Central Europe"]]$plot
figures_out[["Sub-Saharan Africa"]]$plot
figures_out[["Australia"]]$plot


