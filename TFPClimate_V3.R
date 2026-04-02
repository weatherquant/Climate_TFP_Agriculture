## AGRICULTURE TFP - CLIMATE ANALYSIS
# Core time-series econometrics
library(urca)      # Johansen cointegration, ca.jo
library(vars)      # VAR / VECM tools
library(tseries)   # ADF tests

# Data handling
library(readr)
library(dplyr)
library(tidyr)

# Plotting (optional)
library(ggplot2)

# 1.1 Load datasets
tfp    <- read_csv("TFP_indexes.csv")
output <- read_csv("Agricultural_Growth.csv")
inputs <- read_csv("Ag_Indexes.csv")
temp   <- read_csv("Global_Land_Temp.csv")
spei   <- read_csv("Main_Spei.csv")

# 1.2 Harmonise years and regions
regions <- c("N. America", "C. Europe", "SSA", "Australia")

# Restrict to common sample (1961–2021)
years <- 1961:2021

tfp    <- tfp    |> filter(Year %in% years)
output <- output |> filter(Year %in% years)
inputs <- inputs |> filter(Year %in% years)
temp   <- temp   |> filter(Year %in% years)
spei   <- spei   |> filter(Year %in% years)

temp.w <- temp %>%
  pivot_wider(names_from = Location, values_from = Value)

spei.w <- spei %>%
  pivot_wider(names_from = Country, values_from = SPEI)

colnames(tfp)<- c("Year","SSA","CentralEurope","NorthAmerica","Australia")
colnames(temp.w)<- c("Year","Global","Australia","Canada","NorthAmerica","SSA","Brazil","CentralEurope")
colnames(spei.w)<- c("Year","Australia","CentralEurope","NorthAmerica","Brazil","Canada","SSA","Nigeria")
colnames(inputs)<- c("Country","Region","Subregion","Inc","Year","Attribute","Value")

inputs.w <- inputs %>%
  pivot_wider(names_from = Attribute, values_from = Value) %>%
  dplyr::select(!c(Region,Subregion,Inc)) 

inputs.w <- as.data.table(inputs.w)
inputs.w <- inputs.w[Country %in% c("Australia","United States", "Europe, Central", "SSA, Total"),]
inputs <- melt.data.table(inputs.w, id.vars = c("Country","Year"))
inputs$Country[inputs$Country == "United States"] <- "NorthAmerica"
inputs$Country[inputs$Country == "SSA, Total"] <- "SSA"
inputs$Country[inputs$Country == "Europe, Central"] <- "CentralEurope"

temp$Location[temp$Location == "US"] <- "NorthAmerica"
spei$Country[spei$Country == "US"] <- "NorthAmerica"
spei$Country[spei$Country == "Central_Eur"] <- "CentralEurope"
spei$Country[spei$Country == "Sudan"] <- "SSA"

# 2. Pre‑estimation checks (ADF tests)
adf_test <- function(x, name){
  cat("\nADF test:", name, "\n")
  print(adf.test(na.omit(x)))
}

# ADF stationary checks
# Australia
adf_test(tfp$Australia, "TFP Australia")
adf_test(temp$Location=="Australia", "Temperature Australia")
adf_test(spei$Country=="Australia", "SPEI Australia")

# US
adf_test(tfp$NorthAmerica, "TFP Nth America")
adf_test(temp$Location=="NorthAmerica", "Temperature Nth America")
adf_test(spei$Country=="NorthAmerica", "SPEI Nth America")

# Central Europe
adf_test(tfp$CentralEurope, "TFP Central Europe")
adf_test(temp$Location=="CentralEurope", "Temperature Central Europe")
adf_test(spei$Country=="CentralEurope", "SPEI Central Europe")

# SSA
adf_test(tfp$SSA, "TFP SSA")
adf_test(temp$Location=="SSA", "Temperature SSA")
adf_test(spei$Country=="SSA", "SPEI SSA")

# 3. VECM estimation — TFP models (Table 1)
estimate_vecm_tfp <- function(region){
  
  df <- data.frame(
    TFP  = log(tfp[[region]]),
    TEMP = temp.w[[region]],
    SPEI = spei.w[[region]]
  )
  
  df <- ts(df, start = 1961, frequency = 1)
  
  # Lag length selection
  lag_sel <- VARselect(df, lag.max = 8, type = "const")
  p <- max(2,lag_sel$selection["AIC(n)"])
  
  # Johansen test
  joh <- ca.jo(df,
               type = "trace",
               ecdet = "const",
               K = p)
  
  # Estimate VECM (r = cointegrating rank)
  vecm <- cajorls(joh, r = 2)
  
  return(list(johansen = joh,
              vecm = vecm))
}

# 3.1a Extract coefficients
extract_ect <- function(vecm_obj){
  
  coef_mat <- summary(vecm_obj$vecm$rlm)#$coefficients
  
  coef_mat <- coef_mat[["Response TFP.d"]][["coefficients"]]
  
  ect_rows <- coef_mat[grep("ect", rownames(coef_mat)), ]
  
  return(ect_rows[, c("Estimate", "Std. Error", "Pr(>|t|)")])
}

# 3.2 Run for all regions
extract_ect(estimate_vecm_tfp("Australia"))
extract_ect(estimate_vecm_tfp("NorthAmerica"))
extract_ect(estimate_vecm_tfp("CentralEurope"))
extract_ect(estimate_vecm_tfp("SSA"))

##############################################
#### TO REPRODUCE TABLE 1:

regions <- c("Australia", "NorthAmerica", "CentralEurope", "SSA")

extract_ect_table <- function(vecm_obj) {
  
  coef_mat <- summary(vecm_obj$vecm$rlm)#$coefficients
  
  coef_mat <- coef_mat[["Response TFP.d"]][["coefficients"]]
  
  ects <- coef_mat[grep("^ect", rownames(coef_mat)), ]
  
  stars <- ifelse(
    ects[, "Pr(>|t|)"] < 0.01, "***",
    ifelse(ects[, "Pr(>|t|)"] < 0.05, "**",
           ifelse(ects[, "Pr(>|t|)"] < 0.10, "*", ""))
  )
  
  formatted <- paste0(
    sprintf("%.3f", ects[, "Estimate"]),
    stars,
    " (",
    sprintf("%.3f", ects[, "Std. Error"]),
    ")"
  )
  
  names(formatted) <- rownames(ects)
  
  return(formatted)
}

# Run TFP models
tfp_models <- lapply(regions, estimate_vecm_tfp)
names(tfp_models) <- regions

# Extract ECTs
tfp_table <- lapply(tfp_models, extract_ect_table)

# Assemble Table 1
Table1_TFP <- data.frame(
  Region = c("Australia",
             "North America",
             "Central Europe",
             "Sub-Saharan Africa"),
  ECT1 = c(tfp_table$Australia["ect1"],
           tfp_table$NorthAmerica["ect1"],
           tfp_table$CentralEurope["ect1"],
           tfp_table$SSA["ect1"]),
  ECT2 = c(tfp_table$Australia["ect2"],
           tfp_table$NorthAmerica["ect2"],
           tfp_table$CentralEurope["ect2"],
           tfp_table$SSA["ect2"])
)

print(Table1_TFP)

# 4. Production‑component models (Table 2)
# Repeat the same logic, variable by variable.
# 4.1 Generic production ECM function

# Ensure proper types
inputs <- inputs %>%
  mutate(
    Year = as.integer(Year),
    value = as.numeric(value)
  )

estimate_vecm_component <- function(varname, region,
                                    start_year = 1961,
                                    end_year   = 2021,
                                    max_lag    = 6) {
  
  # --- Extract production component ---
  y_df <- inputs %>%
    filter(
      Country  == region,
      variable == varname,
      Year >= start_year,
      Year <= end_year
    ) %>%
    arrange(Year)
  
  if (nrow(y_df) == 0) {
    stop(paste("No data for", varname, "in", region))
  }
  
  # --- Extract climate data ---
  t_df <- temp %>%
    filter(
      Location == region,
      Year >= start_year,
      Year <= end_year
    ) %>%
    arrange(Year)
  
  s_df <- spei %>%
    filter(
      Country == region,
      Year >= start_year,
      Year <= end_year
    ) %>%
    arrange(Year)
  
  # --- Merge into single data frame ---
  df <- data.frame(
    Y    = log(y_df$value),
    TEMP = t_df$Value,
    SPEI = s_df$SPEI
  )
  
  # --- Convert to time series object ---
  df_ts <- ts(df, start = start_year, frequency = 1)
  
  # --- Lag selection (AIC) ---
  lag_sel <- VARselect(df_ts, lag.max = max_lag, type = "const")
  p <- max(2,lag_sel$selection["AIC(n)"])
  
  # --- Johansen cointegration test ---
  joh <- ca.jo(
    df_ts,
    type  = "trace",
    ecdet = "const",
    K     = p
  )
  
  # --- Estimate VECM ---
  vecm <- cajorls(joh, r = 2)
  
  return(list(
    johansen = joh,
    vecm     = vecm,
    lags     = p
  ))
}

#############
# REPRODUCE TABLE 2:

#########################################
## LIVESTOCK ############################
regions <- c("Australia", "NorthAmerica", "CentralEurope", "SSA")

livestock_models <- list()

for (r in regions) {
  livestock_models[[r]] <- estimate_vecm_component(
    varname = "Livestock_Q",
    region  = r
  )
  
  cat("\n============================\n")
  cat("Livestock ECM –", r, "\n")
  cat("============================\n")
  
  print(summary(livestock_models[[r]]$vecm$rlm))
}

#########################################
## FERTILISER ############################
regions <- c("Australia", "NorthAmerica", "CentralEurope", "SSA")

fert_models <- list()

for (r in regions) {
  fert_models[[r]] <- estimate_vecm_component(
    varname = "Fertilizer_Q",
    region  = r
  )
  
  cat("\n============================\n")
  cat("Fertiliser ECM –", r, "\n")
  cat("============================\n")
  
  print(summary(fert_models[[r]]$vecm$rlm))
}

#############################################################
# compute coefficients for the ratio indexes (cap-lab and crop-past)
# 1. Reshape to wide to compute ratio
inputs <- inputs %>%
  pivot_wider(
    names_from  = variable,
    values_from = value
  )

# 2. Compute capital–labour ratio
inputs <- inputs %>%
  mutate(
    CapLab_Index = Capital_Index / Labor_Index
  ) %>% 
  mutate(
    CropPast_Index = Cropland_Q / Pasture_Q
  )

# 3. Return to long format and append
inputs <- inputs %>%
  pivot_longer(
    cols      = c(CapLab_Index,CropPast_Index),
    names_to  = "variable",
    values_to = "value"
  )

## CAPITAL-LABOUR ############################
regions <- c("Australia", "NorthAmerica", "CentralEurope", "SSA")

caplab_models <- list()

for (r in regions) {
  caplab_models[[r]] <- estimate_vecm_component(
    varname = "CapLab_Index",
    region  = r
  )
  
  cat("\n============================\n")
  cat("Capital-Labour ECM –", r, "\n")
  cat("============================\n")
  
  print(summary(caplab_models[[r]]$vecm$rlm))
}

## CROP-PASTURE ############################
regions <- c("Australia", "NorthAmerica", "CentralEurope", "SSA")

croppast_models <- list()

for (r in regions) {
  croppast_models[[r]] <- estimate_vecm_component(
    varname = "CropPast_Index",
    region  = r
  )
  
  cat("\n============================\n")
  cat("Crop-Pasture ECM –", r, "\n")
  cat("============================\n")
  
  print(summary(croppast_models[[r]]$vecm$rlm))
}

#############
# REPRODUCE TABLE 2:
regions <- c("Australia", "NorthAmerica", "CentralEurope", "SSA")

r<-regions[1]
region<-regions[1]

livestock_models[[r]]$vecm$rlm
fert_models[[r]]$vecm$rlm
caplab_models[[r]]$vecm$rlm
croppast_models[[r]]$vecm$rlm
vecm$rlm
summary(vecm$rlm)

livestock_models[[r]]$vecm$rlm

extract_table2 <- function(rlm_obj, region, component) {
  
  summ <- summary(rlm_obj)
  
  # Identify the production (Y) equation safely
  eq_names <- names(summ)
  
  y_eq <- eq_names[grepl("^Response\\s+Y", eq_names)]
  
  if (length(y_eq) == 0) {
    stop(paste("No Response Y equation found for", component, "in", region))
  }
  
  # Extract coefficient matrix for ΔY equation
  coef_mat <- summ[[y_eq[1]]]$coefficients
  
  # Extract error-correction terms only
  ect_rows <- grep("^ect", rownames(coef_mat), value = TRUE)
  
  if (length(ect_rows) == 0) {
    stop(paste("No ECT terms found for", component, "in", region))
  }
  
  tibble::tibble(
    Region    = region,
    Component = component,
    ECT       = ect_rows,
    Estimate  = as.numeric(coef_mat[ect_rows, "Estimate"]),
    StdError  = as.numeric(coef_mat[ect_rows, "Std. Error"]),
    p_value   = as.numeric(coef_mat[ect_rows, "Pr(>|t|)"])
  )
}

livestock_table <- purrr::map_dfr(
  regions,
  ~ extract_table2(livestock_models[[.x]]$vecm$rlm, .x, "Livestock")
)

fertiliser_table <- purrr::map_dfr(
  regions,
  ~ extract_table2(fert_models[[.x]]$vecm$rlm, .x, "Fertiliser")
)

caplab_table <- purrr::map_dfr(
  regions,
  ~ extract_table2(caplab_models[[.x]]$vecm$rlm, .x, "Capital–Labour")
)

croppast_table <- purrr::map_dfr(
  regions,
  ~ extract_table2(croppast_models[[.x]]$vecm$rlm, .x, "Crop–Pasture")
)

table2 <- dplyr::bind_rows(
  caplab_table,
  croppast_table,
  livestock_table,
  fertiliser_table
)

print(table2)


table2_dom <- table2 %>%
  group_by(Region, Component) %>%
  slice_max(abs(Estimate), n = 1, with_ties = FALSE) %>%
  ungroup()

table2_dom <- table2_dom %>%
  mutate(
    Stars = case_when(
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

table2_final <- table2_dom %>%
  dplyr::select(Component, Region, Cell) %>%
  tidyr::pivot_wider(
    names_from  = Region,
    values_from = Cell
  )

print(table2_final)
###
## End

