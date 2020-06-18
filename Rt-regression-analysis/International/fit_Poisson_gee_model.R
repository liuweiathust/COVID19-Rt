#!/usr/bin/R

#----------------------------------------------------------------
# Project: Lin Lab - Covid19
# Analysis: Country-level Rt regression
# Description: Poisson GEE direct fit of cases, offset(Lambda)
# Author: Hui Li
# Requirements: cleaned country level data
#----------------------------------------------------------------

rm(list=ls())

# data management
library(data.table)
library(dplyr)
library(tidyr)
library(car)

# graphing
library(ggplot2)
library(grid)

# model fitting
library(geepack)
library(lme4)
library(splines)
library(usdm)
library(nlme)

#### ===========================
#### LOAD THE DATA
#### ===========================
# Load code to calculate cumulative incidence.
source("CalculateCumulativeIncidence.R")

# directories setup (if needed)

# load cleaned merged data
dt <- as.data.frame(fread("cleaned_country_level_2020_06_18.csv", 
                          sep=",", header = TRUE, fill = TRUE))

#### ======================
#### HELPER FUNCTIONS
#### ======================
# calculate days since an event with specified start_date
days_since <- function(start_date, dates, nz = TRUE){
  out <- as.numeric(as.Date(dates, format = "%Y-%m-%d") - as.Date(start_date, format = "%Y-%m-%d"))
  if(nz){out[out < 0] <- 0}
  return(out)
}

# Generate spline knots, fixed interval of days apart
seq_days <- function(x, delta){
  out <- seq(min(x) + delta, max(x), delta)
}

# function to sum over sliding windows
lagSum <- function(x, n){
  csx = cumsum(x)
  csx - data.table::shift(csx, n = n, type = 'lag', fill = 0)
}

# function to average over sliding windows
lagMean <- function(x, n){
  csx = cumsum(x)
  (csx - data.table::shift(csx, n = n, type = 'lag', fill = 0))/n
}

# Reverse lower case country names: for plotting aesthetics
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

# Number of ticks
number_ticks <- function(n) {function(limits) pretty(limits, n)}

#### ===============
#### GEE functions
#### ===============
# Fit overdispersed Poisson model
# Betas from this model can be interpreted as 
# differences in log(Rt). 
geeglm_fit_pois <- function(dt, model_formula, gee_id_var, corstr = "ar1",
                            estimate_type = "geese"){
  
  #' @param dt Data object for fitting models.
  #' @param model_formula Model formula for GEE GLM.
  #' @param gee_id_var Cluster variable name.
  #' @param corstr Correlation structure type.
  #' @param estimate_type What estimation method to use: geese is used by default. 
  
  # mold gee-id-var and sort by it
  dt <- as.data.frame(dt)
  dt$gee_id_var <- as.integer(factor(dt[[gee_id_var]], order = TRUE))
  dt <- as.data.table(dt)
  setorder(dt, gee_id_var)
  
  if (estimate_type == "glm") {
    fit <- geepack::geeglm(formula = model_formula,
                           data = dt, 
                           id = gee_id_var, 
                           family = poisson(),
                           corstr = corstr,
                           control = geese.control("epsilon" = 1e-4, 
                                                   "maxit" = 50000, 
                                                   trace = TRUE))
  } else if (estimate_type == "geese") {
    fit <- geepack::geese(formula = model_formula,
                          data = dt, 
                          id = gee_id_var, 
                          family = poisson(),
                          corstr = corstr,
                          control = geese.control("epsilon" = 1e-4, 
                                                  "maxit" = 50000, 
                                                  trace = TRUE))
  }
  
  return(fit)
}

# Wrapper function: Extract regression coefficients and prepare for bar plots
geeglm_coef_data <- function(fit, var_type, gee_obj, variable_name = NULL, estimate_name = "Estimate") {
  
  #' @param fit The model fit object, either geese or glm, specified by gee_obj.
  #' @param var_type Type of variance to use. 
  #' @param gee_obj What type of object is "fit". Slightly different cleaning steps for these two types.
  #' @param variable_name Labels of the variable.
  #' @param estimate_name The name of the coefficient estimates in the regression summary.
  
  if (gee_obj == "glm") {
    # extract coefficient from the fit object
    coef_table <- as.data.table(coef(summary(fit)), keep = "Variable")
    
    # filter out factor variable and time variables
    ind <- which(coef_table$Variable %like% "scale" | 
                   coef_table$Variable %like% "lockdown" | 
                   coef_table$Variable %like% "continent" | 
                   # coef_table$Variable %like% "weekend" |
                   coef_table$Variable %like% "positive")
    
    y_se <- sqrt(diag(fit$geese[[var_type]]))[ind]
    
  } else if (gee_obj == "geese") {
    # extract coefficient from the fit object
    coef_table <- as.data.table(fit$beta, keep = "Variable")
    colnames(coef_table) <- c("Variable", estimate_name)
    
    # filter out factor variable and time variables
    ind <- which(coef_table$Variable %like% "scale" | 
                   coef_table$Variable %like% "lockdown" | 
                   coef_table$Variable %like% "continent" | 
                   # coef_table$Variable %like% "weekend" |
                   coef_table$Variable %like% "positive")
    
    y_se <- sqrt(diag(fit[[var_type]]))[ind]
    
  }
  
  coef_table <- coef_table %>% filter(Variable %like% "scale" |
                                        Variable %like% "lockdown" | 
                                        Variable %like% "continent" | 
                                        # Variable %like% "weekend" |
                                        coef_table$Variable %like% "positive")
  
  # clean up coefficient names and extract corresponding s.e. values 
  y <- coef_table[, estimate_name]
  
  if(is.null(variable_name)){
    coef_name <- gsub(pattern = "scale", replacement="",
                      x = coef_table$Variable)
    coef_name <- gsub(pattern="\\)", replacement="",
                      x = gsub(pattern="\\(", replacement="",
                               x = coef_name))
  }else{
    coef_name <- variable_name
  }
  
  # collect coefficients and s.e.
  plot_dt <- data.frame(name=coef_name, 
                        value=y, 
                        se=y_se)
  return(plot_dt)
}


#### ===========================
#### PREPARE REGRESSION DATA
#### ===========================

# define continuous time variables
dt <- data.table(dt)
dt[,`:=`(
  # days since start of 2020
  days_int = days_since("2020-01-01", date),
  
  # days since lockdown order was placed
  days_s_lockdown = days_since(lockdown_date, date),

  # days since first reported case
  days_s_1p = days_since(min(date[positive > 0]), date),
  
  # days since >50 reported cases
  days_s_50p = days_since(min(date[positive >= 50]), date)
), by = countryName]


# define dummy time variables
dt[,`:=`(
  lockdown = 1*(days_s_lockdown > 0),
  weekday = weekdays(as.Date(date))
),]


# Factorize variables for regression
dt <- as.data.frame(dt)
dt$lockdown <- as.factor(dt$lockdown)
dt$countryName <- as.factor(dt$countryName)
dt$weekday <- as.factor(dt$weekday)

# Add weekend indicator variable
dt$weekend <- ifelse(dt$weekday %in% c("Saturday", "Sunday"), 1, 0)

# sort the data before getting Lambda
setorder(dt, countryName, date)

# Set any negative daily increase values to 0.
dt <- as.data.table(dt)
dt[, daily_inc := positiveIncrease*(positiveIncrease>=0),]

# Calculate the cumulative daily incidence
# use serial interval mean=5.2, sd = 5.5 to 
# match Sheila and Andy's analysis.
dt[, Lam := getLambda(
  x = daily_inc,
  t = days_int, 
  mu = 5.2, sd = 5.5, NTS = 30
), by = countryName]


# To mimick EpiEstim, aggregate within 7-day window:
dt[,`:=`(
  Lam_s7 = lagSum(Lam, 7),
  inc_s7 = lagSum(daily_inc, 7),
  mobility_avg7 = lagMean(mobility_change_residential, 7)
), by = countryName]

# cross sectional data (use latest observation)
dt_country <- (dt) %>%
  group_by(countryName) %>%
  arrange(countryName, desc(date)) %>%
  summarise_all(first)
# 184 countries --> 166 countries


# Quality control of Rt
# Option 1: truncate Rt values (keep 1-99 percentile)
# dt <- (dt) %>%
#   filter(!is.na(mean_rt)) %>%
#   filter(mean_rt < quantile(mean_rt, 0.99) & mean_rt > quantile(mean_rt, 0.01))

# Option 2: Restrict country-dates to total cases > 50 <-- Avoid too much weighting on priors
dt <- dt[dt$positive >= 50, ]

# Remove NA values for RHS (ONLY KEEP TERMS TO BE INCLUDED)
dt <- as.data.frame(dt)
dt_reg <- dt[complete.cases(dt[c("daily_inc", "positive", 
                                 "Lam", 
                                 "days_int", "days_s_50p",
                                 "weekday", 
                                 "unique_continent",
                                 "population", 
                                 "population_density", 
                                 "median_age", 
                                 "pct_diabetes", 
                                 "cvd_death_rate", 
                                 "Physicians_per_1000", 
                                 "hospital_beds_per_100k", 
                                 "Health_exp_pct_GDP_2016",
                                 "gdp_per_capita_PPP",
                                 # "Life_expectancy_at_birth_both",
                                 "days_s_lockdown",
                                 "lockdown",
                                 "stringencyIndex",
                                 "mobility_avg7"
                                  )]), ]

# 82 countries, 6594 observations

# Compare calculations to EpiEstim.
# Calculated Rt is positively biased
dt_reg <- as.data.frame(dt_reg %>%
                group_by(countryName) %>%
                arrange(countryName, date) %>%
                mutate(days_within_country = sequence(n())))
                            
# png(file = paste0( "Rt_compare.png"))
# plot(I((inc_s7)/(Lam_s7+1)) ~ mean_rt, 
#      subset(dt_reg, days_within_country > 7 & Lam_s7 > 50 & positive >= 50),
#      xlab = "EpiEstim Rt", 
#      ylab = "Hand-Calculated Rt"); 
# abline(0, 1, col='red')
# dev.off()

# # Investigate the deviations and calibrate the distance
# dt_reg$rt_diff <- abs((dt_reg$inc_s7 / (dt_reg$Lam_s7+1)) - dt_reg$mean_rt)
# aa <- dt_reg[dt_reg$days_within_country > 1 & dt_reg$Lam_s7 > 50 & dt_reg$positive >= 50, ]
# denom_count <- dim(aa)[1]
# 1 - dim(aa[aa$rt_diff > 0.1, ])[1] / denom_count
# 1 - dim(aa[aa$rt_diff > 0.5, ])[1] / denom_count
# 1 - dim(aa[aa$rt_diff > 1, ])[1] / denom_count

# 96% within 0.1 
# 98.8% within 0.5 
# 99.3% within 1

# subset to tolerance within 0.5
dt_reg <- dt_reg[dt_reg$days_within_country > 7 & dt_reg$positive >= 50, ]

# 80 countries, 6026 observations

#### ===============
#### MODEL FITTING
#### ===============
# Fixed Effects
fe_eqn = daily_inc ~ 
  offset(log1p(Lam)) +
  log1p(positive) +
  factor(weekday) +
  factor(unique_continent) + 
  scale(population) +
  scale(population_density) +
  scale(median_age) +
  scale(pct_diabetes) +
  scale(cvd_death_rate)+ 
  scale(Physicians_per_1000) +
  scale(hospital_beds_per_100k) +
  scale(Health_exp_pct_GDP_2016) +
  scale(gdp_per_capita_PPP) +
  # scale(Life_expectancy_at_birth_both) +
  scale(days_s_lockdown) +
  factor(lockdown) +
  scale(stringencyIndex) + 
  scale(mobility_avg7)

# check the date range
as.Date("2020-01-01")+min(dt_reg$days_int) #03-15
as.Date("2020-01-01")+max(dt_reg$days_int) #06-07

# check VIF of the model
check_vif <- vif(lm(fe_eqn, data = dt_reg))

###### Model 1. Uncentered time (Calendar date)

# knots values for B-spline terms
day_knots <- seq_days(dt_reg$days_int, 30)

# clean up final equation for fitting
eqn = as.formula(paste0(fe_eqn[2], " ~ ", fe_eqn[3],
                        "+ bs(days_int, knots = day_knots)"))

# run geeglm
fit_uncenter <- geeglm_fit_pois(dt = dt_reg, 
                       model_formula = eqn,
                       gee_id_var = "countryName", 
                       corstr = "ar1", 
                       estimate_type = "glm")

###### Model 2. Centered time (Days since outbreak, 50 cases reported within country)

# knots values for B-spline terms
day_knots <- seq_days(dt_reg$days_s_50p, 30)

# clean up final equation for fitting
eqn = as.formula(paste0(fe_eqn[2], " ~ ", fe_eqn[3],
                        "+ bs(days_s_50p, knots = day_knots)"))

# run geeglm
fit_center <- geeglm_fit_pois(dt = dt_reg, 
                                model_formula = eqn,
                                gee_id_var = "countryName", 
                                corstr = "ar1", 
                                estimate_type = "glm")

#### ======================
#### REGRESSION PLOTS
#### ======================
# manually enter covariate names
variable_name <- c(
   # "Weekend or not",
   "Log1p(Total cases)",
   "Ref-Africa: Asia", 
   "Ref-Africa: Europe", 
   "Ref-Africa: North America",
   "Ref-Africa: Oceania", 
   "Ref-Africa: South America",
   "Ref-Africa: Trans Asia and Europe",
   "Population", 
   "Population density",
   "Median age",
   "Diabetes prevalence", 
   "CVD Death rate",
   "# of Physicians per 1000", 
   "# of hospital beds per 100k",
   "Health expenditure as % of GDP",
   "GDP per capita (PPP)",
   # "Life expectancy at birth",
   "Days since lockdown",
   "Lockdown or not",
   "Government policy stringency",
   "Mobility change (lag 7 average)"
   )

# plot for both S.E. options
variance_types = c("vbeta", "vbeta.naiv")
tags <- c("san", "naive")
se_labels <- c("Sandwich", "Model-based")

# plot for these models
results <- list(fit_uncenter, fit_center)
model_labels <- c("Uncentered Time (Calendar Date)", 
                  "Centered Time (Days since 50 cases reported in the country)")
model_tags <- c("fit_uncenter_", "fit_center_")

for (j in 1:2) {
  for (i in 1:2) {
    # 1. Plot results from geeglm
    plot_dt <- geeglm_coef_data(results[[j]], var_type = variance_types[i], variable_name = variable_name, gee_obj = "glm")
    p1 <- ggplot(plot_dt) +
      geom_bar(aes(x=name, y=value), stat="identity", fill="royalblue", alpha=0.7) +
      geom_errorbar(aes(x=name, ymin=value-(1.96*se), ymax=value+(1.96*se), width=0.2),
                    colour="firebrick2", alpha=0.9, size=0.8) +
      coord_flip() +
      labs(x = 'Variables', y = expression("Coefficient" %+-% "1.96 SE"), 
           title = 'Poisson GEE Regression Coefficients',
           subtitle = paste0("80 countries, 6026 observations. March 15 - June 7. \nCountry clusters, AR1 correlation structure. ", 
                             se_labels[i], " SE. \n", model_labels[j])) +
      theme(plot.title=element_text(hjust=0.5, size=30),
            plot.subtitle=element_text(size=15),
            axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0), hjust=0.5, size=15),
            axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0), vjust=0.5, size=15),
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=15),
            legend.title=element_text(size=15),
            legend.text=element_text(size=15))
    ggsave(paste0("pois_gee_", model_tags[j], tags[i], "_se_reg_coef_bar.png"), p1, width=15, height=10)
  }
}

