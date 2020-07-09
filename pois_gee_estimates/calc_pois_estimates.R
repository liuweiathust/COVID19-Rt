library(dplyr)
library(splines)
library(geepack)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)

source("../preprocess_data/preprocess_jhu.R")
source("../helpers/helper_functions.R")
source('../estimate_rt/estimate_rt_master.R')

na_min <- function(x) {
  ret <- NA
  if (length(x) > 0) {
    ret <- min(x)
  }
  return(ret)
}

prep_data <- function(dt_input, byvar = "UID") {
  dt <- copy(dt_input)
  dt[, `:=` (
    # days since start of 2020
    days_cal = days_since("2020-01-01", date),
    # days since >50 reported cases
    days_50case = days_since(na_min(date[positive >= 50]), date)),
    by = byvar]

  dt[, `:=` (
    weekday = weekdays(date),
    # set negative daily increase values to 0
    positiveIncrease_trunc = positiveIncrease * (positiveIncrease >= 0))]

  # sort before getting lambda
  setorderv(dt, cols = c(byvar, "date"))
  dt[, Lambda := getLambda(x = positiveIncrease_trunc, t = days_cal,
                        mu = 5.2, sd = 5.5, NTS = 30), by = byvar]
  dt[,`:=` (Lambda_sum7 = lagSum(Lambda, 7),
            positiveIncrease_sum7 = lagSum(positiveIncrease_trunc, 7)),
     by = byvar]
  return(dt)
}


read_jhu <- function(level = c("County", "State", "Global"), start_date) {
  level <- match.arg(level)
  jhu_url <- paste0("https://raw.githubusercontent.com/lin-lab/COVID-data-cleaning/master/jhu_data/cleaned_data/JHU_COVID-19_", level, ".csv")

  dat <- fread(jhu_url)
  dat[, date := as.Date(date, format = "%Y-%m-%d")]
  ret <- dat[date >= start_date, ]
  return(ret)
}

fit_rt <- function(dt, days_touse = c("days_cal", "days_50case"), n_knots = 30,
                   corstr = "ar1", conf_level = 0.95) {
  days_touse <- match.arg(days_touse)
  day_knots <- seq_days(dt[[days_touse]], n_knots)
  if (days_touse == "days_cal") {
    model_formula <- positiveIncrease_sum7 ~ offset(log1p(Lambda_sum7)) +
      factor(weekday) + bs(days_cal, knots = day_knots)
  } else if (days_touse == "days_50case") {
    model_formula <- positiveIncrease_sum7 ~ offset(log1p(Lambda_sum7)) +
      factor(weekday) + bs(days_50case, knots = day_knots)
  }

  setorder(dt, date)
  fit <- geepack::geeglm(formula = model_formula, data = dt,
                         id = UID, family = poisson(link = "log"),
                         corstr = corstr)
  x <- model.matrix(~ bs(dt$days_cal, knots = day_knots))
  coef_names <- names(coef(fit))
  coef_idx <- startsWith(coef_names, "bs") | coef_names == "(Intercept)"
  coef_use <- coef(fit)[coef_idx]
  log_Rt <- x %*% coef_use
  log_Rt_se <- sqrt(diag(x %*% tcrossprod(vcov(fit)[coef_idx, coef_idx], x)))

  zstar <- qnorm(1 - (1 - conf_level) / 2)
  ret <- data.frame(rt = exp(log_Rt),
                    ci_lower = exp(log_Rt - zstar * log_Rt_se),
                    ci_upper = exp(log_Rt + zstar * log_Rt_se),
                    date = dt$date)
  return(ret)
}

start_date <- as.Date("2020-03-01", format = "%Y-%m-%d")
jhu_counties <- read_jhu("County", start_date = start_date) %>%
  prep_data()
jhu_states <- read_jhu("State", start_date = start_date) %>%
  prep_data()
jhu_global <- read_jhu("Global", start_date = start_date) %>%
  prep_data()

ma_data <- jhu_states[stateName == "Massachusetts" &
                      date >= as.Date("2020-03-19", format = "%Y-%m-%d")]

ma_data <- jhu_global[Combined_Key == "Equatorial Guinea"]

rt_poisgee_50_ma <- fit_rt(ma_data, days_touse = "days_50case") %>%
  mutate(fit = "GEE, days since 50th case")
rt_poisgee_cal_ma <- fit_rt(ma_data, days_touse = "days_cal") %>%
  mutate(fit = "GEE, calendar day")
rt_epiestim_ma <- estimate_rt_EpiEstim(ma_data$date, ma_data$positiveIncrease,
                                       mean_serial = 5.2, std_serial = 5.1) %>%
  dplyr::select(rt = mean_rt, ci_lower, ci_upper, date = interval_end) %>%
  mutate(fit = "EpiEstim")

compare_df <- rbind(rt_poisgee_50_ma, rt_poisgee_cal_ma, rt_epiestim_ma)

ggplot(compare_df, aes(x = date, y = rt, ymin = ci_lower, ymax = ci_upper,
                       color = fit, fill = fit)) +
  geom_line() + geom_ribbon() +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("Date") + ylab("Rt") +
  ggtitle("Rt for MA") +
  theme_cowplot()
