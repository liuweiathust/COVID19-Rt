library(dplyr)
library(splines)
library(magrittr)
library(data.table)
library(purrr)

source("../preprocess_data/preprocess_jhu.R")
source("../helpers/helper_functions.R")

read_jhu <- function(level = c("County", "State", "Global"), start_date) {
  level <- match.arg(level)
  jhu_url <- paste0("https://raw.githubusercontent.com/lin-lab/COVID-data-cleaning/master/jhu_data/cleaned_data/JHU_COVID-19_", level, ".csv")

  dat <- data.table::fread(jhu_url)
  dat[, date := as.Date(date, format = "%Y-%m-%d")]
  ret <- dat[date >= start_date, ]
  return(ret)
}

# Weighted shift to "correct" for random lags
random_lag <- function(y, lags, weights) {
  stopifnot(length(weights) == length(lags))
  stopifnot(min(weights) >= 0)
  out <- rep(0, length(y))

  weights <- weights / sum(weights)

  for (i in seq_along(weights)) {
    if (weights[i] > 0) {
      out <- out + weights[i] * data.table::shift(y, n = -lags[i], fill = NA)
    }
  }
  return(out)
}

# Weighted shift to "correct" for random Gamma-distributed lags
gamma_lag <- function(y, max_days = 14, shape = 5, rate = 1) {
  l <- 0:max_days
  w <- pgamma(l + 1, shape, rate) - pgamma(l, shape, rate)
  w <- w / sum(w)
  return(random_lag(y, l, w))
}


# Estimate Rt
fit_poisson <- function(date, new_counts = NULL,
                        population = NULL, days_per_knot = 30,
                        min_positive = 50, adjust_weekday = TRUE,
                        family = quasipoisson(), SI_mean = 5.2, SI_sd = 5.5,
                        SI_NTS = 30, lagged_counts = TRUE, min_days = 14,
                        conf_level = 0.95) {

  stopifnot(0 < conf_level && conf_level < 1)
  stopifnot(SI_mean > 0 && SI_sd > 0)
  stopifnot(SI_NTS > 0)

  if (min(new_counts) < 0) {
    cat(sprintf("Warning: Set %d/%d negative case counts to 0.\n",
                sum(new_counts < 0), length(new_counts)))
    new_counts[new_counts < 0] <- 0
  }

  full_data <- data.table::data.table(new_counts = new_counts, date = date)
  data.table::setorder(full_data, date)

  if (lagged_counts) {
    full_data[, new_counts := gamma_lag(new_counts),]
  }

  # set offset variable: either Lambda from EpiEstim or log(population)
  if (is.null(population)) {
    full_data[, offset_var := log1p(getLambda(
      x = new_counts,
      t = as.integer(date),
      mu = SI_mean, sd = SI_sd, NTS = SI_NTS))]
    full_data[, offset_orig := exp(offset_var) - 1]
  } else {
    full_data[, `:=`(offset_var = log(population),
                     offset_orig = population)]
  }

  data <- full_data[cumsum(new_counts) >= min_positive, ]

  # data.frame to return in case of error
  out_error <- data.table::data.table(date = full_data$date,
                                      outcome_hat = NA, ci_lower = NA,
                                      ci_upper = NA, outcome_obs = NA,
                                      offset_var = NA)

  if (nrow(data) <= min_days) {
    return(out_error)
  }

  weekday_matrix <- function(dates) {
    model.matrix(~weekdays(dates))[, -1]
  }

  seq_days <- function(x, delta, min_interval = 0.5) {
    range_x <- max(x) - min(x)
    if (range_x < delta / 2) {
      return(mean(x))
    } else if (range_x < delta) {
      delta <- delta / 2
    }
    out <- seq(min(x) + delta, max(x), delta)
    if (max(x) - max(out) < delta * min_interval) {
      # take out the last knot if we go over
      out <- head(out, -1)
    }
    return(out)
  }

  if (adjust_weekday) {
    model_formula <- new_counts ~ offset(offset_var) +
      bs(date, knots = seq_days(date, days_per_knot)) + weekday_matrix(date)
  } else {
    model_formula <- new_counts ~ offset(offset_var) +
      bs(date, knots = seq_days(date, days_per_knot))
  }

  fit <- tryCatch({
    stats::glm(model_formula, data = data, family = family)
    },
    warning = function(w) {
      # usually will get a warning if GLM didn't converge
      message("Handled warning: ", conditionMessage(w))
      NULL
    })

  if (is.null(fit)) {
    return(out_error)
  }

  # redefine weekday_matrix so when we predict, we can predict using the avg of
  # the weekdays
  weekday_matrix <- function(dates) {
    out <- model.matrix(~weekdays(dates))[, -1]
    out[] <- 1 / 7
    out
  }

  pred <- stats::predict.glm(fit, newdata = data, type = "link", se.fit = TRUE)


  zstar <- qnorm(1 - (1 - conf_level) / 2)
  out <- data.table(date = data$date,
                    outcome_hat = exp(pred$fit) / data$offset_orig,
                    ci_lower = exp(pred$fit - zstar * pred$se.fit) /
                      data$offset_orig,
                    ci_upper = exp(pred$fit + zstar * pred$se.fit) /
                      data$offset_orig,
                    outcome_obs = data$new_counts,
                    offset_var = data$offset_var)
  return(out)
}


fit_poisson_by_unit <- function(data, days_per_knot = 30, min_positive = 50,
                                adjust_weekday = TRUE, family = quasipoisson(),
                                SI_mean = 5.2, SI_sd = 5.5, SI_NTS = 30,
                                lagged_counts = TRUE, min_days = 14) {

  uniq_uids <- unique(data$UID)
  n_uids <- length(uniq_uids)
  fit_poisson_default <-
    purrr::partial(fit_poisson, days_per_knot = days_per_knot,
                   min_positive = min_positive, adjust_weekday = adjust_weekday,
                   family = family, SI_mean = SI_mean, SI_sd = SI_sd,
                   SI_NTS = SI_NTS, lagged_counts = lagged_counts,
                   min_days = min_days)
  data.table::setkey(data, UID, date)

  fit_df_lst <- list()
  for (ii in 1:n_uids) {
    message(sprintf("Working on location %d of %d", ii, n_uids))
    cur_uid <- uniq_uids[[ii]]
    data_subset <- data[UID == cur_uid]
    rt_fit <- fit_poisson_default(date = data_subset$date,
                                  new_counts = data_subset$positiveIncrease)[,
              .(date, rt = outcome_hat, rt_lower = ci_lower,
                rt_upper = ci_upper, cases_lag = outcome_obs)]

    case_fit <- fit_poisson_default(date = data_subset$date,
                                    new_counts = data_subset$positiveIncrease,
                                    population = data_subset$population)[,
              .(date, case_rate = outcome_hat,
                case_lower = ci_lower, case_upper = ci_upper)]
    death_fit <- fit_poisson_default(date = data_subset$date,
                                     new_counts = data_subset$deathIncrease,
                                     population = data_subset$population)[,
              .(date, death_rate = outcome_hat, death_lower = ci_lower,
                death_upper = ci_upper)]
    tmp <- merge(rt_fit, case_fit, by = "date", all = TRUE)
    fit_df <- merge(tmp, death_fit, by = "date", all = TRUE)
    fit_df$UID <- cur_uid
    fit_df_lst[[ii]] <- fit_df
  }

  all_fitted <- do.call(rbind, fit_df_lst)
  all_fitted[, date_str := format(date, "%Y-%m-%d")]
  all_fitted[, date := NULL]
  data[, date_str := format(date, "%Y-%m-%d")]
  data.table::setkey(data, UID, date_str)
  data.table::setkey(all_fitted, UID, date_str)

  ret <- data[all_fitted, nomatch = NULL]
  ret[, date_str := NULL]
  return(ret)
}

start_date <- as.Date("2020-03-01", format = "%Y-%m-%d")
jhu_counties <- read_jhu("County", start_date = start_date)
jhu_states <- read_jhu("State", start_date = start_date)
jhu_global <- read_jhu("Global", start_date = start_date)

system.time({
  state_fit <- fit_poisson_by_unit(jhu_states)
  county_fit <- fit_poisson_by_unit(jhu_counties)
  global_fit <- fit_poisson_by_unit(jhu_global)
})

fwrite(state_fit, file = "jhu_state_rt_case_death_rate.csv")
fwrite(county_fit, file = "jhu_county_rt_case_death_rate.csv")
fwrite(global_fit, file = "jhu_global_rt_case_death_rate.csv")
