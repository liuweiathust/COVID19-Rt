library(dplyr)
library(magrittr)
library(data.table)

source("../preprocess_data/preprocess_jhu.R")
source("../helpers/helper_functions.R")

today_date <- format(Sys.Date(), "%Y-%m-%d")
jhu_counties <- load_jhu(level='County', start_date = '2020-03-01', end_date = today_date) %>%
  data.table()
jhu_states <- load_jhu(level='State',  start_date = '2020-03-01', end_date = today_date) %>%
  data.table()
jhu_global <- load_jhu(level='Global', start_date = '2020-03-01', end_date = today_date) %>%
  data.table()
