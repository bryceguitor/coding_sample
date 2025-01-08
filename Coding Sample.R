rm(list=ls())

# load libraries
library(ggplot2)
library(dplyr)
library(lubridate)
library(readr)
library(patchwork)
library(tseries)
library(tibble)
library(forecast)
library(FinTS)
library(rugarch)

# import the data
data <- read_csv("/Users/bryceguitor/Desktop/Carleton University MA ECON/Fall 2024/5029/bitcoin_data1.csv")

# convert DATE column to date type
data$DATE <- as.Date(data$DATE, format = "%Y-%m-%d")

# calculate daily log returns
data <- data %>%
  mutate(BTC_log_return = c(NA, diff(log(PRICE))))

# aggregate daily data into weekly data
data_weekly <- data %>%
  mutate(Week = floor_date(DATE, unit = "week")) %>%
  group_by(Week) %>%
  summarize(
    PRICE = last(PRICE), 
    BTC_log_return = sum(BTC_log_return, na.rm = TRUE)
  ) %>%
  ungroup()

# shared y-axis ranges for price and log returns
price_ylim <- range(c(data$PRICE, data_weekly$PRICE), na.rm = TRUE)  # for price
log_return_ylim <- range(c(data$BTC_log_return, data_weekly$BTC_log_return), na.rm = TRUE)  # for log returns

# function to create a plot
create_plot <- function(data, x_var, y_var, title, x_label, y_label, y_limits) {
  ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_line(color = "black") +
    ggtitle(title) +
    xlab(x_label) +
    ylab(y_label) +
    ylim(y_limits) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
}

# daily plots
btc_daily_price_plot <- create_plot(
  data = data, 
  x_var = "DATE", 
  y_var = "PRICE", 
  title = "Daily Bitcoin Price", 
  x_label = "Date", 
  y_label = "Price", 
  y_limits = price_ylim
)

btc_daily_log_return_plot <- create_plot(
  data = data, 
  x_var = "DATE", 
  y_var = "BTC_log_return", 
  title = "Daily Bitcoin Log Return", 
  x_label = "Date", 
  y_label = "Log Return", 
  y_limits = log_return_ylim
)

# weekly plots
btc_weekly_price_plot <- create_plot(
  data = data_weekly, 
  x_var = "Week", 
  y_var = "PRICE", 
  title = "Weekly Bitcoin Price", 
  x_label = "Date", 
  y_label = "Price", 
  y_limits = price_ylim
)

btc_weekly_log_return_plot <- create_plot(
  data = data_weekly, 
  x_var = "Week", 
  y_var = "BTC_log_return", 
  title = "Weekly Bitcoin Log Return", 
  x_label = "Date", 
  y_label = "Log Return", 
  y_limits = log_return_ylim
)

# combine daily and weekly plots
combined_daily_plot <- btc_daily_price_plot / btc_daily_log_return_plot
combined_weekly_plot <- btc_weekly_price_plot / btc_weekly_log_return_plot
print(combined_daily_plot)
print(combined_weekly_plot)


# function to generate summary statistics
generate_summary <- function(data, column_name, label) {
  cat(paste0("\nBase R Summary for Bitcoin ", label, " Log Returns:\n"))
  base_summary <- summary(data[[column_name]])
  print(base_summary)
  
  additional_summary <- data %>%
    summarize(
      Min = min(.data[[column_name]], na.rm = TRUE),
      Mean = mean(.data[[column_name]], na.rm = TRUE),
      Median = median(.data[[column_name]], na.rm = TRUE),
      Max = max(.data[[column_name]], na.rm = TRUE),
      SD = sd(.data[[column_name]], na.rm = TRUE)
    )
  
  cat(paste0("\nAdditional Statistics for Bitcoin ", label, " Log Returns:\n"))
  print(additional_summary)
}

# generate summaries
generate_summary(data, "BTC_log_return", "Daily")
generate_summary(data_weekly, "BTC_log_return", "Weekly")


# ADF test on daily and weekly log returns
adf_daily <- adf.test(na.omit(data$BTC_log_return))
adf_weekly <- adf.test(na.omit(data_weekly$BTC_log_return))  

# format ADF statistic and p-value 
statistic_daily <- formatC(adf_daily$statistic, format = "f", digits = 6)
statistic_weekly <- formatC(adf_weekly$statistic, format = "f", digits = 6)
p_value_daily <- formatC(adf_daily$p.value, format = "e", digits = 6)
p_value_weekly <- formatC(adf_weekly$p.value, format = "e", digits = 6)

# create tibble with ADF results
adf_results <- tibble(
  Series = c("Daily Log Returns", "Weekly Log Returns"),
  Statistic = c(statistic_daily, statistic_weekly),
  P_Value = c(p_value_daily, p_value_weekly),
)

# print results
cat("\nADF Test Results for Bitcoin Log Returns (High Precision):\n")
print(adf_results)


# find optimal lags for ARMA daily
cat("Optimal ARMA order for Bitcoin Daily Log Returns:\n")
btc_daily_arma <- auto.arima(
  na.omit(data$BTC_log_return), 
  max.p = 5, max.q = 5, ic = "aic", seasonal = FALSE 
)
print(btc_daily_arma)

# find optimal lags for ARMA weekly
cat("\nOptimal ARMA order for Bitcoin Weekly Log Returns:\n")
btc_weekly_arma <- auto.arima(
  na.omit(data_weekly$BTC_log_return), 
  max.p = 5, max.q = 5, ic = "aic", seasonal = FALSE
)
print(btc_weekly_arma)


# test for ARCH effects in daily
cat("\nARCH Test for Bitcoin Daily Log Returns:\n")
arch_test_daily <- ArchTest(na.omit(data$BTC_log_return), lags = 5) 
print(arch_test_daily)

# test for ARCH effects in weekly
cat("\nARCH Test for Bitcoin Weekly Log Returns:\n")
arch_test_weekly <- ArchTest(na.omit(data_weekly$BTC_log_return), lags = 5)
print(arch_test_weekly)


# function to fit a GARCH model
fit_garch_model <- function(data, column_name, model, garch_order = c(1, 1), arma_order = c(0, 0), 
                            distribution = "std", submodel = NULL) {
  spec <- ugarchspec(
    variance.model = list(model = model, garchOrder = garch_order, submodel = submodel),
    mean.model = list(armaOrder = arma_order, include.mean = TRUE),
    distribution.model = distribution
  )
  fit <- ugarchfit(spec = spec, data = na.omit(data[[column_name]]))
  return(fit)
}

# fit models
daily_garch_fit <- fit_garch_model(data, "BTC_log_return", model = "sGARCH")
daily_egarch_fit <- fit_garch_model(data, "BTC_log_return", model = "eGARCH")
daily_gjr_garch_fit <- fit_garch_model(data, "BTC_log_return", model = "fGARCH", submodel = "GJRGARCH")

# conditional variances
data <- data %>%
  mutate(
    Daily_GARCH_Volatility = c(rep(NA, length(data$BTC_log_return) - length(sigma(daily_garch_fit))), sigma(daily_garch_fit)),
    Daily_EGARCH_Volatility = c(rep(NA, length(data$BTC_log_return) - length(sigma(daily_egarch_fit))), sigma(daily_egarch_fit)),
    Daily_GJR_GARCH_Volatility = c(rep(NA, length(data$BTC_log_return) - length(sigma(daily_gjr_garch_fit))), sigma(daily_gjr_garch_fit))
  )

# daily plot
daily_volatility_plot <- ggplot(data, aes(x = DATE)) +
  geom_line(aes(y = Daily_GARCH_Volatility, color = "GARCH")) +
  geom_line(aes(y = Daily_EGARCH_Volatility, color = "EGARCH")) +
  geom_line(aes(y = Daily_GJR_GARCH_Volatility, color = "GJR-GARCH")) +
  ggtitle("Daily Volatility Estimates for Bitcoin") +
  xlab("Date") +
  ylab("Volatility") +
  scale_color_manual(values = c(
    "GARCH" = "black",
    "EGARCH" = "#D3D3D3",
    "GJR-GARCH" = "#999999"
  )) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.95, 0.95), 
    legend.justification = c("right", "top"),  
    legend.background = element_rect(fill = "white", color = "white") 
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )


# fit models for weekly data
weekly_garch_fit <- fit_garch_model(data_weekly, "BTC_log_return", model = "sGARCH")
weekly_egarch_fit <- fit_garch_model(data_weekly, "BTC_log_return", model = "eGARCH")
weekly_gjr_garch_fit <- fit_garch_model(data_weekly, "BTC_log_return", model = "fGARCH", submodel = "GJRGARCH")

# conditional variances
data_weekly <- data_weekly %>%
  mutate(
    Weekly_GARCH_Volatility = c(rep(NA, length(data_weekly$BTC_log_return) - length(sigma(weekly_garch_fit))), sigma(weekly_garch_fit)),
    Weekly_EGARCH_Volatility = c(rep(NA, length(data_weekly$BTC_log_return) - length(sigma(weekly_egarch_fit))), sigma(weekly_egarch_fit)),
    Weekly_GJR_GARCH_Volatility = c(rep(NA, length(data_weekly$BTC_log_return) - length(sigma(weekly_gjr_garch_fit))), sigma(weekly_gjr_garch_fit))
  )

# weekly plot
weekly_volatility_plot <- ggplot(data_weekly, aes(x = Week)) +
  geom_line(aes(y = Weekly_GARCH_Volatility, color = "GARCH")) +
  geom_line(aes(y = Weekly_EGARCH_Volatility, color = "EGARCH")) +
  geom_line(aes(y = Weekly_GJR_GARCH_Volatility, color = "GJR-GARCH")) +
  ggtitle("Weekly Volatility Estimates for Bitcoin") +
  xlab("Date") +
  ylab("Volatility") +
  scale_color_manual(values = c(
    "GARCH" = "black",
    "EGARCH" = "#D3D3D3",
    "GJR-GARCH" = "#999999"
  )) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.95, 0.95),  
    legend.justification = c("right", "top"),  
    legend.background = element_rect(fill = "white", color = "white") 
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

# display plots
print(daily_volatility_plot)
print(weekly_volatility_plot)


# function to split data into training and testing sets
split_data <- function(data, train_ratio = 0.8) {
  split_index <- floor(train_ratio * nrow(data))
  list(
    training_set = data[1:split_index, ],
    testing_set = data[(split_index + 1):nrow(data), ]
  )
}

# function to calculate realized volatility proxy
calculate_proxy <- function(data, column_name, mean_value, proxy_name) {
  data %>%
    mutate(!!proxy_name := (get(column_name) - mean_value)^2) %>%
    filter(!is.na(!!sym(proxy_name)))
}

# function to calculate MAE
calculate_mae <- function(test_data, proxy_column, model_column) {
  mean(abs(test_data[[proxy_column]] - test_data[[model_column]]), na.rm = TRUE)
}

# function to calculate RMSE
calculate_rmse <- function(test_data, proxy_column, model_column) {
  sqrt(mean((test_data[[proxy_column]] - test_data[[model_column]])^2, na.rm = TRUE))
}

# function to display model coefficients
display_model <- function(model, label) {
  cat(paste0("\n", label, " Model Coefficients:\n"))
  show(model)
}

# process daily data
daily_split <- split_data(data)
mean_log_return <- mean(daily_split$training_set$BTC_log_return, na.rm = TRUE)
testing_set <- calculate_proxy(daily_split$testing_set, "BTC_log_return", mean_log_return, "Realized_Volatility_Proxyd")


# base R summary and additional statistics for daily realized volatility proxy
daily_vol_summary <- summary(testing_set$Realized_Volatility_Proxyd)
generate_summary(testing_set, "Realized_Volatility_Proxyd", "Daily Realized Volatility Proxy")


# define models and corresponding volatility columns
models <- c("GARCH" = "Daily_GARCH_Volatility", 
            "EGARCH" = "Daily_EGARCH_Volatility", 
            "GJR-GARCH" = "Daily_GJR_GARCH_Volatility")

# loop through models for MAE and RMSE (daily)
cat("\nDaily Data Metrics:\n")
for (model in names(models)) {
  mae <- calculate_mae(testing_set, "Realized_Volatility_Proxyd", models[model])
  rmse <- calculate_rmse(testing_set, "Realized_Volatility_Proxyd", models[model])
  cat("MAE for Daily", model, "Model:", mae, "\n")
  cat("RMSE for Daily", model, "Model:", rmse, "\n")
}

# display model summaries (daily)
daily_fits <- list("GARCH" = daily_garch_fit, 
                   "EGARCH" = daily_egarch_fit, 
                   "GJR-GARCH" = daily_gjr_garch_fit)

for (model in names(daily_fits)) {
  display_model(daily_fits[[model]], paste("Daily", model))
}

# process weekly data
weekly_split <- split_data(data_weekly)
mean_log_return_weekly <- mean(weekly_split$training_set$BTC_log_return, na.rm = TRUE)
testing_set_weekly <- calculate_proxy(weekly_split$testing_set, "BTC_log_return", mean_log_return_weekly, "Realized_Volatility_Proxyw")


# base R summary and additional statistics for daily realized volatility proxy
weekly_vol_summary <- summary(testing_set_weekly$Realized_Volatility_Proxyw)
generate_summary(testing_set_weekly, "Realized_Volatility_Proxyw", "Weekly Realized Volatility Proxy")


# define models and corresponding volatility columns
weekly_models <- c("GARCH" = "Weekly_GARCH_Volatility", 
                   "EGARCH" = "Weekly_EGARCH_Volatility", 
                   "GJR-GARCH" = "Weekly_GJR_GARCH_Volatility")

# loop through models for MAE and RMSE (weekly)
cat("\nWeekly Data Metrics:\n")
for (model in names(weekly_models)) {
  mae <- calculate_mae(testing_set_weekly, "Realized_Volatility_Proxyw", weekly_models[model])
  rmse <- calculate_rmse(testing_set_weekly, "Realized_Volatility_Proxyw", weekly_models[model])
  cat("MAE for Weekly", model, "Model:", mae, "\n")
  cat("RMSE for Weekly", model, "Model:", rmse, "\n")
}

# display model summaries (weekly)
weekly_fits <- list("GARCH" = weekly_garch_fit, 
                    "EGARCH" = weekly_egarch_fit, 
                    "GJR-GARCH" = weekly_gjr_garch_fit)

for (model in names(weekly_fits)) {
  display_model(weekly_fits[[model]], paste("Weekly", model))
}





