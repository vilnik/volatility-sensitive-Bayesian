rm(list=ls())

library("quantmod")
library("xts")
library("expm")
library("plyr")
library("rmgarch")
library("MASS")

source("calculations.R")
source("portfolio_VaRs.R")

# Parameters controlling what to run
COMPARE_PORTFOLIOS_BY_SIMULATION = TRUE
COMPARE_PORTFOLIOS_EMPIRICALLY =  TRUE

# Parameters
alpha = 0.975 # VaR level
n = 500 # Rolling window size
k = 20 # Number of stocks
year = 2019 # Empirical year
simulation_data_type = "SIMULATION_NORMAL" # Simulation scenario "SIMULATION_NORMAL", "SIMULATION_NORMAL_SWITCHING" or "SIMULATION_MGARCH"
simulation_size = 1000

# Hyperparameters to control volatility sensitivity method
conjugate_n_r = 3
conjugate_h = 4
conjugate_l = -4
conjugate_hyperparameters = list("conjugate_n_r" = conjugate_n_r, "conjugate_h" = conjugate_h, "conjugate_l" = conjugate_l)

# Input stock prices as a list of xts objects e.g., list("AAPL" = <xts object with daily market caps (no specific column name required)>, ...)
all_stocks_data = readRDS(file = 'all_stocks_data.rds')

# ======================================================================================
# Portfolio VaR evaluation
# ======================================================================================

Portfolio_VaR_evaluation <- function(data_type, all_stocks_data, n, k, alpha, conjugate_hyperparameters, year, simulation_size)
{
  if (data_type == "SIMULATION_NORMAL" | data_type == "SIMULATION_NORMAL_SWITCHING" | data_type == "SIMULATION_MGARCH")
  {
    stocks_data_for_comparison = vector(mode = 'list')
    stocks = list()
    stocks_adjusted = list()
    number_of_dates = n + simulation_size + 1
    
    for (i in 1: k)
    {
      stock <- paste("STOCK", i, sep = "")
      stock_adjusted <- paste(stock, "Adjusted", sep = ".")
      stocks = c(stocks, stock)
      stocks_adjusted = c(stocks_adjusted, stock_adjusted)
    }
    
    k_stocks_data = all_stocks_data[sample.int(length(names(all_stocks_data)), k)]
    k_stocks_names = names(k_stocks_data)
    date_indexes = index(k_stocks_data[[k_stocks_names[1]]])
    date_indexes_n_adjusted = date_indexes[(n + 1): length(date_indexes)]
    date_index_for_parameter_estimation_in_simulation = sample(date_indexes_n_adjusted, 1)
    
    if (data_type == "SIMULATION_NORMAL")
    {
      simulation_parameters = Conventional_parameters(k_stocks_data, date_index_for_parameter_estimation_in_simulation, n, x_bar = NULL, S = NULL)
      k_log_returns = rmvnorm(n = number_of_dates, mean = simulation_parameters$conventional_mu, sigma = simulation_parameters$conventional_sigma)
    }
    
    if (data_type == "SIMULATION_NORMAL_SWITCHING")
    {
      simulation_parameters = Conventional_parameters(k_stocks_data, date_index_for_parameter_estimation_in_simulation, n, x_bar = NULL, S = NULL)
      i = 0
      k_log_returns = c()
      while (i < number_of_dates)
      {
        u = runif(1, min = 0, max = 1)
        days = min(sample(2:5, 1), number_of_dates - i)
        if (u <= 0.05)
        {
          volatility_factor = replicate(k, runif(1, 0.5, 0.75))
        } else if (u <= 0.95)
        { 
          volatility_factor = replicate(k, 1)
        } else
        {
          volatility_factor = replicate(k, runif(1, 1.5, 3))
        }
        k_log_returns_window = rmvnorm(n = days, mean = simulation_parameters$conventional_mu, sigma = (diag(volatility_factor) %*% simulation_parameters$conventional_sigma %*% diag(volatility_factor)))
        k_log_returns = rbind(k_log_returns, k_log_returns_window)
        i = i + days
      }
    }
    if (data_type == "SIMULATION_MGARCH")
    {
      k_historical_log_returns = All_log_returns(k_stocks_data)
      k_historical_log_returns = k_historical_log_returns[-1, ]
      dates = index(k_historical_log_returns)
      date_end = which(dates == date_index_for_parameter_estimation_in_simulation)
      k_historical_log_returns = k_historical_log_returns[(date_end - n + 1):date_end, ]
      
      uspec1 = ugarchspec(mean.model = list(armaOrder = c(0, 0)), variance.model = list(garchOrder = c(1, 1), model = 'sGARCH'), distribution.model = 'norm')
      spec1 = dccspec(uspec = multispec(replicate(k, uspec1)), dccOrder = c(1, 1), distribution = 'mvnorm')
      fit_dcc = dccfit(spec1, data = k_historical_log_returns, fit.control = list(scale = TRUE), solver = "gosolnp")
      
      dcc_sim = dccsim(fit_dcc, n.sim = number_of_dates, rseed = 100)
      k_log_returns_dcc = dcc_sim@msim$simX[[1]]
      
      k_log_returns = k_log_returns_dcc
    }
    
    for (i in c(1:k))
    {
      stock_adjusted_prices <- c()
      # All stock prices are initially 100
      stock_adjusted_prices[1] <- 100
      for (j in 2:length(k_log_returns[, i]))
      {
        stock_adjusted_prices[j] <- stock_adjusted_prices[j - 1]*exp(k_log_returns[j, i])
      }
      sim_stock_data <- xts(rev(stock_adjusted_prices), order.by = Sys.Date() - 1:number_of_dates)
      names(sim_stock_data) <- stocks_adjusted[i]
      stocks_data_for_comparison[[stocks[[i]]]] <- sim_stock_data
    }
  } else if (data_type == "EMPIRICAL")
  {
    all_stocks_data_for_current_year = list()
    start_date = as.Date(paste0(year - ceiling(n / 250), "-01-01"))
    end_date = as.Date(paste0(year, "-12-31"))
    for (stock_i in c(1:length(names(all_stocks_data))))
    {
      all_stocks_data_for_current_year[[names(all_stocks_data)[stock_i]]] = all_stocks_data[[names(all_stocks_data)[stock_i]]][start_date <= index(all_stocks_data[[names(all_stocks_data)[stock_i]]]) & index(all_stocks_data[[names(all_stocks_data)[stock_i]]]) <= end_date]
    }
    
    stocks_data_for_comparison = all_stocks_data_for_current_year[sample.int(length(names(all_stocks_data_for_current_year)), k)]
    k_stocks_names = names(stocks_data_for_comparison)
  }
  
  # Calculate VaR
  portfolio_daily_VaRs_result = Calculate_portfolio_daily_VaRs(stocks_data_for_comparison, n, k, alpha, conjugate_hyperparameters, year, simulation_size)

  return(list("portfolio_daily_VaRs_result" = portfolio_daily_VaRs_result, "k_stocks_names" = k_stocks_names, "stocks_data_for_comparison" = stocks_data_for_comparison))
}

# ======================================================================================
# Simulation comparison
# ======================================================================================

if (COMPARE_PORTFOLIOS_BY_SIMULATION)
{
  simulation_out = Portfolio_VaR_evaluation(simulation_data_type, all_stocks_data, n, k, alpha, conjugate_hyperparameters, year = NA, simulation_size)
}

# ======================================================================================
# Empirical comparison
# ======================================================================================

if (COMPARE_PORTFOLIOS_EMPIRICALLY)
{
  empirical_out = Portfolio_VaR_evaluation("EMPIRICAL", all_stocks_data, n, k, alpha, conjugate_hyperparameters, year, simulation_size = NA)
}
