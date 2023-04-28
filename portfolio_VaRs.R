# ======================================================================================
# Calculate portfolio daily VaRs
# ======================================================================================

# Portfolio daily VaRs
Calculate_portfolio_daily_VaRs <- function(stocks_data, n, k, alpha, conjugate_hyperparameters, year, simulation_size)
{
  # Stock names
  stocks = names(stocks_data)
  
  # Define vectors and variables used to store results for portfolio
  empirical_portfolio_expected_returns = c()
  conjugate_portfolio_expected_returns = c()
  conventional_portfolio_expected_returns = c()
  dcc_portfolio_expected_returns = c()
  
  portfolio_actual_returns = c()
  portfolio_weights = rep(1 / k, k)
  
  # VaRs for different methods
  empirical_portfolio_alpha_VaRs = c()
  conjugate_portfolio_alpha_VaRs = c()
  conventional_portfolio_alpha_VaRs = c()
  dcc_portfolio_alpha_VaRs = c()

  all_dates = index(stocks_data[[names(stocks_data)[1]]])
  
  if (!is.na(year))
  {
    start_date = as.Date(paste0(year, "-01-01"))
    end_date = as.Date(paste0(year, "-12-31"))
    
    start_date_i = which(all_dates >= start_date)[1] - 1
    end_date_i = which(all_dates <= end_date)[length(which(all_dates <= end_date))] - 1
    number_of_dates = end_date_i - start_date_i + 1
  } else
  {
    start_date_i = n + 1
    end_date_i = start_date_i + simulation_size - 1
    number_of_dates = end_date_i - start_date_i + 1
  }

  dcc_fit = Dcc_fit(stocks_data, number_of_dates)
  dcc_forecast = rcov(dccforecast(dcc_fit, n.ahead = 1, n.roll = number_of_dates))
  all_log_returns = All_log_returns(stocks_data)
  
  # Loop over all dates and calculate portfolio values
  return_dates = c()
  for (date_i in c(start_date_i: end_date_i))
  {
    date = all_dates[date_i]
    print(paste0("Date ", date))
    date_plus_one = all_dates[date_i + 1]
    return_dates = c(return_dates, date_plus_one)
    S = S_matrix(date, stocks_data, n)
    x_bar = X_bar(date, stocks_data, n)
    
    # Parameters
    empirical_parameters = Bayesian_parameters(stocks_data, date, n, k, "CONJUGATE_EMPIRICAL", x_bar = x_bar, S = S)
    conjugate_parameters = Bayesian_parameters(stocks_data, date, n, k, "CONJUGATE_VS", x_bar = x_bar, S = S, conjugate_hyperparameters = conjugate_hyperparameters, weights = portfolio_weights)
    conventional_parameters = Conventional_parameters(stocks_data, date, n, x_bar = x_bar, S = S)

    dcc_forecast_cov_date_i = date_i - start_date_i + 1
    if (!names(dcc_forecast[dcc_forecast_cov_date_i]) == date)
    {
      stop("Error! DCC forecast date does not match actual date!")
    }
      
    dcc_parameters = Dcc_parameters(stocks_data, dcc_forecast[[dcc_forecast_cov_date_i]][,,1], date, x_bar = x_bar)
    
    # Rho alpha
    empirical_alpha_VaR_rho_alpha = Bayesian_rho_alpha("VaR", alpha, empirical_parameters$bayesian_d_kn)
    conjugate_alpha_VaR_rho_alpha = Bayesian_rho_alpha("VaR", alpha, conjugate_parameters$bayesian_d_kn)
    conventional_alpha_VaR_rho_alpha = Conventional_rho_alpha("VaR", alpha)
    
    # Portfolio VaR
    empirical_portfolio_alpha_VaR = Bayesian_Q(empirical_parameters, empirical_alpha_VaR_rho_alpha, portfolio_weights)
    empirical_portfolio_alpha_VaRs = c(empirical_portfolio_alpha_VaRs, empirical_portfolio_alpha_VaR)
    
    conjugate_portfolio_alpha_VaR = Bayesian_Q(conjugate_parameters, conjugate_alpha_VaR_rho_alpha, portfolio_weights)
    conjugate_portfolio_alpha_VaRs = c(conjugate_portfolio_alpha_VaRs, conjugate_portfolio_alpha_VaR)
    
    conventional_portfolio_alpha_VaR = Conventional_Q(conventional_parameters, conventional_alpha_VaR_rho_alpha, portfolio_weights)
    conventional_portfolio_alpha_VaRs = c(conventional_portfolio_alpha_VaRs, conventional_portfolio_alpha_VaR)
    
    dcc_portfolio_alpha_VaR = Dcc_Q(dcc_parameters, conventional_alpha_VaR_rho_alpha, portfolio_weights)
    dcc_portfolio_alpha_VaRs = c(dcc_portfolio_alpha_VaRs, dcc_portfolio_alpha_VaR)
    
    # Portfolio expected returns
    empirical_portfolio_expected_return = t(portfolio_weights) %*% empirical_parameters$bayesian_x_bar
    empirical_portfolio_expected_returns = c(empirical_portfolio_expected_returns, empirical_portfolio_expected_return)
    
    conjugate_portfolio_expected_return = t(portfolio_weights) %*% conjugate_parameters$bayesian_x_bar
    conjugate_portfolio_expected_returns = c(conjugate_portfolio_expected_returns, conjugate_portfolio_expected_return)
    
    conventional_portfolio_expected_return = t(portfolio_weights) %*% conventional_parameters$conventional_mu
    conventional_portfolio_expected_returns = c(conventional_portfolio_expected_returns, conventional_portfolio_expected_return)
    
    dcc_portfolio_expected_return = t(portfolio_weights) %*% dcc_parameters$dcc_mu
    dcc_portfolio_expected_returns = c(dcc_portfolio_expected_returns, dcc_portfolio_expected_return)
    
    # Portfolio actual returns
    portfolio_actual_return = t(portfolio_weights) %*% as.numeric(all_log_returns[date_plus_one])
    portfolio_actual_returns = c(portfolio_actual_returns, portfolio_actual_return)
  }
  
  return(list("stocks" = stocks, "return_dates" = return_dates,
              "empirical_portfolio_expected_returns" = empirical_portfolio_expected_returns, "conjugate_portfolio_expected_returns" = conjugate_portfolio_expected_returns, "conventional_portfolio_expected_returns" = conventional_portfolio_expected_returns, "dcc_portfolio_expected_returns" = dcc_portfolio_expected_returns,
              "portfolio_actual_returns" = portfolio_actual_returns, "portfolio_weights" = portfolio_weights,
              "empirical_portfolio_alpha_VaRs" = empirical_portfolio_alpha_VaRs, "conjugate_portfolio_alpha_VaRs" = conjugate_portfolio_alpha_VaRs, "conventional_portfolio_alpha_VaRs" = conventional_portfolio_alpha_VaRs, "dcc_portfolio_alpha_VaRs" = dcc_portfolio_alpha_VaRs))
}
