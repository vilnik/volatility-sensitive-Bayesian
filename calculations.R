# ======================================================================================
# Calculate log returns
# ======================================================================================

# Log returns for a single stock
Log_returns = function(stock, stock_data)
{
  column = paste(stock, "Adjusted", sep=".")
  t = diff(log(stock_data[,column]))
  names(t) <- paste(stock, "LogReturns", sep=".")
  tclass(t) <- "Date"
  return(t)
}

# Mean log returns for a single stock
Mean_log_returns = function(stock, log_returns, n)
{
  t = rollapply(log_returns, n, mean)
  names(t) <- paste(stock, "MeanLogReturns", sep=".")
  tclass(t) <- "Date"
  return(t)
}

# Compute all log returns
All_log_returns = function(stocks_data)
{
  all_log_returns = xts()
  tclass(all_log_returns) <- "Date"
  for (stock in names(stocks_data))
  {
    all_log_returns = merge(all_log_returns, Log_returns(stock, stocks_data[[stock]])) 
  }
  return(all_log_returns)
}

# Compute all mean log returns
All_mean_log_returns = function(stocks_data, n) 
{
  all_log_returns = All_log_returns(stocks_data)
  all_mean_log_returns = xts()
  tclass(all_mean_log_returns) <- "Date"
  for (stock in names(stocks_data))
  {
    all_mean_log_returns = merge(all_mean_log_returns, Mean_log_returns(stock, all_log_returns[, paste(stock, 'LogReturns', sep='.')], n)) 
  }
  return(all_mean_log_returns)
}

# ======================================================================================
# Calculate x_bar and S matrix
# ======================================================================================

# Calculate x_bar vector at some date
X_bar = function(date, stocks_data, n)
{
  matrix(All_mean_log_returns(stocks_data, n)[date,], ncol = 1)
}

# Calculate S matrix at some date
S_matrix = function(date, stocks_data, n, x_bar = NULL)
{
  all_log_returns = All_log_returns(stocks_data)
  if (is.null(x_bar))
  {
    x_bar = X_bar(date, stocks_data, n)
  }
  all_dates = index(all_log_returns)
  date_end_i = which(all_dates == date)
  S = matrix(0, dim(all_log_returns)[2], dim(all_log_returns)[2])
  for (date_i in ((date_end_i - n + 1):date_end_i))
  {
    S_addition = matrix(all_log_returns[all_dates[date_i],], ncol = 1) - x_bar
    S = S + S_addition %*% t(S_addition)
  }
  return(S)
}

# ======================================================================================
# Calculate parameters
# ======================================================================================

# Calculate Bayesian parameters
Bayesian_parameters = function(stocks_data, date, n, k, prior, x_bar = NULL, S = NULL, conjugate_hyperparameters = NULL, weights = NULL)
{
  if (is.null(S))
  {
    S = S_matrix(date, stocks_data, n)
  }
  if (is.null(x_bar))
  {
    x_bar = X_bar(date, stocks_data, n)
  }
  
  if (prior == "JEFFREYS")
  {
    bayesian_x_bar = x_bar
    bayesian_S = S
    bayesian_S_inv = solve(bayesian_S)
    bayesian_d_kn = n - k
    bayesian_r_kn = (n + 1) / (n * (n - k))
  }
  else if (prior == "CONJUGATE_EMPIRICAL")
  {
    d_0 = n
    # Reference 'Bayesian inference of the multi-period optimal portfolio for an exponential utility'
    S_0 = (d_0 - k - 1) / n * S
    r_0 = n
    m_0 = x_bar
    bayesian_x_bar = (n * x_bar + r_0 * m_0) / (n + r_0)
    bayesian_S = S + S_0 + n * r_0 * ((m_0 - bayesian_x_bar) %*% t(m_0 - bayesian_x_bar)) / (n + r_0)
    bayesian_S_inv = solve(bayesian_S)
    bayesian_d_kn = n + d_0 - 2 * k
    bayesian_r_kn = (n + r_0 + 1) / ((n + r_0) * (n + d_0 - 2 * k))
  }
  else if (prior == "CONJUGATE_VS")
  {
    n_recent = conjugate_hyperparameters$conjugate_n_r
    high_volatility_exponent = conjugate_hyperparameters$conjugate_h
    low_volatility_exponent = conjugate_hyperparameters$conjugate_l

    S_recent = S_matrix(date, stocks_data, n_recent, x_bar = x_bar)
    x_bar_recent = X_bar(date, stocks_data, n_recent)
      
    # Calculate recent and full sample estimates of sigma
    conventional_sigma_hat = 1 / (n - 1) * S
    conventional_sigma_hat_recent = 1 / n_recent * S_recent
    
    volatility_recent = diag(conventional_sigma_hat_recent)^(1/2)
    volatility_long_term = diag(conventional_sigma_hat)^(1/2)
    
    volatility_recent_long_term_ratio = volatility_recent / volatility_long_term
    volatility_recent_long_term_ratio_diag = diag(volatility_recent_long_term_ratio)
    
    covariance_matrix_recent = (volatility_recent_long_term_ratio_diag) %*% conventional_sigma_hat %*% (volatility_recent_long_term_ratio_diag)
    covariance_matrix_long_term = conventional_sigma_hat
    
    portfolio_variance_recent = as.numeric(t(weights) %*% covariance_matrix_recent %*% weights)
    portfolio_variance_long_term = as.numeric(t(weights) %*% covariance_matrix_long_term %*% weights)

    d_0 = max(k + 2, n * (max(1, portfolio_variance_recent / portfolio_variance_long_term))^high_volatility_exponent * (max(1, portfolio_variance_long_term / portfolio_variance_recent))^low_volatility_exponent)
    S_0 = ((d_0 - k - 1) * (n - 1) / n * covariance_matrix_recent)
    r_0 = n
    m_0 = x_bar
    bayesian_x_bar = (n * x_bar + r_0 * m_0) / (n + r_0)
    bayesian_S = S + S_0 + n * r_0 * ((m_0 - bayesian_x_bar) %*% t(m_0 - bayesian_x_bar)) / (n + r_0)
    
    if (rcond(bayesian_S) < .Machine$double.eps)
    {
      print("Sigma matrix from conjugate adaptive is near singular.")
      stop("Near singular matrix.") 
    }
    
    bayesian_S_inv = solve(bayesian_S)
    bayesian_d_kn = n + d_0 - 2 * k
    bayesian_r_kn = (n + r_0 + 1) / ((n + r_0) * (n + d_0 - 2 * k))
  }
  return(list("bayesian_x_bar" = bayesian_x_bar, "bayesian_S" = bayesian_S, "bayesian_S_inv" = bayesian_S_inv, "bayesian_d_kn" = bayesian_d_kn, "bayesian_r_kn" = bayesian_r_kn))
}

 
# Calculate conventional parameters
Conventional_parameters = function(stocks_data, date, n, x_bar = NULL, S = NULL)
{
  if (is.null(x_bar))
  {
    x_bar = X_bar(date, stocks_data, n)
  }

  if (is.null(S))
  {
    S = S_matrix(date, stocks_data, n)
  }
  
  conventional_sigma_hat = 1 / (n - 1) * S
  
  conventional_sigma_hat_inv = solve(conventional_sigma_hat)
  conventional_mu_hat = x_bar
  return(list("conventional_mu" = conventional_mu_hat, "conventional_sigma" = conventional_sigma_hat, "conventional_sigma_inv" = conventional_sigma_hat_inv))
}

# Calculate DCC-GARCH parameters
Dcc_fit <- function(stocks_data, out_of_sample)
{
  log_returns = All_log_returns(stocks_data)
  # Have to remove first element in all_log_returns since it is NA
  log_returns = log_returns[-1, ]
  k = length(names(log_returns))
  
  # DCC model
  uspec1 = ugarchspec(mean.model = list(armaOrder = c(0, 0)), variance.model = list(garchOrder = c(1, 1), model = 'sGARCH'), distribution.model = 'norm')
  spec1 = dccspec(uspec = multispec(replicate(k, uspec1)), dccOrder = c(1, 1), distribution = 'mvnorm')
  fit_dcc = dccfit(spec1, data = log_returns, out.sample = out_of_sample, solver = "gosolnp")
  return(fit_dcc)
}

# DCC-GARCH(1,1) parameters
Dcc_parameters <- function(stocks_data, dcc_forecast_H, date, x_bar = NULL)
{
  if (is.null(x_bar))
  {
    x_bar = X_bar(date, stocks_data, n)
  }
  k = length(names(stocks_data))
  dcc_mu_hat = x_bar
  dcc_sigma_hat = matrix(dcc_forecast_H, ncol = k, nrow = k, byrow = TRUE)
  dcc_sigma_hat_inv = solve(dcc_sigma_hat)
  return(list("dcc_mu" = dcc_mu_hat, "dcc_sigma" = dcc_sigma_hat, "dcc_sigma_inv" = dcc_sigma_hat_inv))
}

# ======================================================================================
# Calculate VaR and CVaR 
# ======================================================================================
  
# Calculate Q (VaR or CVaR) of a Bayesian portfolio
Bayesian_Q = function(bayesian_parameters, bayesian_rho_alpha, weights)
{
  bayesian_x_bar = bayesian_parameters$bayesian_x_bar
  bayesian_S = bayesian_parameters$bayesian_S
  bayesian_r_kn = bayesian_parameters$bayesian_r_kn
  return(as.numeric(-t(weights) %*% bayesian_x_bar + bayesian_rho_alpha * bayesian_r_kn^(1 / 2) * (t(weights) %*% bayesian_S %*% weights)^(1 / 2)))
}

# Calculate Q (VaR or CVaR) of a conventional portfolio
Conventional_Q = function(conventional_parameters, conventional_rho_alpha, weights)
{
  conventional_mu = conventional_parameters$conventional_mu
  conventional_sigma = conventional_parameters$conventional_sigma
  return(as.numeric(-t(weights) %*% conventional_mu + conventional_rho_alpha * (t(weights) %*% conventional_sigma %*% weights)^(1 / 2)))
}

# Calculate Q (VaR or CVaR) of a DCC-GARCH portfolio
Dcc_Q = function(dcc_parameters, dcc_rho_alpha, weights)
{
  dcc_mu = dcc_parameters$dcc_mu
  dcc_sigma = dcc_parameters$dcc_sigma
  return(as.numeric(-t(weights) %*% dcc_mu + dcc_rho_alpha * (t(weights) %*% dcc_sigma %*% weights)^(1 / 2)))
}
# ======================================================================================
# Calculate M matrix and s (M is also sometimes referred to as Q)
# ======================================================================================

# Calculate Bayesian M matrix
Bayesian_M_matrix = function(bayesian_parameters)
{
  bayesian_S_inv = bayesian_parameters$bayesian_S_inv
  ones_vector = matrix(1, dim(bayesian_S_inv)[1], 1)
  return(bayesian_S_inv - (bayesian_S_inv %*% ones_vector %*% t(ones_vector) %*% bayesian_S_inv) / as.numeric((t(ones_vector) %*% bayesian_S_inv %*% ones_vector)))
}

# Calculate conventional M matrix
Conventional_M_matrix = function(conventional_parameters)
{
  conventional_sigma_inv = conventional_parameters$conventional_sigma_inv
  ones_vector = matrix(1, dim(conventional_sigma_inv)[1], 1)
  return(conventional_sigma_inv - conventional_sigma_inv %*% ones_vector %*% t(ones_vector) %*% conventional_sigma_inv / as.numeric(t(ones_vector) %*% conventional_sigma_inv %*% ones_vector))
}

# Calculate DCC_GARCH M matrix
Dcc_M_matrix = function(dcc_parameters)
{
  dcc_sigma_inv = dcc_parameters$dcc_sigma_inv
  ones_vector = matrix(1, dim(dcc_sigma_inv)[1], 1)
  return(dcc_sigma_inv - dcc_sigma_inv %*% ones_vector %*% t(ones_vector) %*% dcc_sigma_inv / as.numeric(t(ones_vector) %*% dcc_sigma_inv %*% ones_vector))
}

# Calculate Bayesian s
Bayesian_s = function(bayesian_parameters, bayesian_M)
{
  bayesian_x_bar = bayesian_parameters$bayesian_x_bar
  return(as.numeric(t(bayesian_x_bar) %*% bayesian_M %*% bayesian_x_bar))
}

# Calculate conventional s
Conventional_s = function(conventional_parameters, conventional_M)
{
  conventional_mu = conventional_parameters$conventional_mu
  return(as.numeric(t(conventional_mu) %*% conventional_M %*% conventional_mu))
}

# Calculate DCC-GARCH s
Dcc_s = function(dcc_parameters, dcc_M)
{
  dcc_mu = dcc_parameters$dcc_mu
  return(as.numeric(t(dcc_mu) %*% dcc_M %*% dcc_mu))
}
# ======================================================================================
# Calculate rho_alpha
# ======================================================================================

# Calculate Bayesian rho_alpha
Bayesian_rho_alpha = function(VaRorCVaR, alpha, bayesian_d_kn)
{
  d_alpha = qt(alpha, bayesian_d_kn)
  if (VaRorCVaR == "VaR")
  {
    return(d_alpha)
  } else if (VaRorCVaR == "CVaR")
  {
    return(1 / (1 - alpha) * exp(lgamma((bayesian_d_kn + 1) / 2) - lgamma(bayesian_d_kn / 2)) / (pi * bayesian_d_kn)^(1 / 2) * bayesian_d_kn / (bayesian_d_kn - 1) * ((1 + d_alpha^2 / bayesian_d_kn)^(-(bayesian_d_kn - 1) / 2)))
  }
}

# Calculate conventional rho_alpha
Conventional_rho_alpha = function(VaRorCVaR, alpha)
{
  n_alpha = qnorm(alpha)
  if (VaRorCVaR == "VaR")
  {
    return(n_alpha)
  } else if (VaRorCVaR == "CVaR")
  {
    return(exp(-(n_alpha^2) / 2) / ((2 * pi)^(1 / 2) * (1 - alpha)))  
  }
}


