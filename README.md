# volatility-sensitive-Bayesian
This repository contains code related to the research paper 'Volatility Sensitive Bayesian Estimation of Portfolio VaR and CVaR' by T. Bodnar, V. Niklasson and E. Thorsén. The purpose of the code is to show how the methods and numerical study in this paper have been implemented in R.

The file main.R controls what to run.

The file portfolio_VaRs.R is used by main.R in order to compute daily portfolio VaRs.

The file calculations.R handles different computations in portfolio_VaRs.R.

To run the main file, you must first download stock prices and load them at the specified place in the main file (see further instructions in the main file).
