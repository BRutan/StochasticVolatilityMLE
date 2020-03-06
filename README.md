# StochasticVolatilityMLE

Black-Scholes-Merton model assumes the underlying has constant volatility. In reality, the model implied volatility used to set 
option prices exhibits mean reversion, clustering and correlation with changes in the asset.

"Financial Modeling in a Fast Mean-Reverting Stochastic Volatility Environment" (Fouque, Papanicolaou, Sircar 2/2000) implements the 
following correlated stochastic volatility Geometric Brownian Motion for the underlying asset and an Ornstein-Uhlenbeck process for 
the log implied volatility:
