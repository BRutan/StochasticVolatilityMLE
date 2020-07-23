# StochasticVolatilityMLE

Black-Scholes-Merton model assumes the underlying has constant volatility. In reality, the model implied volatility used to set 
option prices exhibits mean reversion, clustering and correlation with changes in the asset price.

"Financial Modeling in a Fast Mean-Reverting Stochastic Volatility Environment" (Fouque, Papanicolaou, Sircar; Feb 2000) implements the 
following correlated stochastic volatility Geometric Brownian Motion for the underlying asset and an Ornstein-Uhlenbeck process for 
the log implied volatility:

![alt text](https://raw.githubusercontent.com/BRutan/StochasticVolatilityMLE/master/Assets/SDEs.png)

