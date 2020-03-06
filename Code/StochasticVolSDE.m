%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% StochasticVolSDE.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wrapper for generic stochastic volatility
% stochastic differential equation.

function [asset_new, vol_new] = StochasticVolSDE(assetSDE, volSDE, assetdata, voldata)

vol_new = volSDE(voldata);
asset_new = assetSDE(assetdata, vol_new);

end