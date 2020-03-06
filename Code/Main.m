%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main method for Financial Econometrics II Final Project.

import ConvertSymbolicFunction.*;
import GenerateSDEPaths.*;
import GenerateSDEStepData.*;
import GenerateStochasticVolPaths.*;
import HypTest.*;
import MLEEstimation.*;
import SDELogLikelihood.*;
import TransitionDensity.*;
import VolProcess.*;

%%%%%%%%%%%%
% MLE Data:
%%%%%%%%%%%%
% Pull local S&P 500 index price, implied volatility data:
assetData = readtable('^GSPC.csv');
ivData = readtable('^VIX.csv');
% Divide implied volatility and corresponding asset data into
% periods exhibiting obvious mean reversion:
startdate = '2/15/2018';
enddate = '2/26/2020';
daysPerYear = 365;
% Assume 1.8% annual dividend yield.
divYield = .018;
tenor = (datenum(enddate) - datenum(startdate)) / daysPerYear;
ivPeriod = ivData(ivData.Date >= startdate & ivData.Date <= enddate, {'Date', 'AdjClose'});
periodDates = ivPeriod.Date;
assetData = table2array(assetData(assetData.Date >= startdate & assetData.Date <= enddate, {'AdjClose'}));
assetLogReturns = log(assetData(2:length(assetData))) - log(assetData(1:length(assetData) - 1));
ivData = arrayfun(@(x) x / 100, ivPeriod.AdjClose);
ivLogReturns = log(ivData(2:length(ivData))) - log(ivData(1:length(ivData)-1));
ivLogs = log(ivData);
actualdIVLogs = diff(ivLogs);
% Calculate parameters estimated from market data:
asset_mu = mean(assetLogReturns) * daysPerYear - divYield; 
asset_sigma = sqrt(var(assetLogReturns) * daysPerYear);
asset_norm_devs = zeros(length(assetLogReturns), 1);
for i = 1 : length(assetLogReturns)
   asset_norm_devs(i) = (assetLogReturns(i) - asset_mu / daysPerYear) / (exp(ivLogs(i)) * sqrt(1 / daysPerYear));
end
iv_mu = mean(ivLogs);
rho = corrcoef(assetLogReturns,actualdIVLogs);
rho = rho(1,2);
%%%%%%%%%%%%
% Implied Volatility Process:
%%%%%%%%%%%%
% Generate Ornstein-Uhlenbeck process using estimated parameters, to be MLE
% estimated:
% (x_1: Current Y level).
% (x_2: dW_t (Wiener process for underlying stock).
% (x_3: dZ_t (Wiener process for underlying volatility)) 
% (x_4: dt).
syms x_1 x_2 x_3 x_4 p_1 p_2
z(x_1, x_2)=sym(rho)*x_1+(1-sym(rho)^2)^sym(.5)*x_2;
sde_iv(x_1,x_2,x_3,x_4,p_1,p_2)=(p_1*(sym(iv_mu)-x_1)+p_2*z(x_2, x_3))*x_4;
pdf_iv_var = @(params) params(2) ^ 2 / (2 * params(1));
pdf_iv = @(data, params) normpdf(data(1), iv_mu, sqrt(pdf_iv_var(params))); 
pdf_iv_var_sym(p_1, p_2) = p_2^sym(2) / (sym(2) * p_1);
pdf_iv_sym(x_1, p_1, p_2) = 1 / sym(sqrt(2 * pi) * pdf_iv_var_sym(p_1, p_2) ^ sym(.5)) * sym(exp(1)) ^ ((x_1 - sym(iv_mu))^sym(2)) / (sym(2) * pdf_iv_var_sym(p_1, p_2)); 
dZ_iv_handle = @(data, params) (data(1) - params(1) * (iv_mu - data(2)) * data(3)) / (params(2) * sqrt(data(3))); 
% sde_iv_handle = @(data, params) data(1) + (params(1) * (iv_mu - data(1)) + params(2) * data(3)) * dZ_iv([data(1) data(2)], params);
% Calibrate implied volatility stochastic process to data:
% (find alpha, beta):
startParams = [50 5]; 
n = length(ivLogs);
estimData = [ivLogs(2:n) actualdIVLogs arrayfun(@(x) tenor/daysPerYear, zeros(n - 1, 1)) asset_norm_devs];
sde_iv_handle = @(data, params) VolProcess(data, [params iv_mu rho]);
lFunc = @(params) SDELogLikelihood(sde_iv_handle, estimData, pdf_iv, params);
[mleParams, llikelihood] = MLEEstimation(estimData, pdf_iv, startParams, lFunc);
sde_y_fixed = @(data) sde_iv_handle(data, mleParams);
iv_forecasts = GenerateSDEStepData(sde_y_fixed, estimData);
% Convert to % implied volatilities:
iv_forecasts = arrayfun(@(x) exp(x), iv_forecasts);
actuals = ivData(1:n - 1);
% Chart forecasted and actual data, print as data:
fig = figure;
plot(periodDates(2:n), actuals, periodDates(2:n), iv_forecasts);
title("VIX Actual vs MLE Forecasts");
xlabel("Date");
ylabel("VIX");
print("VIX Actual vs MLE Forecasts", "-dpng");
close(fig);
dateStrs = arrayfun(@(x) datestr(x), periodDates(2:n), 'UniformOutput', false);
results = array2table(dateStrs, 'VariableNames', {'Date'});
results = [results array2table(actuals) array2table(iv_forecasts)];
results.Properties.VariableNames = {'Date', 'Actual_IV', 'Projected_IV'};
writetable(results, "ImpliedVolatilityForecastsVsActual.csv");
% Generate and plot a few implied volatility sample paths:
dZ_t = @(iid, dt) rho * iid(1) * sqrt(dt) + sqrt(1 - rho * rho) * iid(2) * sqrt(dt);
sde_iv_path = @(logiv, iids, dt) logiv + (mleParams(1) * (iv_mu - logiv) * dt + mleParams(2) * dZ_t(iids, dt));
icdf = @(unif) norminv(unif, 0, 1);
startData = [ivLogs(1)];
optParams.fullpaths = true;
optParams.numiids = 2;
paths = GenerateSDEPaths(sde_iv_path, icdf, 1, 1, 1000, startData, optParams);
paths = [paths(:, 1) arrayfun(@(x) exp(x), paths(:, 2))];
fig = figure;
plot(paths(:,1), paths(:,2));
title("VIX Simulated Path");
xlabel("Step");
ylabel("VIX");
print("VIX Simulated Path", "-dpng");
close(fig);
% Calculate statistical significance of results:
hypStruct.alpha = .05;
hypStruct.nullparams = [0 0];
hypStruct.data = estimData;
hypStruct.sympdf = pdf_iv_sym;
hypStruct.symbols = [p_1 p_2];
resultStruct = HypTest(mleParams, hypStruct);
resultStruct.estimators = mleParams;
resultStruct.likelihood = llikelihood;
resultStruct.nulls = hypStruct.nullparams;
resultStruct.hypothesisparams = hypStruct;
resultStruct.paramNames = ["alpha" "beta"];
PrintMLESummary(resultStruct, "HypothesisTesting-IVSDE.csv");
%%%%%%%%%%%%
% Underlying Stock Process:
%%%%%%%%%%%%
% Create modified GBM sde for underlying stock process
% dXt/Xt = mu dt + f(Yt) dWt
asset_sde = @(data, iids, dt) data(1) * (1 + asset_mu * dt + exp(data(2)) * iids(1) * sqrt(dt)); 
%%%%%%%%%%%%%%
% Generate large number of implied volatility sample paths,
% underlying sample paths given volatility sample paths,
% to price an option following stochastic volatility geometric brownian motion:
%%%%%%%%%%%%%% 
% Price ^GSPC option expiring 6/19/2020:
expDate = '6/19/2020';
today = '3/4/2020';
priceToday = 3130.12;
ivToday = .2203124351;
ivLogToday = log(ivToday);
rf = 0.68 / 100;
strike = 3125;
tenor = (datenum(expDate) - datenum(today)) / daysPerYear;
midPrice = (142.1999969 + 144.3999939) / 2;
blackScholesPrice = BlackScholes(rf, priceToday, tenor, divYield, ivToday, strike, true);
numpaths = 3000;
numsteps = 2500;
[assetpaths, volpaths] = GenerateStochasticVolPaths(asset_sde, sde_iv_path, priceToday, ivLogToday, icdf, tenor, numsteps, numpaths);
% Plot asset and volatility path:
fig = figure;
tiledlayout(2,1);
nexttile;
plot(assetpaths(:,1), assetpaths(:,2));
title("Asset Simulated Path");
xlabel("Step #");
ylabel("Asset Price");
nexttile;
plot(assetpaths(:,1), arrayfun(@(x) exp(x), volpaths(:,2)));
xlabel("Step #");
ylabel("IV");
title("IV Simulated Path");
print("IV And Asset Simulated Path", "-dpng");
close(fig);
filename = "GSPC Option Pricing.xlsx";
discPayoffFunc = @(price) max(price - strike, 0) * exp(-rf * tenor);
discountedPayoffs = arrayfun(@(val) discPayoffFunc(val), transpose(assetpaths(numsteps, :)));
optionPrice = sum(discountedPayoffs) / length(discountedPayoffs);
optionPricingTable = array2table([tenor, ivToday, priceToday, divYield, rf, strike, midPrice, blackScholesPrice, optionPrice]);
optionPricingTable.Properties.VariableNames = {'T','Sigma_0','S_0','DivYield','RiskFreeRate','Strike','MarketMid','BlackScholesValue','StochModelValue'};
writetable(optionPricingTable, filename, 'Sheet', 1, 'Range', 'A1');