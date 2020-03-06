%%%%%%%%%%%%%%%%%%%%%%%%%
% HypTest.m
%%%%%%%%%%%%%%%%%%%%%%%%
% Test hypothesis regarding model parameters.

function resultStruct = HypTest(estParams, hypStruct, errMat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expecting following inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: Estimated parameters using MLE
% hypStruct: Struct containing hypothesis testing parameters:
% - .alpha: 1 - confidence%, must be in (0, .5].
% - .nullparams: 1xn of null values for all estimators, with n
% corresponding to number of parameters.
% - (Optional):
% - .data: Data for generating standard error matrix
% - .sympdf: Symbolic version of data pdf.
% - .symbols: Symbols used in sympdf.
% Optional:
% errMat: nxn standard error matrix, with n = number of parameters.

if nargin < 2
   error("Require estParams, hypStruct as parameters.")
end

alpha = hypStruct.alpha;
nulls = hypStruct.nullparams;
% Generate covariance matrix using observed data, or use provided one:
if nargin ~= 3 || isempty(errMat)
    errMat = StdErrMatrix(hypStruct.data, hypStruct.sympdf, hypStruct.symbols, estParams);
    n = length(hypStruct.data);
end
% Under two-tailed hypothesis testing, calculate all pvalues of observing the maximum 
% likelihood estimators of the distribution.
results.confIntervals = zeros(length(nulls), 2);
results.pvals = zeros(length(nulls), 1);
results.stderrs = zeros(length(nulls), 1);
results.t_obs = zeros(length(nulls), 1);
results.acceptreject = zeros(length(nulls), 1);

for i = 1:length(nulls)
    std_err = sqrt(abs(errMat(i, i)) / n);
    if isempty(n)
        t_crit = abs(norminv(alpha / 2));
    else
        t_crit = abs(tinv(alpha / 2, n - length(nulls)));
    end
    t_obs = (estParams(i) - nulls(i)) / std_err;
    results.stderrs(i) = std_err;
    results.t_obs(i) = t_obs;
    interval = [estParams(i) - std_err * t_crit, estParams(i) + std_err * t_crit];
    results.confIntervals(i, 1) = interval(1);
    results.confIntervals(i, 2) = interval(2);
    if t_obs <= 0
        results.pvals(i) = tcdf(t_obs, n - length(nulls)) * 2;
    else
        results.pvals(i) = (1 - tcdf(t_obs, n - length(nulls))) * 2;    
    end
    if t_crit < abs(results.t_obs(i))
        results.acceptreject(i) = 1;
    else
        results.acceptreject(i) = 0;
    end
end

resultStruct = results;

end

function matr = StdErrMatrix(data, density, symbols, obsParams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expecting following inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data: mxn matrix containing generated data from distribution we seek to estimate.
% m > 1 for markov processes.
% density: Symbolic density function.
% symbols: variables of density function.
% obsParams: Hypothesis testing parameters or observed parameters of distribution.

% Generate the Covariance Matrix:
n = length(data);
covMatr = zeros(length(obsParams), length(obsParams));
dims = size(covMatr);
hessianF = hessian(log(density), sym(symbols));
hessStr = string(hessianF);
hessStr = split(hessStr, ';');
hessStr = split(hessStr, ',');
cleanedHessians = arrayfun(@(str) lower(regexprep(str, "(\[|\])", '')),hessStr);
for row = 1 : dims(1)
   for col = 1 : dims(2)
       func = ConvertSymbolicFunction(cleanedHessians(row, col));
       covMatr(row, col) = 0;
       for obs = 1 : length(data)
           covMatr(row, col) = covMatr(row, col) + func(data, obsParams) / n;
       end
       covMatr(row, col) = -1 / covMatr(row, col);
   end
end

matr = covMatr;

end