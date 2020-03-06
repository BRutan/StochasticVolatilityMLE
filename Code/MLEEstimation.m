%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MLEEstimation.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% * Maximum likelihood estimation for parameters given
% probability density function.

function [params, llikelihood] = MLEEstimation(data, pdf, startParams, logFunc, geq, leq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expecting following inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data: nxm matrix of randomly sampled data from some distribution
% (doubles) to estimate distribution parameter using MLE.
% pdf: Symbolic version of PDF data was sampled from.
% symbols: 1xm Vector of symbols used in pdfSym.
% startParams: 1xm vector of MLE start values for estimators.
% logFunc: Optional log likelihood function override.

if nargin < 3
   error('Require data, pdf, pdfSym, hypStruct, symbols as arguments.')
elseif nargin < 4
    % Use prespecified log function (assumes constant parameters):
   logFunc = @(params) LogLikelihood(params, data, pdf);
end

if nargin < 5
    % Maximize the function with the given data set (unconstrained):
    [estimators, logLikelihood] = Maximize(logFunc, startParams);
elseif nargin == 5
    % Use leq constraint only:
    [estimators, logLikelihood] = Maximize(logFunc, startParams, geq);
else
    % Use leq and geq constraints:
    [estimators, logLikelihood] = Maximize(logFunc, startParams, geq, leq);
end

params = estimators;
llikelihood = logLikelihood;

end

function [result, functionVal] = Maximize(func, start, geq, leq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expecting following inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% func: Function object to maximize. Must take only one parameter,
% the arguments in vector form to maximize the function.
% start: 1 x m array containing bounds to maximize function over
% (m must correspond to number of parameters being estimated). 

minFunc = @(x)-func(x);
if nargin < 3
    [x1, x2] = fminsearch(minFunc, start);
elseif nargin == 3
    % Use constraints with minimization method:
    [x1, x2] = fmincon(minFunc, start, [], [], [], [], geq);
elseif nargin == 4
    % Use constraints with minimization method:
    [x1, x2] = fmincon(minFunc, start, [], [], [], [], geq, leq);
end
result = x1;
functionVal = -x2;

end

function [likelihood] = LogLikelihood(estimator_vec, data, densityFunc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expecting following inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimator_vec: 1xk vector of estimators.
% data: n x m data set drawn from densityFunc. m > 1 for markov data.
% densityFunc: Symbolic PDF. 
out = 0;
dims = size(data);
numRows = dims(1);
useSheets = 0;
if length(dims) > 2
   useSheets = 1;
end
for row = 1:numRows
    if useSheets == 0
        x = data(row, :);
    else
        x = data(row, :, :);
    end
    out = out + log(densityFunc(x, estimator_vec));
end

likelihood = out;

end

