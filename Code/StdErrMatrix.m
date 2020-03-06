%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% StdErrMatrix.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate estimated standard error matrix
% based upon data and density function 
% using Taylor Series expansion.

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
hessianF = hessian(density, sym(symbols));
hessStr = string(hessianF);
hessStr = split(hessStr, ';');
hessStr = split(hessStr, ',');
cleanedHessians = arrayfun(@(str) lower(regexprep(str, "(\[|\])", '')),hessStr);
for row = 1 : dims(1)
   for col = 1 : dims(2)
       func = ReplaceSymbols(cleanedHessians(row, col));
       covMatr(row, col) = 0;
       for obs = 1 : length(data)
           covMatr(row, col) = covMatr(row, col) + func(data, obsParams) / n;
       end
       covMatr(row, col) = -1 / covMatr(row, col);
   end
end

matr = covMatr;

end