%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GenerateSDEData.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate data from SDE using data
% set.

function forecasts = GenerateSDEStepData(sde, priorData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expecting following parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sde: Unary function with handle @(data), i.e. parameters
% are fixed.
% priorData: nxm array containing prior data to generate next
% step data with.

n = length(priorData);
forecasts = zeros(n, 1);
for i = 1 : n
    forecasts(i) = sde(priorData(i, :));
end

end