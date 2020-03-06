%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GenerateStochasticVolPaths.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate large number of asset and 
% stochastic volatility paths.


function [assetpaths, volpaths] = GenerateStochasticVolPaths(assetSDE, volSDE, assetStartData, volStartData, icdf, numYears, numsteps, numpaths)

dt = numYears / numsteps;

assetpaths = [transpose(1:numsteps) zeros(numsteps, numpaths)];
volpaths = [transpose(1:numsteps) zeros(numsteps, numpaths)];

for path = 1 : numpaths    
    rng shuffle;
    iids = arrayfun(@(x) icdf(x), rand(numsteps,2));
    for step = 1 : numsteps
        if step ~= 1
            volpaths(step, path + 1) = volSDE(volpaths(step - 1, path + 1), iids(step,:), dt);
            assetpaths(step, path + 1) = assetSDE([assetpaths(step - 1, path + 1) volpaths(step - 1, path + 1)], iids(step,:), dt);
        else
            volpaths(step, path + 1) = volSDE(volStartData, iids(step,:), dt);
            assetpaths(step, path + 1) = assetSDE([assetStartData volStartData], iids(step,:), dt);
        end
    end
end

end