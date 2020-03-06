%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SDELogLikelihood.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log likelihood function for stochastic differential
% equation based upon observed data and 
% transition density function.

function likelihood = SDELogLikelihood(sde, data, pdf, params)

% prior_data: [currlevel actualcurrchg dt]

likelihood = 0;

for i = 1 : length(data)
   [x_new dx] = sde(data(i, :), params);
   likelihood = likelihood + log(pdf(x_new, params));
end

end