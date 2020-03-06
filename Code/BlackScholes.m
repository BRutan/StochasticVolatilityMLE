%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BlackScholes.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function price = BlackScholes(rf, s_0, tenor, q, sig, strike, call)

d_1 = (log(s_0 / strike) + (rf - q + sig * sig / 2) * tenor) / (sig * sqrt(tenor));
d_2 = d_1 - sig * sqrt(tenor);

price = s_0 * exp(-q * tenor) * normcdf(d_1) - strike * exp(-rf * tenor) * normcdf(d_2);

if call == false
   price = price + strike * exp(-rf * tenor) - s_0 * exp(-q * tenor); 
end


end