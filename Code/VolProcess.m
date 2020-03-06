


function [logiv_new, dlogiv] = VolProcess(data, params)

% data = [curr_logiv actual_dlogiv dt dWt/dt]
% params = [alpha beta mu rho]

curr_logiv = data(1);
actual_dlogiv = data(2);
dt = data(3);
dWt = data(4) * sqrt(dt);
alpha = params(1);
beta = params(2);
mu = params(3);
rho = params(4);
dZt = (actual_dlogiv - alpha * (mu - curr_logiv) * dt) / beta;
dZt = dZt - rho * dWt;
dZt = dZt / sqrt(1 - rho * rho);
dlogiv = alpha * (mu - curr_logiv) * dt;
dlogiv = dlogiv + beta * (rho * dWt + sqrt(1 - rho * rho) * dZt);
logiv_new = curr_logiv + dlogiv;

end