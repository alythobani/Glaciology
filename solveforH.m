function [t, H] = solveforH(h0, parameters, lambda, pressurestruct)
%solve for H = [h, dh/dh0, dh/dlambda_k]
%t will be a column of some number N times, and H will be a matrix where
% each column is a vector of N values of a certain variable (h, dh_dh0, or
% dh_dlambda_k)

fout = @(t,H) setupodes(t, H, lambda, parameters, pressurestruct);

if nargin > 3 %if lambda is passed as a parameter
    H_0 = [h0; 1; zeros(parameters.n_lambda, 1)]; %initial value for H
else %only solving for h and h_h0
    H_0 = [h0; 1]; %initial value for H
end

[t, H] = ode45(fout, parameters.tspan, H_0);

end