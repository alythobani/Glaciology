function [t, H] = solveforH(h0, parameters, lambda, pressurestruct, velocity)
%solve for H = [h, dh/dh0, dh/dlambda_k]
%t will be a column of some number N times, and H will be a matrix where
% each column is a vector of N values of a certain variable (h, dh_dh0, or
% dh_dlambda_k)

%{
INPUTS
-h0: initial value for h
-parameters is a structure containing the fields:
    -> n_lambda: number of parameters in the parameter vector lambda (i.e, 
          the number of parameters we wish to solve for)
    -> lambda_k (for k = 1 to 4): a structure that has fields:
        -> index: the index of lambda_k in the parameter vector lambda
            -e.g. h_r = lambda(1) = lambda_1 so h_r.index = 1
        -> in_lambda: boolean for whether lambda_k is a parameter in lambda
            -i.e. whether or not we are optimizing the value of lambda_k
                while we minimize the cost function in the gauss-newton method
        -> value: if in_lambda is false, then we have a value for lambda_k
        -> expectedvalue: an expected value for lambda_k based on known
            information about the parameter lambda_k
        -> expectedvalueweight: how confident we are in our "expected
            value" for lambda_k
    -> n_SE: number of switching events
    -> tspan: range of time we are iterating over, in the form of tspan =
            [t_init, t_end]
-lambda is a vector of values for the parameters in [h_r; logk; p_ice; n_G]
    that we wish to solve for
-pressurestruct is a structure containing the fields:
    -> time: vector of times to go along with pressures and velocities
    -> pressure1: pressure in the first borehole (sensor1)
    -> pressure2: pressure in the second borehole (sensor2)
-velocity is an array of glacier velocity values for each time in
    pressurestruct.time
%}

fout = @(t,H) setupodes(t, H, lambda, parameters, pressurestruct, velocity);

if nargin > 3 %if lambda is passed as a parameter
    H_0 = [h0; 1; zeros(parameters.n_lambda, 1)]; %initial value for H
else %only solving for h and h_h0
    H_0 = [h0; 1]; %initial value for H
end

[t, H] = ode45(fout, parameters.tspan, H_0);

end