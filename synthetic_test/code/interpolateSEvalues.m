function [h_SE, h_h0_SE, h_lambda_SE ] = ...
    interpolateSEvalues( Beta_s, parameters, switchingevents, ...
    pressurestruct, velocity )
% Retrieve interpolated values of h, h_h0 and h_lambda at switching event
% times.

%{
INPUTS
-Beta_s = [lambda_s; h0_s; hc_s] is the current values of parameters we're
    solving for
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
-switchingevents is an array of switching event times
-pressurestruct is a structure containing the fields:
    -> time: vector of times to go along with pressures and velocities
    -> pressure1: pressure in the first borehole (sensor1)
    -> pressure2: pressure in the second borehole (sensor2)
-velocity is an array of glacier velocity values for each time in
    pressurestruct.time
%}

lambda = Beta_s(1:parameters.n_lambda); %w is the only parameter in lambda
h0_s = Beta_s(parameters.n_lambda+1);
hc_s = Beta_s(parameters.n_lambda+2);
[t, H] = solveforH(h0_s, parameters, lambda, pressurestruct, velocity);

% interpolate h at switching times
% This will return a vector, h_SE, of values of h evaluated at
% switching times
h = H(:, 1);
h_SE = interp1(t, h, switchingevents, 'spline');


% interpolate h_h0 at switching times
% This will return a vector, h_h0_SE, of values of dh/dh0 evaluated at
% switching times
h_h0 = H(:, 2);
h_h0_SE = interp1(t, h_h0, switchingevents, 'spline');

% interpolate h_lambda at switching times
% This will return a n_SE x n_lambda matrix, h_lambda_SE, whose kth column
% will be a vector of values of dh/dlambda_k evaluated at switching times
h_lambda = H(:, 3:end);
h_lambda_SE = interp1(t, h_lambda, switchingevents, 'spline');

end

