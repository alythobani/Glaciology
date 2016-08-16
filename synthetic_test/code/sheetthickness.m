function [F, dFdh, dFdlambda] = ...
    sheetthickness( t, H, lambda, parameters, pressurestruct, velocity )
%Return [F; dF/dh; dF/dlambda] which are used to calculate
% the derivatives with respect to time of h, h_h0 and h_lambda_k.
% This function is called in setupodes.m

%{
INPUTS
-t is the vector of times at which we will solve for the components of H
-H is a matrix of values of h, dh/dh0 and dh/dlambda_k to be solved for by
    integration
-lambda is a vector of values for the parameters in [h_r; logk; p_ice; n_G]
    that we wish to solve for
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
-pressurestruct is a structure containing the fields:
    -> time: vector of times to go along with pressures and velocities
    -> pressure1: pressure in the first borehole (sensor1)
    -> pressure2: pressure in the second borehole (sensor2)
-velocity is an array of glacier velocity values for each time in
    pressurestruct.time
%}

h = H(1);
h_h0 = H(2);
if (length(H) > 2)
    h_lambda = H(3:end);
    dFdlambda = zeros(parameters.n_lambda, 1);
end

if parameters.h_r.in_lambda
    h_r = lambda(parameters.h_r.index);
else
    h_r = parameters.h_r.value;
end

if parameters.logk.in_lambda
    logk = lambda(parameters.logk.index);
else
    logk = parameters.logk.value;
end

if parameters.p_ice.in_lambda
    p_ice = lambda(parameters.p_ice.index);
else
    p_ice = parameters.p_ice.value;
end

if parameters.n_G.in_lambda
    n_G = lambda(parameters.n_G.index);
else
    n_G = parameters.n_G.value;
end

a = parameters.tspan(1)-5;
b = parameters.tspan(2)+5;

p_w_t = interpolate_pressure_at_t(parameters, pressurestruct, t);

%interpolate velocity vector at times in vector t
velocity_indices_in_range = find(pressurestruct.time > a & pressurestruct.time < b);
velocity_times_in_range = pressurestruct.time(velocity_indices_in_range);
u = velocity(velocity_indices_in_range);
u_t = interp1(velocity_times_in_range, u, t, 'spline');

% dh/dt = F
F = u_t.*h_r - exp(logk).*h.*(abs(p_ice - p_w_t)).^(n_G-1).*(p_ice - p_w_t);

if nargout > 1 %if we want to solve for h_h0
    dFdh = -exp(logk).*(abs(p_ice - p_w_t)).^(n_G-1).*(p_ice - p_w_t);
end

if nargout > 2 %if we want to solve for h_h0 as well as each h_lambda_k
    if (parameters.h_r.in_lambda)
        dFdlambda(parameters.h_r.index) = u_t; %dF/dh_r
    end
    
    if (parameters.logk.in_lambda)
        dFdlambda(parameters.logk.index) = ... %dF/dlogk
            -exp(logk).*h.*(abs(p_ice - p_w_t)).^(n_G-1).*(p_ice-p_w_t);
    end

    if (parameters.p_ice.in_lambda)
        dFdlambda(parameters.p_ice.index) = ... %dF/dp_ice
            -n_G.*exp(logk).*h.*(abs(p_ice - p_w_t)).^(n_G - 1);
    end

    if (parameters.n_G.in_lambda)
        dFdlambda(parameters.n_G.index) = ... %dF/dn_G
            -exp(logk).*h.*(abs(p_ice - p_w_t)).^(n_G-1).*(p_ice - ...
            p_w_t).*log(abs(p_ice-p_w_t));
    end
end
end

