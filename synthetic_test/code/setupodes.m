function fout = setupodes( t, H, lambda, parameters, pressurestruct, velocity)
%Set up ODEs for H = [h; dh_dh0; dh_dlambda_k]
% fout(1) = dh/dt
% fout(2) = d/dt (h_h0)
% fout(3...kmax+2) = d/dt (h_lambda_k)
% We will pass fout into ode45(...) to solve for each of its components.
% This function is called in solveforH.m

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
    fout = zeros(parameters.n_lambda + 2, 1);
else
    fout = zeros(2,1);
end

[F, dFdh, dFdlambda] = sheetthickness(t, H, lambda, parameters, ...
    pressurestruct, velocity);

fout(1) = F; %dh/dt = F

fout(2) = dFdh*h_h0; %d/dt (h_h0) = dF/dh * h_h0

if (length(H) > 2)
    fout(3:end) = dFdlambda + ...
        dFdh*h_lambda; %d/dt(h_lambda_k) = dF/dlambda_k + (dF/dh)*h_lambda_k
end

end


