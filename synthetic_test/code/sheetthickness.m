function [F, dFdh, dFdlambda] = ...
    sheetthickness( t, H, lambda, parameters, pressurestruct )
%Return [F; dF/dh; dF/dlambda] which are used to calculate
% the derivatives with respect to time of h, h_h0 and h_lambda_k.
% This function is called in setupodes.m

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
% disp('pressure diff:');
% disp(p_ice - p_w_t);

% assuming we have a vector called velocity, and a vector called
% velocitytimes

%sample velocity data (made up)
%TODO: get real velocity data
u_t = .1;

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

