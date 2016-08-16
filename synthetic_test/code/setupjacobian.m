function  J_r  = setupjacobian( h_SE, h_h0_SE, h_lambda_SE, parameters, Beta )
%Set up of the Jacobian matrix for use in Gauss-Newton method.
%   J_r(i,j) = dr_i(beta) / dbeta_j

n_SE = parameters.n_SE; % number of switching events
n_lambda = parameters.n_lambda; %number of parameters in lambda
n_C = n_SE + n_lambda;

h_c = Beta(parameters.n_lambda + 2);

J_r = zeros(n_lambda + n_SE, n_lambda + 2); %initialize Jacobian matrix


%dr_i/dlambda_k = dh(t_i)/dlambda_k for 1 <= i <= n_SE
J_r(1:n_SE, 1:n_lambda) = h_lambda_SE;


%dr_i/dlambda_(i-n_SE) = sqrt(sigma_(i-n_SE)) for n_SE < i <= n_C, and
% dr_i/dlambda_k = 0 for n_SE < i <= n_C and k != i - n_SE
%sigma_1 = parameters.h_r.expectedvalueweight;
sigma_2 = parameters.logk.expectedvalueweight;
sigma_3 = parameters.p_ice.expectedvalueweight;
sigma_4 = parameters.n_G.expectedvalueweight;
J_r(n_SE + 1 : n_C, 1:n_lambda) = diag([...%sqrt(sigma_1) ...
    sqrt(sigma_2) ...
    sqrt(sigma_3) ...
    sqrt(sigma_4)]);


%dr_i/dh0 = dh(t_i)/dh0 for 1 <= i <= n_SE
J_r(1:n_SE, n_lambda + 1) = h_h0_SE;


%dr_i/dh_c = (h_c - h(t_i)) / |h(t_i) - h_c|
J_r(1:n_SE, n_lambda + 2) = (h_c - h_SE)./abs(h_SE - h_c);

end

