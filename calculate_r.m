function [ R ] = calculate_r( Beta, h_SE, parameters )
%For a given vector Beta = [lambda; h0; hc], calculate
%   R = [r_1(Beta); r_2(Beta); ... ; r_n_C(Beta)] where n_C = n_SE + n_lambda
% We have that r_i(Beta) = {
%   |h(t_i) - h_c| if 1 <= i <= n_SE
%   sqrt(sigma_(i-n_SE))*(lambda(i-n_SE)-lambda_0(i-n_SE)) if n_SE < i <= n_C
% }.

h_c = Beta(parameters.n_lambda + 2);
n_C = parameters.n_SE + parameters.n_lambda;
R = zeros(n_C,1); %initialize R vector
R(1:parameters.n_SE) = abs(h_SE - h_c);

Sigma_vector = [...%parameters.h_r.expectedvalueweight; ...
    parameters.logk.expectedvalueweight; ...
    parameters.p_ice.expectedvalueweight; ...
    parameters.n_G.expectedvalueweight]; %vector of weights
lambda = Beta(1:parameters.n_lambda); %lambda = [h_r; k; p_ice; n_G]
lambda_0 = [...%parameters.h_r.expectedvalue; ...
    parameters.logk.expectedvalue; ...
    parameters.p_ice.expectedvalue; ...
    parameters.n_G.expectedvalue]; %expected values for lambda

R(parameters.n_SE + 1 : end) = sqrt(Sigma_vector).*(lambda - lambda_0);







end

