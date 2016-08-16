function [ C ] = calculatecostfunction( Beta, h_SE, parameters )
%Compute C(Beta) for a given Beta = [lambda; h0; hc]

h_c = Beta(parameters.n_lambda + 2);
h_SE = h_SE - h_c;
h_SE = (abs(h_SE)).^2;
C1 = sum(h_SE);


lambda = Beta(1:parameters.n_lambda); %lambda = [h_r; k; p_ice; n_G]
lambda_0 = [...%parameters.h_r.expectedvalue; ...
    parameters.logk.expectedvalue; ...
    parameters.p_ice.expectedvalue; ...
    parameters.n_G.expectedvalue]; %expected parameter vector for lambda
Sigma_vector = [...%parameters.h_r.expectedvalueweight ...
    parameters.logk.expectedvalueweight ...
    parameters.p_ice.expectedvalueweight ...
    parameters.n_G.expectedvalueweight];
Sigma = diag(Sigma_vector); %matrix of weights
C2 = (lambda - lambda_0)'*Sigma*(lambda - lambda_0);

C = C1 + C2;




end

