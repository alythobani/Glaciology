function cost = gaussnewton(parametersarray, switchesarray, pressurearray)

%{
This is the implementation of the Gauss-Newton method, which we are using
to optimize parameters of the cost function C (as outlined in the project
write-up) to minimize C.
%}

lambda_guess = [... %parameters.h_r.expectedvalue; ...
    parameters.logk.expectedvalue; ...
    parameters.p_ice.expectedvalue; ...
    parameters.n_G.expectedvalue];
h0_guess = 3;
h_c_guess = 2;

Beta_0 = [lambda_guess; h0_guess; h_c_guess];
Beta_s = Beta_0;

%retrieve values of h, h_h0 and h_lambda interpolated at switching events
[h_SE, h_h0_SE, h_lambda_SE] = ...
    interpolateSEvalues(Beta_s, parameters, ...
    switchesstruct, pressurestruct);

tol = 1; %tolerance

cost = calculatecostfunction(Beta_0, h_SE, parameters); % C(Beta_0)

max_iter = 50; % maximum number of G-N iterations allowed 

iterations = 0;

disp(['cost after ' num2str(iterations) ' iterations: ' num2str(cost)]);

while (and(cost >= tol, iterations < max_iter))
    
    [h_SE, h_h0_SE, h_lambda_SE] = ...
        interpolateSEvalues(Beta_s, parameters, switchesstruct, pressurestruct);
    
    J_r = setupjacobian(h_SE, h_h0_SE, h_lambda_SE, parameters, Beta_s); %J_r^(s)
    
    Beta_s = Beta_s - ... 
        (inv(J_r'*J_r))*J_r'*calculate_r(Beta_s, h_SE, parameters); %Beta_(s)
    
    cost = calculatecostfunction(Beta_s, h_SE, parameters);
    
    iterations = iterations + 1;
    
    disp(['cost after ' num2str(iterations) ' iterations: ' num2str(cost)]);
    
end

plot_stuff(Beta_s, parameters, switchesstruct, pressurestruct);

disp(['Cost function equals ' num2str(cost) ', after ' ...
        num2str(iterations) ' iterations.']);
disp('Optimized value of Beta: ');
disp(Beta_s);

if (cost > tol)
    disp(['Error: Algorithm ran more than ' num2str(max_iter) ' iterations.']);
end
    
end