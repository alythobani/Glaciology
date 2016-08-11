% Version 5, where k is now logk, and the formula for dh/dt = F is now
% changed

%BEFORE RUNNING: execute the commands:
%   switchesstruct = load('testCase.mat');
%   pressurestruct = load('Interpolated data 2015 (2).mat');
%   velocity = load( ... ); <- import u(t) time series

%{
This is the implementation of the Gauss-Newton method, which we are using
to optimize parameters of the cost function C (as outlined in the project
write-up) to minimize C.
%}

%{
parameters is a structure that has the following fields:
    -> n_lambda: number of parameters in the parameter vector lambda
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
    
%}

parameters.sensor1 = 'S15P83';
parameters.sensor2 = 'S14P18';
parameters.p_w_sensor = parameters.sensor1; %or sensor2, or 'both'

parameters.n_lambda = 3; %number of parameters we want to solve for

parameters.h_r.index = false; %h_r = lambda(1)
parameters.h_r.in_lambda = ...
    false; %h_r is not a parameter (in lambda) we seek to optimize
parameters.h_r.value = (6.1344e-19)*(4e5)^3/0.1;
parameters.h_r.expectedvalue = (6.1344e-19)*(2e5)^3/0.1; %lambda_0(1)
parameters.h_r.expectedvalueweight = 0; %sigma_1 

parameters.logk.index = 1; %logk = lambda(2)
parameters.logk.in_lambda = true; %k is a parameter (in lambda) we seek to optimize
parameters.logk.expectedvalue = log(6.1344e-19); %lambda_0(2) 
parameters.logk.expectedvalueweight = .5; %sigma_2 

parameters.p_ice.index = 2; %p_ice = lambda(3)
parameters.p_ice.in_lambda = ...
    true; %p_ice is a parameter (in lambda) we seek to optimize
parameters.p_ice.expectedvalue = ...
    916.7*9.80665*pressurestruct.interp_data.(parameters.sensor1).position{1,1}.thickness;
% ^ lambda_0(3) = ice overburden pressure = (density of ice)*(acceleration
%    due to gravity)*(ice thickness)
parameters.p_ice.expectedvalueweight = .5; %sigma_3

parameters.n_G.index = 3; %n_g = lambda(4)
parameters.n_G.in_lambda = true; %n_G is a parameter (in lambda) we seek to optimize
parameters.n_G.expectedvalue = 3; %lambda_0(4) 
parameters.n_G.expectedvalueweight = .5; %sigma_4 

switchesstruct.(parameters.sensor1).(parameters.sensor2).times = ...
    switchesstruct.(parameters.sensor1).(parameters.sensor2).times(1:7);
switchingevents = ...
    switchesstruct.(parameters.sensor1).(parameters.sensor2).times;
parameters.n_SE = length(switchingevents); %number of switching events

parameters.tspan = [min(switchingevents) - 2, max(switchingevents) + 2];

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
    

