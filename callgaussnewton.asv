function cost = callgaussnewton(  )
% Set up everything to call gaussnewton, then call it. Return the end cost
% after gaussnewton has converged.

%BEFORE RUNNING: execute the commands:
%   switchesstruct1 = load('testCase.mat');
%   pressurestruct1 = load('Interpolated data 2015 (2).mat');
%   velocity1 = load( ... ); <- import u(t) time series

%{
parameters is an array of structures, each of which has the following fields:
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
    
%}

parametersarray(1).sensor1 = 'S15P83';
parametersarray(1).sensor2 = 'S14P18';

% which sensor's pressure data will we use (or will we use both somehow?)
parametersarray(1).p_w_sensor = parametersarray(1).sensor1; %or sensor2, or 'both'

parametersarray(1).n_lambda = 3; %number of parameters we want to solve for

parametersarray(1).h_r.index = false; %h_r = lambda(1)
parametersarray(1).h_r.in_lambda = ...
    false; %h_r is not a parameter (in lambda) we seek to optimize
parametersarray(1).h_r.value = (6.1344e-19)*(4e5)^3/0.1;
parametersarray(1).h_r.expectedvalue = (6.1344e-19)*(2e5)^3/0.1; %lambda_0(1)
parametersarray(1).h_r.expectedvalueweight = 0; %sigma_1 

parametersarray(1).logk.index = 1; %logk = lambda(2)
parametersarray(1).logk.in_lambda = true; %k is a parameter (in lambda) we seek to optimize
parametersarray(1).logk.expectedvalue = log(6.1344e-19); %lambda_0(2) 
parametersarray(1).logk.expectedvalueweight = .5; %sigma_2 

parametersarray(1).p_ice.index = 2; %p_ice = lambda(3)
parametersarray(1).p_ice.in_lambda = ...
    true; %p_ice is a parameter (in lambda) we seek to optimize
parametersarray(1).p_ice.expectedvalue = ...
    916.7*9.80665*pressurestruct1.interp_data.(parametersarray(1).sensor1).position{1,1}.thickness;
% ^ lambda_0(3) = ice overburden pressure = (density of ice)*(acceleration
%    due to gravity)*(ice thickness)
parametersarray(1).p_ice.expectedvalueweight = .5; %sigma_3

parametersarray(1).n_G.index = 3; %n_g = lambda(4)
parametersarray(1).n_G.in_lambda = true; %n_G is a parameter (in lambda) we seek to optimize
parametersarray(1).n_G.expectedvalue = 3; %lambda_0(4) 
parametersarray(1).n_G.expectedvalueweight = .5; %sigma_4 

% only use switching times 1-7 for this test case
switchesstruct1.(parametersarray(1).sensor1).(parametersarray(1).sensor2).times = ...
    switchesstruct1.(parametersarray(1).sensor1).(parametersarray(1).sensor2).times(1:7);

switchingevents1 = ...
    switchesstruct1.(parametersarray(1).sensor1).(parametersarray(1).sensor2).times;

parametersarray(1).n_SE = length(switchingevents1); %number of switching events

parametersarray(1).tspan = [min(switchingevents1) - 2, max(switchingevents1) + 2];

cost = gaussnewton(parametersarray, switchesarray, pressurearray);



end

