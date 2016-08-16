function cost = callgaussnewton(  )
% Set up everything to call gaussnewton, then call it. Return the end cost
% after gaussnewton has converged.

%BEFORE RUNNING: execute the commands:
%   examplestruct1 = load('example1.mat');
%   examplestruct2 = load('example2.mat');
%   etc.

switchesarray(1).switchingevents = examplestruct1.data_struct.tswitch;
pressurearray(1).time = examplestruct1.data_struct.time;
pressurearray(1).pressure1 = examplestruct1.data_struct.pressure_sensor1;
pressurearray(1).pressure2 = examplestruct1.data_struct.pressure_sensor2;
thickness1 = examplestruct1.data_struct.thickness;
velocityarray(1).velocity = examplestruct1.data_struct.velocity;

switchesarray(2).switchingevents = examplestruct2.data_struct.tswitch;
time2 = examplestruct2.data_struct.time;
pressurearray(2).pressure1 = examplestruct2.data_struct.pressure_sensor1;
pressurearray(2).pressure2 = examplestruct2.data_struct.pressure_sensor2;
thickness2 = examplestruct2.data_struct.thickness;
velocityarray(2).velocity = examplestruct2.data_struct.velocity;

%{
parametersarray is an array of structures, each of which has the following fields:
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

parametersarray(1).sensor1 = 1; %will normally be something like 'S15P83'
parametersarray(1).sensor2 = 2; %will normally be something like 'S14P18'

% which sensor's pressure data will we use (or will we use both somehow?)
parametersarray(1).p_w_sensor = parametersarray(1).sensor1; %or sensor2, or 'both'

parametersarray(1).n_lambda = 3; %number of parameters (in lambda) we want to solve for

parametersarray(1).h_r.index = false; %h_r = lambda(1)
parametersarray(1).h_r.in_lambda = ...
    false; %h_r is not a parameter (in lambda) we seek to optimize
parametersarray(1).h_r.value = (6.1344e-19)*(4e5)^3/0.1;
parametersarray(1).h_r.expectedvalue = (6.1344e-19)*(4e5)^3/0.1; %lambda_0(1)
parametersarray(1).h_r.expectedvalueweight = 0; %sigma_1 

parametersarray(1).logk.index = 1; %logk = lambda(2)
parametersarray(1).logk.in_lambda = true; %k is a parameter (in lambda) we seek to optimize
parametersarray(1).logk.expectedvalue = log(6.1344e-19); %lambda_0(2) 
parametersarray(1).logk.expectedvalueweight = .5; %sigma_2 

parametersarray(1).p_ice.index = 2; %p_ice = lambda(3)
parametersarray(1).p_ice.in_lambda = ...
    true; %p_ice is a parameter (in lambda) we seek to optimize
parametersarray(1).p_ice.expectedvalue = ...
    916.7*9.80665*thickness1;
% ^ lambda_0(3) = ice overburden pressure = (density of ice)*(acceleration
%    due to gravity)*(ice thickness)
parametersarray(1).p_ice.expectedvalueweight = .5; %sigma_3

parametersarray(1).n_G.index = 3; %n_g = lambda(4)
parametersarray(1).n_G.in_lambda = true; %n_G is a parameter (in lambda) we seek to optimize
parametersarray(1).n_G.expectedvalue = 3; %lambda_0(4) 
parametersarray(1).n_G.expectedvalueweight = .5; %sigma_4 

%number of switching events
parametersarray(1).n_SE = length(switchesarray(1).switchingevents);

parametersarray(1).tspan = [pressurearray(1).time(1), pressurearray(1).time(end)];

cost = gaussnewton(parametersarray, switchesarray, pressurearray, ...
    velocityarray, length(parametersarray));



end

