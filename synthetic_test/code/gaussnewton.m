function cost = gaussnewton(parametersarray, switchesarray, pressurearray, n)

%{
This is the implementation of the Gauss-Newton method, which we are using
to optimize parameters of the cost function C (as outlined in the project
write-up) to minimize C.
%}

%{
INPUTS
-parametersarray: array of parameter structures, one for each pair of boreholes
    -each parameter structure has the fields:
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
-switchesarray: array of switching event structures, one for each pair of
        boreholes
    -each switching event structure has the fields:
        -> switchingevents: array of switching events
-pressurearray: array of pressure structures, one for each pair of boreholes
    -each pressure structure has the fields:
        -> pressure1: pressure in the first borehole (sensor1)
        -> pressure2: pressure in the second borehole (sensor2)
-n: number of pairs of boreholes being used
%}

for i = 1:n
    lambda_guess(i).vector = [];
    if (parametersarray(i).h_r.in_lambda)
        lambda_guess(i).vector = [lambda_guess(i).vector; ...
            parametersarray(i).h_r.expectedvalue];
    end
    if (parametersarray(i).logk.in_lambda)
        lambda_guess(i).vector = [lambda_guess(i).vector; ...
            parametersarray(i).logk.expectedvalue];
    end
    if (parametersarray(i).p_ice.in_lambda)
        lambda_guess(i).vector = [lambda_guess(i).vector; ...
            parametersarray(i).p_ice.expectedvalue];
    end
    if (parametersarray(i).n_G.in_lambda)
        lambda_guess(i).vector = [lambda_guess(i).vector; ...
            parametersarray(i).n_G.expectedvalue];
    end
    h0_guess = 3;
    h_c_guess = 2;
    Beta_0(i).vector = [lambda_guess(i).vector; h0_guess; h_c_guess];
    Beta_s(i).vector = Beta_0(i).vector;
    
    %retrieve values of h, h_h0 and h_lambda interpolated at switching events
    [h_SE(i).vector, h_h0_SE(i).vector, h_lambda_SE(i).vector] = ...
        interpolateSEvalues(Beta_s, parametersarray(i), ...
        switchesarray(i), pressurearray);
    
end

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