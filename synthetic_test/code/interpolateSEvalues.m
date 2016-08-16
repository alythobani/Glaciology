function [h_SE, h_h0_SE, h_lambda_SE ] = ...
    interpolateSEvalues( Beta_s, parameters, switchesstruct, pressurestruct )
% Retrieve interpolated values of h, h_h0 and h_lambda at switching event
% times.

switchingevents = ...
    switchesstruct.(parameters.sensor1).(parameters.sensor2).times;

lambda = Beta_s(1:parameters.n_lambda); %w is the only parameter in lambda
h0_s = Beta_s(parameters.n_lambda+1);
hc_s = Beta_s(parameters.n_lambda+2);
[t, H] = solveforH(h0_s, parameters, lambda, pressurestruct);

h = H(:,1);
% 
% %PLOT 1%
% subplot(2,2,1);
% y_hc = hc_s*ones(size(t));
% plot(t,h,t,y_hc,'--');
% hold on;
% title('h vs t');
% legend('h', 'h_c');
% plot_switchingevents(parameters, switchesstruct);
% hold off;
% 
% 
% %PLOT 2: Pressures and switching events%
% subplot(2,2,2);
% plot_pressures_and_switchingevents(parameters, switchesstruct, pressurestruct);
% 
% 
% %PLOT 3: Cost function with respect to hc and logk
% subplot(2,2,[3 4]);
% logk_s = Beta_s(parameters.logk.index);
% plot_cost_function_wrt_2_vars(switchingevents, pressurestruct, Beta_s, ...
%     parameters, parameters.logk.index, parameters.n_lambda+2, ...
%     [logk_s - log(100), logk_s + log(100)], [hc_s - 10, hc_s + 10], 1);


% interpolate h at switching times
% This will return a vector, h_SE, of values of h evaluated at
% switching times
h = H(:, 1);
h_SE = interp1(t, h, switchingevents, 'spline');


% interpolate h_h0 at switching times
% This will return a vector, h_h0_SE, of values of dh/dh0 evaluated at
% switching times
h_h0 = H(:, 2);
h_h0_SE = interp1(t, h_h0, switchingevents, 'spline');

% interpolate h_lambda at switching times
% This will return a n_SE x n_lambda matrix, h_lambda_SE, whose kth column
% will be a vector of values of dh/dlambda_k evaluated at switching times
h_lambda = H(:, 3:end);
h_lambda_SE = interp1(t, h_lambda, switchingevents, 'spline');

end

