function plot_stuff(Beta_s, parameters, switchesstruct, pressurestruct)
% Plot stuff.

switchingevents = ...
    switchesstruct.(parameters.sensor1).(parameters.sensor2).times;

lambda = Beta_s(1:parameters.n_lambda); %w is the only parameter in lambda
h0_s = Beta_s(parameters.n_lambda+1);
hc_s = Beta_s(parameters.n_lambda+2);
[t, H] = solveforH(h0_s, parameters, lambda, pressurestruct);

h = H(:,1);


%PLOT 1: h vs t, with switching events overlayed, and h_c overlayed%
subplot(2,1,1);
y_hc = hc_s*ones(size(t));
plot(t,h,t,y_hc,'--');
hold on;
title('h vs t');
legend('h', 'h_c');
plot_switchingevents(parameters, switchesstruct);
limit1 = axis;
hold off;


%PLOT 2: Pressures and switching events%
subplot(2,1,2);
hold on;
plot_pressures_and_switchingevents(parameters, switchesstruct, pressurestruct);
limit2 = axis;
limit2(1:2) = limit1(1:2);
axis(limit2);
hold off;



%PLOT 3: Cost function with respect to hc and logk
figure;
logk_s = Beta_s(parameters.logk.index);
plot_cost_function_wrt_2_vars(switchingevents, pressurestruct, Beta_s, ...
    parameters, parameters.logk.index, parameters.n_lambda+2, ...
    [logk_s - log(100), logk_s + log(100)], [hc_s - 10, hc_s + 10], .1);


end

