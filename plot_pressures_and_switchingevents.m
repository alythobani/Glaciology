function plot_pressures_and_switchingevents( parameters, switchesstruct, pressurestruct )
% Plot p1 and p2, as well as switching events: black for switch off,
% magenta for switch on.


p1 = pressurestruct.interp_data.(parameters.sensor1).pressure;
t1 = pressurestruct.interp_data.(parameters.sensor1).time;
p2 = pressurestruct.interp_data.(parameters.sensor2).pressure;
t2 = pressurestruct.interp_data.(parameters.sensor2).time;
min1 = min(t1);
min2 = min(t2);
max1 = max(t1);
max2 = max(t2);
a = max(min1, min2);
b = min(max1, max2);
p1 = p1(t1 >= a & t1 <= b);
t1 = t1(t1 >= a & t1 <= b);
p2 = p2(t2 >= a & t2 <= b);
t2 = t2(t2 >= a & t2 <= b);
plot(t1, p1, t2, p2);
hold on;
title('pressure vs time in both boreholes');
legend('p1', 'p2');

plot_switchingevents(parameters, switchesstruct);

hold off;


end

