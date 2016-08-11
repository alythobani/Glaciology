function plot_switchingevents( parameters, switchesstruct )
% Plot switching events: black for switch off, magenta for switch on.

switchingevents = ...
    switchesstruct.(parameters.sensor1).(parameters.sensor2).times;

y = get(gca,'ylim');

for i = 1:length(switchingevents)
    st = switchingevents(i);
    if (switchesstruct.(parameters.sensor1).(parameters.sensor2).switchTypes(i) == -1) %switch off
        line([st st],y,'Color','k');
    else %switch on
        line([st st],y,'Color','m');
    end
end


end

