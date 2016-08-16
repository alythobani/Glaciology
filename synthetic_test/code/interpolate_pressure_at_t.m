function [ p_w_t ] = interpolate_pressure_at_t( parameters, pressurestruct, t )
%Generate vector of pressure values at times in vector t, derived from the
% two vectors of pressure values from each borehole sensor

if (strcmp(parameters.p_w_sensor, 'both'))
    % use arithmetic mean of pressure1 and pressure2
    pressure1 = pressurestruct.pressure1;
    pressure2 = pressurestruct.pressure2;
    pressure = (pressure1 + pressure2)/2;
elseif (strcmp(parameters.p_w_sensor, parameters.sensor1))
    pressure = pressurestruct.pressure1;
else
    pressure = pressurestruct.pressure2;
end

a = parameters.tspan(1)-5;
b = parameters.tspan(2)+5;

pressure_indices_in_range = find(pressurestruct.time > a & pressurestruct.time < b);
pressure_times_in_range = pressurestruct.time(pressure_indices_in_range);
p_w = pressure(pressure_indices_in_range);
p_w_t = interp1(pressure_times_in_range, p_w, t, 'spline'); %p_w interpolated at t vector

end

