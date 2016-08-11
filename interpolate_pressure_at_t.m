function [ p_w_t ] = interpolate_pressure_at_t( parameters, pressurestruct, t )
%Generate vector of pressure values at times in vector t, derived from the
% two vectors of pressure values from each borehole sensor

if (strcmp(parameters.p_w_sensor, 'both'))
    % TODO: find a way to combine both sensors' pressure data sets, and
    % generate 'pressure' and 'pressuretimes' vectors
    pressure1 = pressurestruct.interp_data.(parameters.sensor1).pressure;
    pressuretimes1 = pressurestruct.interp_data.(parameters.sensor1).time;
    pressure2 = pressurestruct.interp_data.(parameters.sensor2).pressure;
    pressuretimes2 = pressurestruct.interp_data.(parameters.sensor2).time;
elseif (strcmp(parameters.p_w_sensor, parameters.sensor1))
    pressure = pressurestruct.interp_data.(parameters.sensor1).pressure;
    pressuretimes = pressurestruct.interp_data.(parameters.sensor1).time;
else
    pressure = pressurestruct.interp_data.(parameters.sensor2).pressure;
    pressuretimes = pressurestruct.interp_data.(parameters.sensor2).time;
end

a = parameters.tspan(1)-5;
b = parameters.tspan(2)+5;

pressureSEindices = find(pressuretimes > a & pressuretimes < b);
pressureSEtimes = pressuretimes(pressureSEindices);
p_w = pressure(pressureSEindices);
p_w_t = interp1(pressureSEtimes, p_w, t); %p_w interpolated at t vector

end

