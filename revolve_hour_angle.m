function [h] = revolve_hour_angle(hour)
    h = hour - floor(hour/24)*24;
end
