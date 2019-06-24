function [altitude azimuth] = altitude_azimuth(sidereal_time, right_ascension, declination, latitude)
%
% usage altitude_azimuth(sidreal_time, right_ascension, declination, latitude)
%


	hour_angle = sidereal_time - right_ascension
	hour_angle = revolve_hour_angle(hour_angle)
	hour_angle = hour_angle * 15
	x = cosd(hour_angle)*cosd(declination)
	y = sind(hour_angle)*cosd(declination)
	z = sind(declination)
	x_horizon = x * sind(latitude) - z * cosd(latitude)
	y_horizon = y
	z_horizon = x * cosd(latitude) + z * sind(latitude)
	azimuth = atan2d(y_horizon,x_horizon) + 180 
	altitude = asind(z_horizon)
end
