function [sunrise, sunset] = sunrise_set (year, month, day, latitude, longitude)
%
%First we must decide which altitude we're interested in:
%
% h = 0 degrees: Center of Sun's disk touches a mathematical horizon
% h = -0.25 degrees: Sun's upper limb touches a mathematical horizon
% h = -0.583 degrees: Center of Sun's disk touches the horizon; atmospheric refraction accounted for
% h = -0.833 degrees: Sun's upper limb touches the horizon; atmospheric refraction accounted for
% h = -6 degrees: Civil twilight (one can no longer read outside without artificial illumination)
% h = -12 degrees: Nautical twilight (navigation using a sea horizon no longer possible)
% h = -15 degrees: Amateur astronomical twilight (the sky is dark enough for most astronomical observations)
% h = -18 degrees: Astronomical twilight (the sky is completely dark)
%
	h = -0.833;
	day_number = day_number(year, month, day);
	[GMST0, RA, Decl] = sun (day_number, latitude, longitude, 0); 
	UT_Sun_in_south = (RA - GMST0 - longitude) / 15.0;
	UT_Sun_in_south = revolve_hour_angle(UT_Sun_in_south);
	cos_lha = (sind(h) - sind(latitude)*sind(Decl))/(cosd(latitude) * cosd(Decl));
	if (cos_lha > 1)
		error("Sun is always below our altitude limit.");
	elseif (cos_lha < -1)
		error("Sun is always above our altitude limit.");
	end
	LHA = acosd(cos_lha)/15;
	
	time = localtime(time);
	
	sr = UT_Sun_in_south - LHA + time.gmtoff/3600;
	ss = UT_Sun_in_south + LHA + time.gmtoff/3600;
	sunrise = datestr(sr/24, 'HH:MM');
	sunset = datestr(ss/24, 'HH:MM');
end
