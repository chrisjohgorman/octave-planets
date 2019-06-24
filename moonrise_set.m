function [moonise, moonset] = moonrise_set (year, month, day, latitude, longitude)
	h = -0.833;
	day_number = day_number(year, month, day);
	[GMST0, RA, Decl] = moon (day_number, latitude, longitude, 0); 
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
	
	mr = UT_Sun_in_south - LHA + time.gmtoff/3600;
	ms = UT_Sun_in_south + LHA + time.gmtoff/3600;
	moonrise = datestr(sr/24, 'HH:MM');
	moonset = datestr(ss/24, 'HH:MM');
end
