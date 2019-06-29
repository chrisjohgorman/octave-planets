function [moonrise, moonset] = moonrise_set (year, month, day, hour, latitude, longitude)

	d = day_number(year, month, day, hour);
	[x1, y1, z1, oblecl, L] = sun_rectangular(d)
	[RA, Decl, r, az, alt] = moon(d, latitude, longitude)
	mpar = asind(1/r)
	GMST0 = (L + 180)
	UT_Moon_in_south = RA/15 - (L+180)/15 - longitude/15.0
	cos_lha = (sind(-mpar) - sind(latitude)*sind(Decl))/(cosd(latitude) * cosd(Decl))
	if (cos_lha > 1)
		error("Moon is always below our altitude limit.");
	elseif (cos_lha < -1)
		error("Moon is always above our altitude limit.");
	end
	LHA = acosd(cos_lha)/15.04107
	
	time = localtime(time);
	
	mr = UT_Moon_in_south - LHA
	ms = UT_Moon_in_south + LHA
	%mr = UT_Moon_in_south - LHA + time.gmtoff/3600;
	%ms = UT_Moon_in_south + LHA + time.gmtoff/3600;
	moonrise = mr;
	moonset = ms;
	%moonrise = datestr(mr/24, 'HH:MM:SS');
	%moonset = datestr(ms/24, 'HH:MM:SS');
end
