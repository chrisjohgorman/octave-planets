function [moonrise, moonset] = mriset (year, month, day, hour, latitude, longitude)

	d = day_number(year, month, day, hour);
	[x1, y1, z1, oblecl, L] = sun_rectangular(d);
	[RA, Decl, r, az, alt] = moon(d, latitude, longitude, hour);
	mpar = asind(1/r);
	GMST0 = (L + 180);
	UT_Planet_in_south = RA/15 - (L+180)/15 - longitude/15.0;
	UT_Planet_in_south = revolve_hour_angle(UT_Planet_in_south);
	cos_lha = (sind(-mpar) - sind(latitude)*sind(Decl))/(cosd(latitude) * cosd(Decl));
	if (cos_lha > 1)
		error("Mercury is always below our altitude limit.");
	elseif (cos_lha < -1)
		error("Mercury is always above our altitude limit.");
	end
	LHA = acosd(cos_lha)/15.04107;
	
	moonrise = UT_Planet_in_south - LHA;
	moonset = UT_Planet_in_south + LHA;
end
