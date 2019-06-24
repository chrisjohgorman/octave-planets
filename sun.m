function [azimuth, altitude, RA, Decl, L] = sun(day_number, latitude, longitude, UT)
	% rotate equitorial coordinates
	[x1, y1, z1, oblecl, L] = sun_rectangular(day_number);
	xequat = x1;
	yequat = y1 * cosd(oblecl) - z1 * sind(oblecl);
	zequat = y1 * sind(oblecl) - z1 * cosd(oblecl);

	RA = atan2d(yequat, xequat);
	RA = revolve_degree(RA);
	% convert RA to hours
	RA = RA / 15;
	Decl = atan2d(zequat, sqrt(xequat^2 + yequat^2));
	
	% calculate GMST0 	
	GMST0 = revolve_degree(L + 180) / 15;

	% calculate SIDTIME and Hour Angle
	SIDTIME = sidtime(day_number, longitude, UT);
	HA = (SIDTIME - RA) * 15;

	% convert HA and Decl to rectangular system
	x2 = cosd(HA) * cosd(Decl);
	y2 = sind(HA) * cosd(Decl);
	z2 = sind(Decl);

	% rotate this along the y2 axis
	xhor = x2 * sind(latitude) - z2 * cosd(latitude);
	yhor = y2;
	zhor = x2 * cosd(latitude) + z2 * sind(latitude);

	% finally calculate azimuth and altitude 
	azimuth = atan2d(yhor, xhor) + 180;
	altitude = atan2d(zhor, sqrt(xhor^2 + yhor^2));
end
