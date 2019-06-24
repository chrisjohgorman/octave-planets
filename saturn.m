function [latitude, longitude, r] = saturn(day_number)
	N = 113.6634 + 2.38980E-5   * day_number; % Long of asc. node
	i =   2.4886 - 1.081E-7     * day_number; % Inclination
	w = 339.3939 + 2.97661E-5   * day_number; % Argument of perihelion
	a = 9.55475;			       	  % Semi-major axis
	e = 0.055546 - 9.499E-9     * day_number; % eccentricity
	M = 316.9670 + 0.0334442282 * day_number; % Mean anomaly
	M = revolve_degree(M);
	Ms = M;
	Mj = 19.8950 + 0.0830853001 * day_number; % Mean anomaly Jupiter
	Mj = revolve_degree(Mj);
	Mu = 142.5905 + 0.011725806 * day_number; % Mean anomaly Uranus
	Mu = revolve_degree(Mu);
	oblecl = 23.4393 - 3.563e-7 * day_number; % obliquity of the eliptic

	E = eccentric_anomaly(M, e, 0.0005);
	% saturn's rectrangular coordinates
	x = a * (cosd(E) - e);
	y = a * sind(E) * sqrt(1 - e*e);
	% convert to distance and true anomaly
	r = sqrt(x*x + y*y);
	v = atan2d(y, x);
	% saturn's position in ecliptic coordinates
	xeclip = r * ( cosd(N) * cosd(v+w) - sind(N) * sind(v+w) * cosd(i));
	yeclip = r * ( sind(N) * cosd(v+w) + cosd(N) * sind(v+w) * cosd(i));
	zeclip = r * sind(v+w) * sind(i);
	% add sun's rectangular coordinates
	[x1, y1, z1, sunoblecl, L] = sun_rectangular(day_number);
	xgeoc = x1 + xeclip;
	ygeoc = y1 + yeclip;
	zgeoc = z1 + yeclip;
	% rotate the equitorial coordinates
	xequat = xgeoc;
	yequat = ygeoc * cosd(oblecl) - zgeoc * sind(oblecl);
	zequat = ygeoc * sind(oblecl) + zgeoc * cosd(oblecl);
	% convert to RA and Decl
	RA = atan2d(yequat, xequat);
	RA = revolve_degree(RA);
	RA = RA / 15;
	Decl = atan2d(zequat, sqrt(xequat*xequat + yequat*yequat));
	R = sqrt(xequat^2+yequat^2+zequat^2);
	% convert to ecliptic longitude and latitude
	longitude = atan2d(yeclip, xeclip);
	longitude = revolve_degree(longitude);
	latitude = atan2d(zeclip, sqrt(xeclip*xeclip + yeclip*yeclip));
	perturbations_in_longitude = 0.812 * sind(2*Mj - 5*Ms - 67.6) ...
				    -0.229 * cosd(2*Mj - 4*Ms - 2) ...
				    +0.119 * sind(Mj - 2*Ms - 3) ...
				    +0.046 * sind(2*Mj - 6*Ms - 69) ... 
				    +0.014 * sind(Mj - 3*Ms + 32); 
	perturbations_in_latitude = -0.020 * cosd(2*Mj - 4*Ms - 2) ...
				    +0.018 * sind(2*Mj - 6*Ms - 49);
	longitude = longitude + perturbations_in_longitude;
	longitude = revolve_degree(longitude);
	latitude = latitude + perturbations_in_latitude;
	latitude = revolve_degree(latitude);
end
