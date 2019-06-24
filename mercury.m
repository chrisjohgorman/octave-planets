function [latitude, longitude, r] = mercury(day_number)
%
%
	N =  48.3313 + 3.24587e-5   * day_number;   % (Long of asc. node)
	i =   7.0047 + 5.00e-8      * day_number;   % (Inclination)
	w =  29.1241 + 1.01444e-5   * day_number;   % (Argument of perihelion)
	a = 0.387098;                               % (Semi-major axis)
	e = 0.205635 + 5.59e-10     * day_number;   % (Eccentricity)
	M = 168.6562 + 4.0923344368 * day_number;   % (Mean anonaly)
	M = revolve_degree(M);
	oblecl = 23.4393 - 3.563e-7 * day_number; % obliquity of the eliptic

	E = eccentric_anomaly(M, e, 0.0005);
        % mercury's rectrangular coordinates
        x = a * (cosd(E) - e);
        y = a * sind(E) * sqrt(1 - e*e);
        % convert to distance and true anomaly
        r = sqrt(x*x + y*y);
        v = atan2d(y, x);
        % mercury's position in ecliptic coordinates
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
	RA = RA/15;
	Decl = atan2d(zequat, sqrt(xequat*xequat + yequat*yequat));
	R = sqrt(xequat^2+yequat^2+zequat^2);
        % convert to ecliptic longitude and latitude
        longitude = atan2d(yeclip, xeclip); 
        longitude = revolve_degree(longitude);
        latitude = atan2d(zeclip, sqrt(xeclip*xeclip + yeclip*yeclip));
end
