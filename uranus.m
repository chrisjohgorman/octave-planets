function [latitude, longitude, r] = uranus(day_number)

	N =  74.0005 + 1.3978E-5    * day_number;  % Long of asc. node
	i =   0.7733 + 1.9E-8       * day_number;  % Inclination
	w =  96.6612 + 3.0565E-5    * day_number;  % Argument of perihelion
	a = 19.18171 - 1.55E-8      * day_number;  % Semi-major axis
	e = 0.047318 + 7.45E-9      * day_number;  % eccentricity
	M = 142.5905 + 0.011725806  * day_number;  % Mean anomaly Uranus
	M = revolve_degree(M);
	Mu = M;
	Ms = 316.9670 + 0.0334442282 * day_number; % Mean anomaly Saturn
	Ms = revolve_degree(Ms);
	Mj =  19.8950 + 0.0830853001 * day_number; % Mean anomaly Jupiter
	Mj = revolve_degree(Mj);
	oblecl = 23.4393 - 3.563e-7 * day_number;  % obliquity of the eliptic

	E = eccentric_anomaly(M, e, 0.0005);
        % uranus's rectrangular coordinates
        x = a * (cosd(E) - e);
        y = a * sind(E) * sqrt(1 - e*e);
        % convert to distance and true anomaly
        r = sqrt(x*x + y*y);
        v = atan2d(y, x);
        % uranus's position in ecliptic coordinates
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
        Decl = atan2d(zequat, sqrt(xequat^2 + yequat^2));
        R = sqrt(xequat^2+yequat^2+zequat^2);
        % convert to ecliptic longitude and latitude
        longitude = atan2d(yeclip, xeclip);
        longitude = revolve_degree(longitude);
        latitude = atan2d(zeclip, sqrt(xeclip^2 + yeclip^2)) 
	perturbations_in_longitude = +0.040 * sind(Ms - 2*Mu + 6) ...
			             +0.035 * sind(Ms - 3*Mu + 33) ...
                              	     -0.015 * sind(Mj - Mu + 20);
	longitude = perturbations_in_longitude + longitude;
	longitude = revolve_degree(longitude);
end
