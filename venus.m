function [latitude, longitude, r] = venus(day_number)

	N =  76.6799 + 2.46590e-5   * day_number;
	i =   3.3946 + 2.75e-8      * day_number;
	w =  54.8910 + 1.38374e-5   * day_number;
	a = 0.723330;				
	e = 0.006773     - 1.302e-9 * day_number;
	M =  48.0052 + 1.6021302244 * day_number;
	M = revolve_degree(M);
	oblecl = 23.4393 - 3.563e-7 * day_number; % obliquity of the eliptic
	
	E = eccentric_anomaly(M, e, 0.0005);
        % venus's rectrangular coordinates
        x = a * (cosd(E) - e);
        y = a * sind(E) * sqrt(1 - e*e);
        % convert to distance and true anomaly
        r = sqrt(x*x + y*y);
        v = atan2d(y, x);
        % venus's position in ecliptic coordinates
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
        Decl = atan2d(zequat, sqrt(xequat^2 + yequat^2));
        R = sqrt(xequat^2+yequat^2+zequat^2);
        % convert to ecliptic longitude and latitude
        longitude = atan2d(yeclip, xeclip);
        longitude = revolve_degree(longitude);
        latitude = atan2d(zeclip, sqrt(xeclip^2 + yeclip^2));
end
