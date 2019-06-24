function [latitude, longitude, r] = mars(day_number)
    	N =  49.5574 + 2.11081e-5   * day_number;
    	i =   1.8497 - 1.78e-8      * day_number;
    	w = 286.5016 + 2.92961e-5   * day_number;
    	a = 1.523688;				
    	e = 0.093405     + 2.516e-9 * day_number;
    	M =  18.6021 + 0.5240207766 * day_number;
	M = revolve_degree(M);
        oblecl = 23.4393 - 3.563e-7 * day_number; % obliquity of the eliptic

	E = eccentric_anomaly(M, e, 0.0005);
        % mars's rectrangular coordinates
        x = a * (cosd(E) - e);
        y = a * sind(E) * sqrt(1 - e*e);
        % convert to distance and true anomaly
        r = sqrt(x*x + y*y);
        v = atan2d(y, x);
        % mars's position in ecliptic coordinates
        xeclip = r * ( cosd(N) * cosd(v+w) - sind(N) * sind(v+w) * cosd(i));
        yeclip = r * ( sind(N) * cosd(v+w) + cosd(N) * sind(v+w) * cosd(i));
        zeclip = r * sind(v+w) * sind(i);
        % add sun's rectangular coordinates
	[x1, y1, z1, sunoblecl, L] = sun_rectangular(day_number);
	xgeoc = x1 + xeclip;
	ygeoc = y1 + yeclip;
	zgeoc = z1 + zeclip;
        % rotate the equitorial coordinates
        xequat = xgeoc;
        yequat = ygeoc * cosd(oblecl) - zgeoc * sind(oblecl);
        zequat = ygeoc * sind(oblecl) + zgeoc * cosd(oblecl);
        % convert to RA and Decl
        RA = atan2d(yequat, xequat);
        RA = revolve_degree(RA)';
	RA = RA / 15;
        Decl = atan2d(zequat, sqrt(xequat*xequat + yequat*yequat));
        % convert to ecliptic longitude and latitude
        longitude = atan2d(yeclip, xeclip);
        longitude = revolve_degree(longitude);
        latitude = atan2d(zeclip, sqrt(xeclip*xeclip + yeclip*yeclip));
end
