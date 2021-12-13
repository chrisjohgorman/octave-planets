function [RA, Dec, rg, azimuth, altitude] = neptune(day_number, latitude, longitude, UT)
    N = 131.7806 + 3.0173E-5    * day_number;
    i =   1.7700 - 2.55E-7      * day_number;
    w = 272.8461 - 6.027E-6     * day_number;
    a = 30.05826 + 3.313E-8     * day_number;
    e = 0.008606 + 2.15E-9      * day_number;
    M = 260.2471 + 0.005995147  * day_number;
    M = revolve_degree(M);
    oblecl = 23.4393 - 3.563e-7 * day_number; % obliquity of the eliptic
    
    E = eccentric_anomaly(M, e, 0.0005);
    % neptune's rectrangular coordinates
    x = a * (cosd(E) - e);
    y = a * sind(E) * sqrt(1 - e*e);
    % convert to distance and true anomaly
    r = sqrt(x*x + y*y);
    v = atan2d(y, x);
    % neptune's position in ecliptic coordinates
    xeclip = r * ( cosd(N) * cosd(v+w) - sind(N) * sind(v+w) * cosd(i));
    yeclip = r * ( sind(N) * cosd(v+w) + cosd(N) * sind(v+w) * cosd(i));
    zeclip = r * sind(v+w) * sind(i);
    % sun's rectangular coordinates
    [x1, y1, z1, sunoblecl, L, lonsun, rs] = sun_rectangular(day_number);
    % convert to ecliptic longitude and latitude
    lon = atan2d(yeclip, xeclip);
    lon = revolve_degree(lon);
    lat = atan2d(zeclip, sqrt(xeclip^2 + yeclip^2));
    % use lat and lon to produce RA and Dec
    xh = r * cosd(lon) * cosd(lat);
    yh = r * sind(lon) * cosd(lat);
    zh = r * sind(lat);
    % convert sun's position
    xs = rs * cosd(lonsun);
    ys = rs * sind(lonsun);
    % convert from heliocentric to geocentric
    xg = xh + xs;
    yg = yh + ys;
    zg = zh;
    % convert to equitorial
    xe = xg;
    ye = yg * cosd(oblecl) - zg * sind(oblecl);
    ze = yg * sind(oblecl) + zg * cosd(oblecl);
    % RA and Decl
    RA = atan2d(ye, xe);
    RA = revolve_degree(RA);
    RA = RA/15;
    Dec = atan2d(ze, sqrt(xe*xe+ye*ye));
    rg = sqrt(xe*xe+ye*ye+ze*ze);
    % convert to azimuth and altitude
    HA = sidtime(day_number, longitude, UT) - RA;
    HA = HA * 15;
    HA = revolve_degree(HA);
    x = cosd(HA)*cosd(Dec);
    y = sind(HA)*cosd(Dec);
    z = sind(Dec);
    xhor = x * sind(latitude) - z * cosd(latitude);
    yhor = y;
    zhor = x * cosd(latitude) + z * sind(latitude);
    azimuth = atan2d(yhor,xhor) + 180;
    altitude = atan2d(zhor, sqrt(xhor^2 + yhor^2));
end
