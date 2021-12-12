function [RA, Dec, rg, azimuth, altitude] = mars(day_number, latitude, longitude, UT)
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
    HA = revolve_hour_angle(HA);
    HA = HA * 15;
    x = cosd(HA)*cosd(Dec);
    y = sind(HA)*cosd(Dec);
    z = sind(Dec);
    xhor = x * sind(latitude) - z * cosd(latitude);
    yhor = y;
    zhor = x * cosd(latitude) + z * sind(latitude);
    azimuth = atan2d(yhor,xhor) + 180;
    altitude = atan2d(zhor, sqrt(xhor^2 + yhor^2));
end
