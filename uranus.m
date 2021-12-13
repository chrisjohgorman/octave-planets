function [RA, Dec, rg, azimuth, altitude] = uranus(day_number, latitude, longitude, UT)
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
    % sun's rectangular coordinates
    [x1, y1, z1, sunoblecl, L, lonsun, rs] = sun_rectangular(day_number);
    % convert to ecliptic longitude and latitude
    lon = atan2d(yeclip, xeclip);
    lon = revolve_degree(lon);
    lat = atan2d(zeclip, sqrt(xeclip^2 + yeclip^2));
    perturbations_in_longitude = +0.040 * sind(Ms - 2*Mu + 6) ...
                        +0.035 * sind(Ms - 3*Mu + 33) ...
                        -0.015 * sind(Mj - Mu + 20);
    lon = perturbations_in_longitude + lon;
    lon = revolve_degree(lon);
    % use lat and lon to produce new RA and Dec
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
