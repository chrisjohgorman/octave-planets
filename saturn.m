function [right_ascension, declination, distance, azimuth, altitude] = saturn(day_number, latitude, longitude, UT)
    N = 113.6634 + 2.38980E-5   * day_number; % Long of asc. node
    i =   2.4886 - 1.081E-7     * day_number; % Inclination
    w = 339.3939 + 2.97661E-5   * day_number; % Argument of perihelion
    a = 9.55475;                      % Semi-major axis
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
    [x1, y1, z1, sunoblecl, L, lonsun, rs] = sun_rectangular(day_number);
    % convert to ecliptic longitude and latitude
    lon = atan2d(yeclip, xeclip);
    lon = revolve_degree(lon);
    lat = atan2d(zeclip, sqrt(xeclip*xeclip + yeclip*yeclip));
    perturbations_in_longitude = 0.812 * sind(2*Mj - 5*Ms - 67.6) ...
                    -0.229 * cosd(2*Mj - 4*Ms - 2) ...
                    +0.119 * sind(Mj - 2*Ms - 3) ...
                    +0.046 * sind(2*Mj - 6*Ms - 69) ... 
                    +0.014 * sind(Mj - 3*Ms + 32); 
    perturbations_in_latitude = -0.020 * cosd(2*Mj - 4*Ms - 2) ...
                    +0.018 * sind(2*Mj - 6*Ms - 49);
    lon = lon + perturbations_in_longitude;
    lat = lat + perturbations_in_latitude;
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
    Dec = atan2d(ze, sqrt(xe*xe+ye*ye));
    rg = sqrt(xe*xe+ye*ye+ze*ze);
    % convert to azimuth and altitude
    right_ascension = revolve_hour_angle(RA/15);
    declination = Dec;
    distance = rg;
    % convert to azimuth and altitude
    hour_angle = sidtime(day_number, longitude, UT) - right_ascension;
    hour_angle = revolve_hour_angle(hour_angle);
    hour_angle = hour_angle * 15;
    x = cosd(hour_angle)*cosd(declination);
    y = sind(hour_angle)*cosd(declination);
    z = sind(declination);
    x_horizon = x * sind(latitude) - z * cosd(latitude);
    y_horizon = y;
    z_horizon = x * cosd(latitude) + z * sind(latitude);
    azimuth = atan2d(y_horizon,x_horizon) + 180;
    altitude = atan2d(z_horizon, sqrt(x_horizon^2 + y_horizon^2));
end
