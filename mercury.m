function [right_ascension, declination, distance, azimuth, altitude] = mercury(day_number, latitude, longitude, UT)
    N =  48.3313 + 3.24587e-5   * day_number;   % (Long of asc. node)
    i =   7.0047 + 5.00e-8      * day_number;   % (Inclination)
    w =  29.1241 + 1.01444e-5   * day_number;   % (Argument of perihelion)
    a = 0.387098;                               % (Semi-major axis)
    e = 0.205635 + 5.59e-10     * day_number;   % (Eccentricity)
    M = 168.6562 + 4.0923344368 * day_number;   % (Mean anonaly)
    M = revolve_degree(M);
    oblecl = 23.4393 - 3.563e-7 * day_number;   % obliquity of the eliptic

    E = eccentric_anomaly(M, e, 0.0005);
    % mercury's rectrangular coordinates
    x = a * (cosd(E) - e);
    y = a * sind(E) * sqrt(1 - e*e);
    % convert to distance and true anomaly
    r = sqrt(x*x + y*y);
    v = atan2d(y, x);
    % mercury's position in ecliptic coordinates (heliocentric)
    xeclip = r * ( cosd(N) * cosd(v+w) - sind(N) * sind(v+w) * cosd(i));
    yeclip = r * ( sind(N) * cosd(v+w) + cosd(N) * sind(v+w) * cosd(i));
    zeclip = r * sind(v+w) * sind(i);

    % add sun's rectangular coordinates
    [x1, y1, z1, sunoblecl, L, lonsun, rs] = sun_rectangular(day_number);
    xgeoc = x1 + xeclip;
    ygeoc = y1 + yeclip;
    zgeoc = z1 + zeclip;
    % rotate the equitorial coordinates
    xequat = xgeoc;
    yequat = ygeoc * cosd(oblecl) - zgeoc * sind(oblecl);
    zequat = ygeoc * sind(oblecl) + zgeoc * cosd(oblecl);
    % convert to right_ascension and declination
    right_ascension = atan2d(yequat, xequat);
    right_ascension = revolve_degree(right_ascension);
    right_ascension = right_ascension/15;
    declination = atan2d(zequat, sqrt(xequat^2 + yequat^2));
    distance = sqrt(xequat^2+yequat^2+zequat^2);
    % convert to ecliptic longitude and latitude
    lon = atan2d(yeclip, xeclip);
    lon = revolve_degree(lon);
    lat = atan2d(zeclip, sqrt(xeclip*xeclip + yeclip*yeclip));
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
