function [right_ascension, declination, distance, azimuth, altitude] = uranus(day_number, latitude, longitude)
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
    zgeoc = z1 + zeclip;
    % rotate the equitorial coordinates
    xequat = xgeoc;
    yequat = ygeoc * cosd(oblecl) - zgeoc * sind(oblecl);
    zequat = ygeoc * sind(oblecl) + zgeoc * cosd(oblecl);
    % convert to right_ascension and declination
    right_ascension = atan2d(yequat, xequat);
    right_ascension = revolve_degree(right_ascension);
    right_ascension = right_ascension / 15;
    declination = atan2d(zequat, sqrt(xequat^2 + yequat^2));
    distance = sqrt(xequat^2+yequat^2+zequat^2);
    % convert to ecliptic longitude and latitude
    lon = atan2d(yeclip, xeclip);
    lon = revolve_degree(lon);
    lat = atan2d(zeclip, sqrt(xeclip^2 + yeclip^2));
    perturbations_in_longitude = +0.040 * sind(Ms - 2*Mu + 6) ...
                        +0.035 * sind(Ms - 3*Mu + 33) ...
                        -0.015 * sind(Mj - Mu + 20);
    lon = perturbations_in_longitude + lon;
    lon = revolve_degree(lon);
    % convert to azimuth and altitude
    hour_angle = sidtime(day_number, longitude, 0) - right_ascension;
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
