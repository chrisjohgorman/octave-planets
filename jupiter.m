function [right_ascension, declination, distance, azimuth, altitude] = jupiter(day_number, latitude, longitude, UT)
    N = 100.4542 + 2.76854e-5   * day_number;  % Long of asc. node
    i =   1.3030 - 1.557e-7     * day_number;  % Inclination
    w = 273.8777 + 1.64505e-5   * day_number;  % Argument of perihelion
    a = 5.20256;                   % Semi-major axis
    e = 0.048498 + 4.469e-9     * day_number;  % eccentricity
    M =  19.8950 + 0.0830853001 * day_number;  % Mean anomaly Jupiter
    M = revolve_degree(M);
    Mj = M;
    Ms = 316.9670 + 0.0334442282 * day_number; % Mean anomaly Saturn
    Ms = revolve_degree(Ms);
    Mu = 142.5905 + 0.011725806 * day_number;  % Mean anomaly Uranus
    Mu = revolve_degree(Mu);
    oblecl = 23.4393 - 3.563e-7 * day_number;  % obliquity of the eliptic
    
    E = eccentric_anomaly(M, e, 0.0005);
    % jupiter's rectrangular coordinates
    x = a * (cosd(E) - e);
    y = a * sind(E) * sqrt(1 - e*e);
    % convert to distance and true anomaly
    r = sqrt(x*x + y*y);
    v = atan2d(y, x);
    % jupiter's position in ecliptic coordinates
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
    right_ascension = right_ascension/15;
    declination = atan2d(zequat, sqrt(xequat^2 + yequat^2));
    distance = sqrt(xequat^2+yequat^2+zequat^2);
    % convert to ecliptic longitude and latitude
    lon = atan2d(yeclip, xeclip);
    lon = revolve_degree(lon);
    lat = atan2d(zeclip, sqrt(xeclip^2 + yeclip^2));
    perturbations_of_longitude = -0.332 * sind(2*Mj - 5*Ms - 67.6) ...
                     -0.056 * sind(2*Mj - 2*Ms + 21) ...
                     +0.042 * sind(3*Mj - 5*Ms + 21) ...
                     -0.036 * sind(Mj - 2*Ms) ...
                     +0.022 * cosd(Mj - Ms) ...
                     +0.023 * sind(2*Mj - 3*Ms + 52) ...
                     -0.016 * sind(Mj - 5*Ms - 69);
    lon = lon + perturbations_of_longitude;
    lon = revolve_degree(lon);
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
