function [right_ascension, declination, distance, azimuth, altitude] = venus(day_number, latitude, longitude, UT)
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
    %declination = 22 + 48/60 + 32/3600;
    distance = sqrt(xequat^2+yequat^2+zequat^2);
    % convert to ecliptic longitude and latitude
    lon = atan2d(yeclip, xeclip);
    lon = revolve_degree(lon);
    lat = atan2d(zeclip, sqrt(xeclip^2 + yeclip^2));
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
