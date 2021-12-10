function [right_ascension, declination, distance, azimuth, altitude] = sun(day_number, latitude, longitude, UT)
    % rotate equitorial coordinates
    [x1, y1, z1, oblecl, L, lon, r] = sun_rectangular(day_number);
    xequat = x1;
    yequat = y1 * cosd(oblecl) - z1 * sind(oblecl);
    zequat = y1 * sind(oblecl) - z1 * cosd(oblecl);

    right_ascension = atan2d(yequat, xequat);
    right_ascension = revolve_degree(right_ascension);
    % convert right_ascension to hours
    right_ascension = right_ascension / 15;
    declination = atan2d(zequat, sqrt(xequat^2 + yequat^2));
    distance = sqrt(xequat^2 + yequat^2 + zequat^2);
    
    % calculate hour angle
    hour_angle = (sidtime(day_number, longitude, UT) - right_ascension) * 15;

    % convert hour_angle and declination to rectangular system
    x2 = cosd(hour_angle) * cosd(declination);
    y2 = sind(hour_angle) * cosd(declination);
    z2 = sind(declination);

    % rotate this along the y2 axis
    xhor = x2 * sind(latitude) - z2 * cosd(latitude);
    yhor = y2;
    zhor = x2 * cosd(latitude) + z2 * sind(latitude);

    % finally calculate azimuth and altitude 
    azimuth = atan2d(yhor, xhor) + 180;
    altitude = atan2d(zhor, sqrt(xhor^2 + yhor^2));
end
