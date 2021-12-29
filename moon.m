function [tRA, tDec, r, azimuth, alt_topoc] =  moon(day_number, latitude, longitude, UT)
    N = 125.1228 - 0.0529538083 * day_number;  % long asc. node
    i = 5.1454;                % inclination
    w = 318.0634 + 0.1643573223 * day_number;  % Arg. of perigree
    a = 60.2666;                   % mean distance
    e = 0.054900;                  % eccentricity
    M = 115.3654 + 13.0649929509 * day_number; % mean anomaly
    M = revolve_degree(M);
    oblecl = 23.4393 - 3.563e-7 * day_number;  % obliquity of the eliptic
    E = eccentric_anomaly(M, e, 0.0005);
    % moon's rectrangular coordinates
    x = a * (cosd(E) - e);
    y = a * sind(E) * sqrt(1 - e*e);
    % convert to distance and true anomaly
    r = sqrt(x*x + y*y);
    v = atan2d(y, x);
    % moon's position in ecliptic coordinates
    xeclip = r * ( cosd(N) * cosd(v+w) - sind(N) * sind(v+w) * cosd(i));
    yeclip = r * ( sind(N) * cosd(v+w) + cosd(N) * sind(v+w) * cosd(i));
    zeclip = r * sind(v+w) * sind(i);
    % convert to ecliptic longitude, latitude and distance
    lon = atan2d(yeclip, xeclip);
    lon = revolve_degree(lon);
    lat = atan2d(zeclip, sqrt(xeclip*xeclip + yeclip*yeclip));

    Sw = 282.9404 + 4.70935e-5   * day_number; % sun's (longitude of 
                           % perihelion)
    Ms = 356.0470 + 0.9856002585 * day_number; % sun's mean anomaly
    Ls = Sw + Ms;
    Ls = revolve_degree(Ls);
    Lm = N + w + M;
    Mm = M;
    D = Lm - Ls;
    F = Lm - N;

    perturbations_in_longitude = -1.274 * sind(Mm - 2*D) ...
                     +0.658 * sind(2*D) ... 
                     -0.186 * sind(Ms) ...
                     -0.059 * sind(2*Mm - 2*D) ...
                     -0.057 * sind(Mm - 2*D + Ms) ... 
                                     +0.053 * sind(Mm + 2*D) ...
                     +0.046 * sind(2*D - Ms) ... 
                     +0.041 * sind(Mm - Ms) ...
                     -0.035 * sind(D) ...
                     -0.031 * sind(Mm + Ms) ...
                     -0.015 * sind(2*F - 2*D) ... 
                     +0.011 * sind(Mm - 4*D);
    perturbations_in_latitude = -0.173 * sind(F - 2*D) ...
                    -0.055 * sind(Mm - F - 2*D) ...
                    -0.046 * sind(Mm + F - 2*D) ...
                    +0.033 * sind(F + 2*D) ...
                    +0.017 * sind(2*Mm + F);
    perturbations_in_distance = -0.58 * cosd(Mm - 2*D) ...
                        -0.46 * cosd(2*D);
    lon = lon + perturbations_in_longitude;
    lat = lat + perturbations_in_latitude;
    r = r + perturbations_in_distance;
    x1 = cosd(lon)*cosd(lat);
    y1 = cosd(lat)*sind(lon);
    z1 = sind(lat);
    x2 = x1;
    y2 = y1 * cosd(oblecl) - z1 * sind(oblecl);
    z2 = y1 * sind(oblecl) + z1 * cosd(oblecl);
    RA = atan2d(y2,x2);
    RA = revolve_degree(RA);
    Dec = atan2d(z2, sqrt(x2*x2 + y2*y2));
    HA = (sidtime(day_number, longitude, UT)-RA / 15) *15 ;
    HA = revolve_degree(HA);
    x = cosd(HA) * cosd(Dec);
    y = sind(HA) * cosd(Dec);
    z = sind(Dec);
    xhor = x * sind(latitude) - z * cosd(latitude);
    yhor = y;
    zhor = x * cosd(latitude) + z * sind(latitude);
    azimuth  = atan2d( yhor, xhor ) + 180;
    altitude = atan2d( zhor , sqrt(xhor^2 + yhor^2));
    mpar = asind(1/r);
    alt_topoc = altitude - mpar * cosd(altitude);
    gclat = latitude - 0.1924 * sind(2*latitude);
    rho   = 0.99833 + 0.00167 * cosd(2*latitude);
    g = atand(tand(gclat) / cosd(HA));
    tRA   = RA ... 
        - mpar * rho * cosd(gclat) * sind(HA) / cosd(Dec);
    tRA = revolve_degree(tRA);
    tDec = Dec ...
        - mpar * rho * sind(gclat) * sind(g - Dec) / sind(g);
end
