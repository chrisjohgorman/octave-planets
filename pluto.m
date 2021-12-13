function [RA, Dec, rg, azimuth, altitude] = pluto(day_number, latitude, longitude, UT)
    S  =   50.03  +  0.033459652 * day_number;
    P  =  238.95  +  0.003968789 * day_number;
    lonecl = 238.9508  +  0.00400703 * day_number ...
            - 19.799 * sind(P)     + 19.848 * cosd(P) ...
             + 0.897 * sind(2*P)    - 4.956 * cosd(2*P) ...
             + 0.610 * sind(3*P)    + 1.211 * cosd(3*P) ...
             - 0.341 * sind(4*P)    - 0.190 * cosd(4*P) ...
             + 0.128 * sind(5*P)    - 0.034 * cosd(5*P) ...
             - 0.038 * sind(6*P)    + 0.031 * cosd(6*P) ...
             + 0.020 * sind(S-P)    - 0.010 * cosd(S-P);
    latecl =  -3.9082 ...
             - 5.453 * sind(P)     - 14.975 * cosd(P) ...
             + 3.527 * sind(2*P)    + 1.673 * cosd(2*P) ...
             - 1.051 * sind(3*P)    + 0.328 * cosd(3*P) ...
             + 0.179 * sind(4*P)    - 0.292 * cosd(4*P) ...
             + 0.019 * sind(5*P)    + 0.100 * cosd(5*P) ...
             - 0.031 * sind(6*P)    - 0.026 * cosd(6*P) ...
                                   + 0.011 * cosd(S-P);
    r     =  40.72 ...
           + 6.68 * sind(P)       + 6.90 * cosd(P) ...
           - 1.18 * sind(2*P)     - 0.03 * cosd(2*P) ...
           + 0.15 * sind(3*P)     - 0.14 * cosd(3*P);

    xh = r * cosd(lonecl) * cosd(latecl);
    yh = r * sind(lonecl) * cosd(latecl);
    zh = r               * sind(latecl);

    [xs, ys, zs, oblecl, L, lonsun, rs] = sun_rectangular(day_number); 

    xg = xh + xs;
    yg = yh + ys;
    zg = zh;

    xe = xg;
    ye = yg * cosd(oblecl) - zg * sind(oblecl);
    ze = yg * sind(oblecl) + zg * cosd(oblecl);

    RA  = atan2d( ye, xe );
    RA  = revolve_degree(RA);
    RA = RA/15; 
    Dec = atan2d( ze, sqrt(xe*xe+ye*ye) );
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
