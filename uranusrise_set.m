function [uranusrise, uranusset] = uranusrise_set (year, month, day, latitude, longitude)
    h = -0.833;
    d = day_number(year, month, day);
    [x1, y1, z1, oblecl, L] = sun_rectangular(d);
    [RA, Decl, r, az, alt] = uranus(d, latitude, longitude);
    GMST0 = (L + 180) / 15;
    UT_Planet_in_south = RA - (L+180)/15 - longitude/15.0;
    UT_Planet_in_south = revolve_hour_angle(UT_Planet_in_south);
    cos_lha = (sind(h) - sind(latitude)*sind(Decl))/(cosd(latitude) * cosd(Decl));
    if (cos_lha > 1)
        error("Uranus is always below our altitude limit.");
    elseif (cos_lha < -1)
        error("Uranus is always above our altitude limit.");
    end
    LHA = acosd(cos_lha)/15.04107;
    time = localtime(time);
    ur = UT_Planet_in_south - LHA + time.gmtoff/3600;
    us = UT_Planet_in_south + LHA + time.gmtoff/3600;
    uranusrise = datestr(ur/24, 'HH:MM');
    uranusset = datestr(us/24, 'HH:MM');
end
