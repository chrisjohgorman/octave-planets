function [marsrise, marsset] = marsrise_set (year, month, day, hour, latitude, longitude)
%
%First we must decide which altitude we're interested in:
%
% h = 0 degrees: Center of Sun's disk touches a mathematical horizon
% h = -0.25 degrees: Sun's upper limb touches a mathematical horizon
% h = -0.583 degrees: Center of Sun's disk touches the horizon; atmospheric refraction accounted for
% h = -0.833 degrees: Sun's upper limb touches the horizon; atmospheric refraction accounted for
% h = -6 degrees: Civil twilight (one can no longer read outside without artificial illumination)
% h = -12 degrees: Nautical twilight (navigation using a sea horizon no longer possible)
% h = -15 degrees: Amateur astronomical twilight (the sky is dark enough for most astronomical observations)
% h = -18 degrees: Astronomical twilight (the sky is completely dark)
%
    h = -0.833;
    d = day_number(year, month, day, hour);
    [x1, y1, z1, oblecl, L] = sun_rectangular(d);
    [RA, Decl, r, az, alt] = mars(d, latitude, longitude, hour);
    GMST0 = (L + 180) / 15;
    UT_Planet_in_south = RA - (L+180)/15 - longitude/15.0;
    UT_Planet_in_south = revolve_hour_angle(UT_Planet_in_south);
    cos_lha = (sind(h) - sind(latitude)*sind(Decl))/(cosd(latitude) * cosd(Decl));
    if (cos_lha > 1)
        error("Mars is always below our altitude limit.");
    elseif (cos_lha < -1)
        error("Mars is always above our altitude limit.");
    end
    LHA = acosd(cos_lha)/15.04107;
    time = localtime(time);
    marsrise = UT_Planet_in_south - LHA;
    for i=1:5
        [marsrise] = mariset(year,month,day,marsrise,latitude,longitude);
    end
    marsrise = marsrise + time.gmtoff/3600;
    marsrise = datestr((marsrise/24), 'HH:MM');
    marsset = UT_Planet_in_south + LHA;
    for i=1:5
        [mrise, marsset] = mariset(year,month,day,marsset,latitude,longitude);
    end
    marsset = marsset + time.gmtoff/3600;
    marsset = datestr((marsset/24), 'HH:MM');
end
