function [x1, y1, z1, oblecl, L] = sun_rectangular(day_number)
	w = 282.9404 + 4.70935e-5 * day_number;   % longitude of perihelion
	a = 1;                                    % mean distance, a.u.
	e = 0.016709 - 1.151e-9 * day_number;     % eccentricity
	M = 356.0470 + 0.9856002585 * day_number; % mean anomaly
	M = revolve_degree(M);
	oblecl = 23.4393 - 3.563e-7 * day_number; % obliquity of the eliptic
	L = w + M;				  % sun's mean longitude
	L = revolve_degree(L);
	% sun's eccentric anomaly
	E = M + (180/pi) * e * sind(M) * (1 + e * cosd(M));
	% sun's rectrangular coordinates
	x = cosd(E) - e;
	y = sind(E) * sqrt(1 - e*e);
	% convert to distance and true anomaly
	r = sqrt(x*x + y*y);
	v = atan2d(y, x);
	% sun's longitude
	lon = v + w;
	lon = revolve_degree(lon);
	% sun's ecliptic rectangular coordinates
	x1 = r * cosd(lon);
	y1 = r * sind(lon);
	z1 = 0;
end
