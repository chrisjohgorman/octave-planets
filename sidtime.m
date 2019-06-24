function [SIDTIME] = sidtime(day_number, longitude, UT)
	[x1, y1, z1, oblecl, L] = sun_rectangular(day_number);
	GMST0 = revolve_degree(L + 180);
	if ~exist('UT', 'var')
		UT = gmtime(time());
		UT = UT.hour + UT.min/60 + UT.sec/3600;
	end	
	SIDTIME = GMST0/15 + UT + longitude/15;
end	
