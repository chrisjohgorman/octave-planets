function [d] = day_number(year, month, day, UT)
	% working for 0 utc, valid from march 1900 to February 2100.
	%d = 367*year - floor(7 * ( year + floor((month+9)/12))  / 4) + floor(275*month/9) + day - 730530;
	% valid for entire Gregorian Calendar.
	d = 367*year - floor(7 * (year + floor((month+9)/12))  / 4) - floor(3 * ((year + floor((month -9)/ 7)) / 100 + 1) / 4)+ floor(275*month/9) + day - 730515;
	if ~exist('UT', 'var')
		UT = gmtime(time());
		d = d + (UT.hour + UT.min/60 + UT.sec/3600)/24;
	else
		d = d + UT/24;
	end
end
