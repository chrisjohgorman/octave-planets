function [new_degree] = revolve_degree(degree)
	new_degree = degree - floor(degree/360)*360;
end
