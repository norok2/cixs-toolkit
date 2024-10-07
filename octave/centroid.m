function [x0] = centroid(x, y)
	% calculate centroid of a function

	x0 = sum(x .* y) / sum(y);
end
