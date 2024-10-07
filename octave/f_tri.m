function y = f_tri(x, x0, L)
	% calculate normalized triangular function centered at x0 with width L
	
	% define default center position
	if (nargin < 2)
		x0 = 0;
	end
	% define default peak width
	if (nargin < 3)
		L = 1;
	end

	y = f_rect(x, x0, L) .* (L ./ 2.0 - abs(x - x0)) .* 4.0 ./ L;
end