function y = f_rect(x, x0, L)
	% calculate normalized rectangular function centered at x0 with width L (with boundary values 1/2L)
	
	% define default center position
	if (nargin < 2)
		x0 = 0;
	end
	% define default peak width
	if (nargin < 3)
		L = 1;
	end

	y = (f_step(x - x0 + L / 2) - f_step(x - x0 - L / 2)) ./ L;
end