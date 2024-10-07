function y = f_box(x, x0, L, delimiter)
	% calculate rectangular function centered at x0 with width L and height 1.0
	
	% define default center position
	if (nargin < 2)
		x0 = 0;
	end
	% define default peak width
	if (nargin < 3)
		L = 1;
	end

	y = zeros(size(x));
	for i=1:length(x)
		if (x(i) < x0 + L / 2) && (x(i) > x0 - L / 2)
			y(i) = 1;
		elseif (x(i) == x0 + L / 2) || (x(i) == x0 - L / 2)
			y(i) = delimiter;
		else
			y(i) = 0;
		end
	end
end