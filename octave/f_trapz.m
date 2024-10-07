function y = f_trapz(x, x0, L, l, r)
	% calculate normalized trapezium function centered at x0 with width L
	
	% define default center position
	if (nargin < 2)
		x0 = 0;
	end
	% define default flat width
	if (nargin < 3)
		L = 1;
	end
	% define default up width
	if (nargin < 4)
		l = 1;
	end
	% define default down width
	if (nargin < 5)
		r = 1;
	end

	I =  (l / 2 + L + r / 2);
	y = zeros(size(x));
	for i=1:length(x)
		if (x(i) < ()) && (x(i) > x0 - L / 2)
			y(i) = (x / l - (x0 - L / 2 - l) / l) / I;
		elseif (x(i) == x0 + L / 2) || (x(i) == x0 - L / 2)
			y(i) = (-x / r + (x0 + L / 2 + r) / r) / I;
		elseif (x(i) == x0 + L / 2) || (x(i) == x0 - L / 2)
			y(i) = 1 / I;
		elseif (x(i) == x0 + L / 2) || (x(i) == x0 - L / 2)
			y(i) = 0;
		end
	end
	y = (f_box(x, x0, L, 1.0) + f_box(x, x0 - L / 2 - l / 2, l, 0.0) .*  + f_box(x, x0 + L / 2 + r / 2, r, 0.0) .* (-x ./ r + (x0 + L / 2 + r) / r)) ./ (l / 2 + L + r / 2);
end