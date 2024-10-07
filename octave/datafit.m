function [pp, yy, chi2] = datafit(x, y, p0, func, method)
	% fit data (x,y) with func(x, p) starting from initial p0 parameters using fminsearch

	if (nargin < 5)
		method = 'fminsearch';
	end

	% remove NaN from data
	y_nan = isnan(y);
	x(y_nan) = [];
	y(y_nan) = [];

	switch method
		case 'fminsearch' % Nelder-Mead Simplex (Matlab / Octave)
			pp = fminsearch(@(p) (sum((y - func(x, p)).^2.0)), p0);
			yy = func(x, pp);
			chi2 = sum((y - func(x, pp)).^2.0 ./ y.^2.0) ./ (length(y) - length(p0) - 1);
		case 'leasqr' % Levenberg-Marquardt (Octave)
			[yy, pp] = leasqr(x, y, p0, func);
			%[yy, pp] = leasqr(x, y, p0, func, 1e-6, 20, sqrt(abs(y)));
		case 'lsqcurvefit' % Levenberg-Marquardt (Matlab)
			pp = lsqcurvefit(@(p, x) (func(x, p)), p0, x, y);
		otherwise
			error('Unknown fitting algorithm');
	end
	yy = func(x, pp);
	chi2 = sum((y - func(x, pp)).^2.0 ./ y.^2.0) ./ (length(y) - length(p0) - 1);
end

