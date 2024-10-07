function I = beam_profile(theta, E, method, p)
	
	% select default method
	if (nargin < 3)
		method = 'rect';
	end
	
	switch method
		case 'expt'
			beam = load('beam.dat');
			I = normalize(beam(:,E));
		case 'rect'
			% f_rect(x, x0, L)
			if (nargin < 4)
				p = [0.0, 1.0];
			end
			I = f_rect(theta, p(1), p(2));
		case 'tri'
			% f_tri(x, x0, L)
			if (nargin < 4)
				p = [0.0, 1.0];
			end
			I = f_tri(theta, p(1), p(2));
		case 'gauss'
			% f_gauss(x, x0, sigma)
			if (nargin < 4)
				p = [0.0, 1.0];
			end
			I = f_gauss(theta, p(1), (p(2))/sqrt(8.0*log(2.0)));
		case 'cauchy'
			% f_cauchy(x, x0, hwhm)
			if (nargin < 4)
				p = [0.0, 1.0];
			end
			I = f_cauchy(theta, p(1), p(2));
		case 'trapz'
			% f_trapz(x, x0, l_flat, l_up, l_down)
			if (nargin < 4)
				p = [0.0, 1.0, 1.0, 1.0];
			end
			I = f_trapz(theta, p(1), p(2), p(3), p(4));
		otherwise
			I = zeros(sizeof(theta));
	end
end