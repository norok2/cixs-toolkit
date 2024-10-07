function y = f_lorentz(x, x0, gamma, I)
	% calculate Lorentz distribution: non-normalized Cauchy distribution, with maximum value I specified
	
	y = I .* gamma.^2.0 ./ ((x - x0).^2.0 + gamma.^2.0);
end