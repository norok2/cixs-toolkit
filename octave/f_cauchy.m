function y = f_cauchy(x, x0, gamma)
	% calculate normalized Cauchy distribution
	
	y = gamma ./ pi ./ ((x - x0).^2.0 + gamma.^2.0);
end