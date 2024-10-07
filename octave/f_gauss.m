function y = f_gauss(x, x0, sigma)
	% calculate normalized Gaussian distribution
	
	y = exp(-(x-x0).^2/(2*sigma^2))/(sigma*sqrt(2*pi));
end