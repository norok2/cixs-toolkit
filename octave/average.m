function [avg, err, stdev_avg, stdev] = average(v, w)
	% computes the weighted average of vector v using weights w and the standard deviation of the average (i.e. the statistical error on average assuming Gaussian distribution)
	
	if (nargin < 2)
		w = ones(size(v));
	end
	
	avg = sum(v .* w) / sum(w);
	stdev = sqrt(sum((v .* w).^2.0) / sum(w) - avg.^2.0);
	stdev_avg = stdev / sqrt(numel(v));
	err = 3.0 * stdev_avg;
end