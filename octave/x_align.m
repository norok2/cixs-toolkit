function dx = x_align(x1, y1, x2, y2, method)
	% search best translation by which f2(x+dx,y) is most similar to f1(x,y) according to specified method assuming f1 is oversampled respect to f2
		
	% select default method
	if (nargin < 5)
		method = 'fminsearch';
	end
	
	% sort functions in ascending order
	[x1, sort_index1] = sort(x1);
	y1 = y1(sort_index1);
	[x2, sort_index2] = sort(x2);
	y2 = y2(sort_index2);
	
	switch method
		case 'centroid'
			dx = centroid(x1, y1) - centroid(x2, y2);
		case 'fminsearch'
			dx = centroid(x1, y1) - centroid(x2, y2);
			dx = fminsearch(@(k) (diff_sq_sum(x1, y1, x2 + k, y2)), dx);
		case 'sequential'
			x0 = centroid(x1, y1) - centroid(x2, y2);
			range = min(abs([x2(end) - x2(1), x1(end) - x1(1)]));
			step = abs(x1(2) - x1(1)) / (length(x1) / length(x2));
			dx = x0 - range / 2;
			chi2 = diff_sq_sum(x1, y1, x2 + dx, y2);
			min_chi2 = chi2;
			min_dx = dx;
			while (dx < x0 + range / 2)
				chi2 = diff_sq_sum(x1, y1, x2 + dx, y2);
				if (chi2 < min_chi2)
					min_chi2 = chi2;
					min_dx = dx;
				end
				dx = dx + step;
			end
			dx = min_dx;
		otherwise
			dx = NaN;
	end
		
end