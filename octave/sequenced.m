function s = sequenced(min, max, n)
	if (min >= max)
		error('Maximum limit must be strictly greater than minimum limit.');
	elseif (n < 2)
		error('Inconsistent density request: density must be greater than 1.');
	elseif (n == 2)
		s = [min, max];
	else
		s = min:((max - min) / (n - 1)):max;
	end
end