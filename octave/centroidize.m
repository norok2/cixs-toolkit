function vv = centroidize(v)
	l = length(v);
	n = 1:l;
	n_0 = centroid(n, v);
	shift = (round(2.0 * n_0) - 2.0 * n_0) / 4.0;
	n_min = max([1, round(2 * n_0) - l]);
	n_max = min([l, round(2.0 * (n_0 - 1)) + 1]);
	vv = interp1(n, v, (n_min:n_max) + shift, 'linear', 'extrap');
end