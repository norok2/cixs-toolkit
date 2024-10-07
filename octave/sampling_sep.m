function v = sampling_sep(v, k)
	remove = [];
	i = 1;
	while (i < length(v))
		j = i + 1;
		while (abs(v(i) - v(j)) < k) && (j < length(v))
			remove = [remove, j];
			j = j + 1;
		end
		i = j;
	end
	v(remove) = [];
end