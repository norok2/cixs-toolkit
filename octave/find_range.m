function id = find_range(v, min, max)
	id = find(v < max);
	v_max = v(id);
	id = find(v_max > min);
end