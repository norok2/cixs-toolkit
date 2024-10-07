function diffsum = diff_sq_sum(x1, y1, x2, y2)
	% calculate the normalized sum of the square of the differences between the two numeric functions (x1,y1) and (x2,y2)

	yy = interp1(x1, y1, x2, 'linear');
	valid = ~isnan(yy);
	weight = sum(abs(valid));
	if (weight == 0)
		diffsum = Inf;
	else
		diffsum = sum(sqrt((y2(valid) - yy(valid)).^2.0)) / weight;
	end
end