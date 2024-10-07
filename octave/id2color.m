function c = id2color(i, method)
	% associate colors to index (useful for plot color cycling)

	if (nargin < 2)
		method = 'full';
	end

	array_colors = color_list(method);
	c = array_colors(mod(round(i) - 1, length(array_colors)) + 1);
end
