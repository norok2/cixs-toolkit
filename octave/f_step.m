function y = f_step(x)
	% calculate step function (with 0-value equal to 1/2)
	
	y = (1 + sign(x)) ./ 2;
end
