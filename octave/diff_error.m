function [dv, Ddv] = diff_error(v, Dv)
	% calculate diff of v and its error
	
	% calculate diff of v
	dv = diff(v);
	
	% calculate its error
	alternate_1 = (-1).^(1:length(Dv));
	Ddv_1 = diff(v + alternate_1 .* Dv); 
	alternate_2 = (-1) .* (-1).^(1:length(Dv));
	Ddv_2 = diff(v + alternate_2 .* Dv);
	Ddv = max(Ddv_1, Ddv_2) - min(Ddv_1, Ddv_2);
end