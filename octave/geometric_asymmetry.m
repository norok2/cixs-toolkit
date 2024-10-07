function alpha_o = geometric_asymmetry(alpha_i, phi)
	alpha_o = atan(tan(alpha_i) ./ cos(phi));
end