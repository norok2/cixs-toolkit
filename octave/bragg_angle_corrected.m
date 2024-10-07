function th_B = bragg_angle_corrected(crystal, E, hkl, energy_unit)
	% calculate Bragg angle of given crystal at hkl reflection for given energy

	% if not specified, assume energy given in 'eV'
	if (nargin < 4)
		energy_unit = 'eV';
	end
	
	lambda = conv_E(energy_unit, 'm', E);
	d = bragg_distance(crystal, hkl);
	if (d > 0)
		[F_0, chi_0] = structure_form_factor(crystal, E, [0 0 0]);
		th_B = asin(lambda * bragg_order(crystal, hkl) / 2 / d);
		gamma_0 = cos(pi/2 - th_B + crystal.asymmetry);
		gamma_h = cos(pi/2 + th_B + crystal.asymmetry);
		gamma = gamma_h ./ gamma_0;
		th_B = th_B - real(chi_0) .* (1 - gamma) ./ (2 .* sin(2 .* th_B));
	else
		th_B = 0;
	end
end