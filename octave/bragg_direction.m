function direction = bragg_direction(crystal, E, hkl, beam, energy_unit)
	% calculate Bragg direction respect to beam

	% if not specified, assume energy given in 'eV'
	if (nargin < 5)
		energy_unit = 'eV';
	end
	
	theta_B = bragg_angle(crystal, E, hkl, energy_unit);
	direction = cos(theta_B) .* beam.direction + sin(theta_B) .* beam.normal;
end