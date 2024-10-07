function DW = debye_waller(crystal, E, hkl, energy_unit)

	% Planck constant in J s
	[h_P, UM_h_P, D_h_P] = physical_constant('PLANCK_CONSTANT');
	% Speed of Light in Vacuum in m / s
	[c_p, UM_c_p, D_c_p] = physical_constant('SPEED_OF_LIGHT_IN_VACUUM');
	% k_B in J / K
	[k_B, UM_k_B, D_k_B] = physical_constant('BOLTZMANN_CONSTANT');
	%  Atomic mass unit in J
	[u_m, UM_u_m, D_u_m] = physical_constant('ATOMIC_MASS_UNIT_JOULE');

	% if not specified, assume energy given in 'eV'
	if (nargin < 4)
		energy_unit = 'eV';
	end
	
	if (crystal.symbol == 'Si')
		% experimental parameter to calculare specific Debye-Waller 
		Theta = 645;
		elem = chemical_properties(crystal.unit_cell(1,5));
		theta_B = bragg_angle(crystal, E, hkl, energy_unit);
		lambda = conv_E(energy_unit, 'm', E);
		B_T = (6 * h_P^2 / (elem.A * u_m) / (k_B * Theta) * c_p^2) .* ((quad(@(x) (x ./ (exp(x) - 1)), 0, Theta ./ crystal.T) ./ (Theta ./ crystal.T).^2) + 1/4);
		DW = B_T * (sin(theta_B)/lambda)^2;
	else
		DW = 0;
	end
end
