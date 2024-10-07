function [res, D_res] = conv_E(iUM, oUM, val, D_val)
	% convert energy values from/to [eV, J, Hz, m^-1, K, m, cm^-1, nm] (needs 'physical_constant')

	% Planck constant in J s
	[h, UM_h, D_h] = physical_constant('PLANCK_CONSTANT');
	% Speed of Light in Vacuum (exact: error is due uncertainty in meter) in m / s
	[c, UM_c, D_c] = physical_constant('SPEED_OF_LIGHT_IN_VACUUM');
	% k_B in J / K
	[k_B, UM_k_B, D_k_B] = physical_constant('BOLTZMANN_CONSTANT');
	% Elementary charge in C
	[q_e, UM_q_e, D_q_e] = physical_constant('ELEMENTARY_CHARGE');
	% Celsius to Kelvin conversion constant (exact): T(°C) = T(K) - T_0
	T_0 = 273.15;
	
	% frequently used combination
	hcq = h * c / q_e;
	
	switch (iUM)
		case {'eV'}
			res = val;
		case {'J'}
			res = val ./ q_e;
		case {'Hz'}
			res = (h / q_e) .* val;
		case {'rad*Hz'}
			res = (h / q_e / 2 / pi) .* val;
		case {'m^-1'}
			res = hcq .* val;
		case {'C', '°C'}
			res = (val + T_0) .* (k_B / q_e);
		case {'K'}
			res = val .* (k_B / q_e);
		case {'m'}
			res = hcq ./ val;
		case {'cm^-1'}
			res = (hcq * 1e2) .* val;
		case {'nm'}
			res = (hcq * 1e9) ./ val;
		otherwise
			error('invalid input unit of measure');
	end
    
	switch oUM
		case {'eV'}
			%res = res;
		case {'J'}
			res = res .* q_e;
		case {'Hz'}
			res = res .* (q_e / h);
		case {'rad*Hz'}
			res = res .* (2 * pi * q_e / h);
		case {'m^-1'}
			res = res ./ hcq;
		case {'C', '°C'}
			res = (res .* (q_e / k_B)) - T_0;
		case {'K'}
			res = res .* (q_e / k_B);
		case {'m'}
			res = hcq ./ res;
		case {'cm^-1'}
			res = res ./ (1e2 * hcq);
		case {'nm'}
			res = (1e9 * hcq) ./ res;
		otherwise
			error('invalid output unit of measure');
	end

	if (nargin >= 4)
		% assuming percentage error much bigger than error on physical constants
		D_res = res .* (D_val ./ val);
	else
		% assuming error only from physical constants
		D_res = res .* max([(D_h / h), (D_c / c), (D_k_B / k_B), (D_q_e / q_e)]);
	end
end
