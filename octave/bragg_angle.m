function th_B = bragg_angle(crystal, E, hkl, energy_unit)
	% calculate Bragg angle of given crystal at hkl reflection for given energy

	% if not specified, assume energy given in 'eV'
	if (nargin < 4)
		energy_unit = 'eV';
	end

	% if not specified, try getting hkl from crystal
	if ((nargin < 3) && (isfield(crystal, 'hkl') > 0))
		hkl = crystal.hkl
	end

	if (exist('hkl') == false)
		error('[h,k,l] direction not specified and not set in ''crystal''.');
	end
	
	lambda = conv_E(energy_unit, 'm', E);
	d = bragg_distance(crystal, hkl);
	if (d > 0)
		th_B = asin(lambda * bragg_order(crystal, hkl) / 2 / d);
	else
		th_B = 0;
	end
end