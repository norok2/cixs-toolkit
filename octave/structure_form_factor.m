function [F, chi] = structure_form_factor(crystal, E, hkl, fourier_sign, im_form_factor_sign, susceptivity_sign)
	% structure form factor depending on atom position in unit cell
	
	if (nargin < 4)
		fourier_sign = -1;
	end
	if (nargin < 5)
		im_form_factor_sign = +1;
	end
	if (nargin < 6)
		susceptivity_sign = -1;
	end

	% classical electron radius (r_e = e^2 / m_e / c^2) in m
	[r_e, UM_r_e, D_r_e] = physical_constant('CLASSICAL_ELECTRON_RADIUS');
	
	% momentum transfer Q
	lambda = conv_E('eV', 'm', E);
	Q = 4 * pi ./ lambda .* sin(bragg_angle(crystal, E, hkl));
	% Debye-Waller
	DW = debye_waller(crystal, E, hkl);
	F = zeros(size(E));
	% structure form factor of Si at hkl direction depends on Q
	for i=1:length(crystal.unit_cell(:,1))
		% atomic form factor (including anomalous dispersion corrections)
		FF = atomic_form_factor(crystal.unit_cell(i,5), E, Q, im_form_factor_sign);
		% occupancy of selected atom
		n = crystal.unit_cell(i,4);
		% vector position of atom
		r = (crystal.unit_cell(i,1:3) - crystal.cell_offset);
		% compute structure form factor
		F = F + (FF .* exp(-DW + 2 .* pi .* 1i .* fourier_sign .* n .* (hkl * r')));
	end
	% calculate correspondent susceptivity
	chi = F .* susceptivity_sign .* r_e .* lambda.^2.0 ./ pi ./ unit_cell_volume(crystal.cell_size, crystal.cell_angle);
end