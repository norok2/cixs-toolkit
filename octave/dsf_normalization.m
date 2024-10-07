function result = dsf_normalization(dE, ddsf, D_ddsf, ndsf, D_ndsf, sample)

	% Planck constant h in J s
	[h_P, UM_h_P, D_h_P] = physical_constant('PLANCK_CONSTANT');
	h_Pb = h_P ./ 2.0 ./ pi; UM_h_Pb = UM_h_P; D_h_Pb = D_h_P;
	[h_Pb_eV, D_h_Pb_eV] = conv_E('J', 'eV', h_Pb, D_h_Pb); UM_h_Pb_eV = 'eV';
	% mass of the electron in in kg
	[m_e, UM_m_e, D_m_e] = physical_constant('ELECTRON_MASS');
	% Elementary charge e in C
	[q_e, UM_q_e, D_q_e] = physical_constant('ELEMENTARY_CHARGE');
	% Speed of Light (Photon) in Vacuum c in m / s
	[c_p, UM_c_p, D_c_p] = physical_constant('SPEED_OF_LIGHT_IN_VACUUM');
	% dielectric constant in vacuum epsilon_0 in F m^-1
	[eh_0, UM_eh_0, D_eh_0] = physical_constant('ELECTRIC_CONSTANT');
	% classical electron radius (r_e = e^2 / m_e / c^2) in m
	[r_e, UM_r_e, D_r_e] = physical_constant('CLASSICAL_ELECTRON_RADIUS');

	% exchanged vector momentum in 1/m
	q_0 = (2.0 * pi / bragg_distance(sample, sample.q)) * normalize(sample.q);
	q_h = (2.0 * pi / bragg_distance(sample, (sample.q + sample.hkl))) * normalize(sample.q + sample.hkl);

	% normalization using: \int_0^\inf {E S(\vect{q}, E) dE} = \frac{\hbar^2 q^2}{2 m_e}
	integr_ddsf = trapz(dE, dE .* ddsf);
	D_integr_ddsf = trapz(dE, dE .* D_ddsf);
	norm_ddsf = (h_Pb.^2.0 .* (q_0 * q_0') / 2.0 / m_e);
	norm_ddsf_eV = conv_E('J', 'eV', norm_ddsf);
	result.alpha_ddsf = integr_ddsf ./ norm_ddsf_eV;
	% normalization check of non-diagonal part of dsf
	integr_ndsf = trapz(dE, dE .* ndsf);
	D_integr_ndsf = trapz(dE, dE .* D_ndsf);
	norm_ndsf = (h_Pb.^2.0 .* (q_0 * q_h') / 2.0 / m_e);
	norm_ndsf_eV = conv_E('J', 'eV', norm_ndsf);
	result.alpha_ndsf = integr_ndsf ./ norm_ndsf_eV;
	result.rho_ratio = result.alpha_ndsf / result.alpha_ddsf; % rho_g / rho_0
	result.D_rho_ratio = (abs(D_integr_ddsf ./ integr_ddsf) + abs(D_integr_ndsf ./ integr_ndsf)) .* result.rho_ratio;

	% normalized dynamic structure factor in eV^-1
	result.ddsf = ddsf ./ result.alpha_ddsf;
	result.D_ddsf = D_ddsf ./ result.alpha_ddsf;
	result.ndsf = ndsf ./ result.alpha_ddsf;
	result.D_ndsf = D_ndsf ./ result.alpha_ddsf;
	
	% normalized imaginary inverse dielectric matrix (dimensionless)
	n_e = 4 * 8 * unit_cell_volume(sample.cell_size, sample.cell_angle); % in m % valence electron density in unit cell
	eps_factor = 4.0 * pi^2 * q_e^2.0 * n_e / eh_0; % J * m^2
	eps_factor_eV = conv_E('J', 'eV', eps_factor);
	result.deps = result.ddsf .* eps_factor_eV ./ (q_0 * q_0');
	result.D_deps = result.D_ddsf .* eps_factor_eV ./ (q_0 * q_0');
	result.neps = result.ndsf .* eps_factor_eV ./ (abs(q_0) * abs(q_h)');
	result.D_neps = result.D_ndsf .* eps_factor_eV ./ (abs(q_0) * abs(q_h)');
end
