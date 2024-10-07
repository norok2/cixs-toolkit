function data = import_literature_expt(data, expt, workdir)
	% import comparison with literature

	% Planck constant h in J s
	[h_P, UM_h_P, D_h_P] = physical_constant('PLANCK_CONSTANT');
	h_Pb = h_P ./ 2.0 ./ pi; UM_h_Pb = UM_h_P; D_h_Pb = D_h_P;
	% mass of the electron in in kg
	[m_e, UM_m_e, D_m_e] = physical_constant('ELECTRON_MASS');

	sample = expt.sample;

	% exchanged vector momentum in 1/m
	q_0 = (2.0 * pi / bragg_distance(sample, sample.q)) * normalize(sample.q);
	q_h = (2.0 * pi / bragg_distance(sample, (sample.q + sample.hkl))) * normalize(sample.q + sample.hkl);

	data.lit.file = [workdir 'literature_dsf_' expt.name(1:6) '.csv'];
	if (~isempty(dir(data.lit.file)))
		data.lit.data = dlmread('import/literature_dsf_hkl456_1.csv', ',', 2, 0);
		data.lit.norm = max(data.val.result.ddsf) ./ max(data.lit.data(:,2));
		data.lit.shift = x_align(data.val.dE, data.val.result.ddsf, data.lit.data(:,1), data.lit.data(:,2) .* data.lit.norm);
		
		% calculate rho_ratio
		dE = data.lit.data(:,1) + data.lit.shift;
		ddsf = data.lit.data(:,2) .* data.lit.norm;
		D_ddsf = abs(ddsf) .* (5.0 / 100.0);
		ndsf = data.lit.data(:,3) .* data.lit.norm;
		D_ndsf = abs(ndsf) .* (5.0 / 100.0);
		integr_ddsf = trapz(dE, dE .* ddsf);
		D_integr_ddsf = trapz(dE, dE .* D_ddsf);
		norm_ddsf = (h_Pb.^2.0 .* (q_0 * q_0') / 2.0 / m_e);
		norm_ddsf_eV = conv_E('J', 'eV', norm_ddsf);
		alpha_ddsf = integr_ddsf ./ norm_ddsf_eV;
		integr_ndsf = trapz(dE, dE .* ndsf);
		D_integr_ndsf = trapz(dE, dE .* D_ndsf);
		norm_ndsf = (h_Pb.^2.0 .* (q_0 * q_h') / 2.0 / m_e);
		norm_ndsf_eV = conv_E('J', 'eV', norm_ndsf);
		alpha_ndsf = integr_ndsf ./ norm_ndsf_eV;
		data.lit.rho_ratio = alpha_ndsf / alpha_ddsf; % rho_g / rho_0
		data.lit.D_rho_ratio = (abs(D_integr_ddsf ./ integr_ddsf) + abs(D_integr_ndsf ./ integr_ndsf)) .* data.lit.rho_ratio;
	else
		data.lit = [];
	end
end
