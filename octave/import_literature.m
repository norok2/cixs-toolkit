function ldata = import_literature(cdata, workdir)
	% import comparison with literature

	ldata.available = [];
	ldata.name = '';
	ldata.file = '';
	ldata.data = [];
	ldata.norm = 1;
	ldata.magnify = 1;
	ldata.shift = 0;
	ldata.rho_ratio = [];
	ldata.D_rho_ratio = [];
	ldata = repmat(ldata, 1, length(cdata));

	% Planck constant h in J s
	[h_P, UM_h_P, D_h_P] = physical_constant('PLANCK_CONSTANT');
	h_Pb = h_P ./ 2.0 ./ pi; UM_h_Pb = UM_h_P; D_h_Pb = D_h_P;
	% mass of the electron in in kg
	[m_e, UM_m_e, D_m_e] = physical_constant('ELECTRON_MASS');

	for i=1:length(ldata)
		ldata(i).name = cdata(i).name;
		ldata(i).file = [workdir 'literature_' 'dsf_' ldata(i).name '.csv'];
		if (~isempty(dir(ldata(i).file)))
			ldata(i).available = true;
			ldata(i).data = dlmread(ldata(i).file, ',', 2, 0);
			% height adjustment
			ldata(i).norm = max(cdata(i).ddsf) ./ max(ldata(i).data(:,2));
			% width adjustment
			[B, x0] = fwhm(ldata(i).data(:,1), ldata(i).data(:,2));
			[BB, xx0] = fwhm(cdata(i).dE, cdata(i).ddsf);
			ldata(i).magnify = xx0 ./ x0;
			% shift adjustment
			ldata(i).shift = x_align(cdata(i).dE, cdata(i).ddsf, ldata(i).data(:,1), ldata(i).data(:,2) .* ldata(i).norm);
			
			% calculate rho_ratio
			q0 = cdata(i).q0;
			qh = cdata(i).qh;
			dE = ldata(i).data(:,1) .* ldata(i).magnify + ldata(i).shift;
			ddsf = ldata(i).data(:,2) .* ldata(i).norm;
			D_ddsf = abs(ddsf) .* (5.0 / 100.0);
			ndsf = ldata(i).data(:,3) .* ldata(i).norm;
			D_ndsf = abs(ndsf) .* (5.0 / 100.0);
			I_d = trapz(dE, dE .* ddsf);
			D_I_d = trapz(dE, dE .* D_ddsf);
			I_n = trapz(dE, dE .* ndsf);
			D_I_n = trapz(dE, dE .* D_ndsf);
			ldata(i).rho_ratio = (I_n / I_d) * ((q0 * q0') / (q0 * qh')); % rho_g / rho_0
			ldata(i).D_rho_ratio = (abs(D_I_d ./ I_d) + abs(D_I_n ./ I_n)) .* ldata(i).rho_ratio;
		else
			ldata(i).available = false;
		end
	end
end
