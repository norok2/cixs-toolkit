function data = compute_dsf_svd(data, expt)
	% compute dynamic structure factor using SVD method

	% Elementary charge e in C
	[q_e, UM_q_e, D_q_e] = physical_constant('ELEMENTARY_CHARGE');
	% classical electron radius (r_e = e^2 / m_e / c^2) in m
	[r_e, UM_r_e, D_r_e] = physical_constant('CLASSICAL_ELECTRON_RADIUS');	
	
	beam = expt.beam;
	sample = expt.sample;
	monochr = expt.monochr;
	n_terms = expt.sample.n_terms;

	% other prefactors
	%CI = r_e^2.0 / 4.0 / pi^2.0 / q_e^2.0;
	CI = 1.0;
	ea = 1; eb = 1; ec = 1; ed = ec;

	% calculate weight matrix
	data.val.svd.aa = mean(data.val.sw.aa);
	data.val.svd.bb = mean(data.val.sw.bb);
	data.val.svd.cc = mean(data.val.sw.cc);
	data.val.svd.dd = mean(data.val.sw.dd);
	data.val.svd.rr = mean(data.val.sw.rr);
	data.val.svd.D_aa = std(data.val.sw.aa);
	data.val.svd.D_bb = std(data.val.sw.bb);
	data.val.svd.D_cc = std(data.val.sw.cc);
	data.val.svd.D_dd = std(data.val.sw.dd);
	data.val.svd.D_rr = std(data.val.sw.rr);
	data.val.svd.W = CI .* [ea .* data.val.svd.aa', eb .* data.val.svd.bb', ec .* data.val.svd.cc', ed .* data.val.svd.dd'];
	switch n_terms
		case 2
			selected_components = [1,2];
			data.val.svd.Wr = [data.val.svd.W(:,1) + data.val.svd.W(:,2), data.val.svd.W(:,3)];
		case 3
			selected_components = [1,3,2];
			data.val.svd.Wr = [data.val.svd.W(:,1), data.val.svd.W(:,3), data.val.svd.W(:,2)];
		case 4
			selected_components = [1,3,2,4];
			data.val.svd.Wr = [data.val.svd.W(:,1), data.val.svd.W(:,3), data.val.svd.W(:,2), data.val.svd.W(:,4)];;
		otherwise
			error('Too many terms in equations for non-diagonal dynamic structure factor determination');
	end

	% compute singular value decomposition
	[data.val.svd.U, data.val.svd.S, data.val.svd.V] = svd(data.val.ss, 'econ'); % D = U * S * V'
	data.val.svd.US = data.val.svd.U * data.val.svd.S;
	% truncate to principal values
	data.val.svd.USr = data.val.svd.US(:,selected_components);
	data.val.svd.Vr = data.val.svd.V(:,selected_components);
	% determine reduced T matrix: D = US * V' = US * TI * T * V'
	data.val.svd.Tr = data.val.svd.Wr' / data.val.svd.Vr';
	% calculate dsf spectra
	data.val.svd.dsf = data.val.svd.USr / data.val.svd.Tr;

	% error calculation
	% upper svd's dsf calculation
	u_ss = data.val.ss + data.val.D_ss;
	[uU, uS, uV] = svd(u_ss, 'econ');
	uUS = uU * uS;
	uUSr = uUS(:,selected_components);
	uVr = uV(:,selected_components);
	uTr = data.val.svd.Wr' / uVr';
	u_dsf = uUSr / uTr;
	% lower svd's dsf calculation
	l_ss = data.val.ss - data.val.D_ss;
	[lU, lS, lV] = svd(l_ss, 'econ');
	lUS = lU * lS;
	lUSr = lUS(:,selected_components);
	lVr = lV(:,selected_components);
	lTr = data.val.svd.Wr' / lVr';
	l_dsf = lUSr / lTr;
	% dsf error calculation
	data.val.svd.D_dsf = abs(u_dsf - l_dsf) / 2.0;
%  	[dU, dS, dV] = svd(data.val.D_ss, 'econ'); % D = U * S * V'
%  	dUS = dU * dS;
%  	% truncate to principal values
%  	dUSr = dUS(:,selected_components);
%  	dVr = dV(:,selected_components);
%  	% determine reduced T matrix: D = US * V' = US * TI * T * V'
%  	dTr = data.val.svd.Wr' / dVr'; % 
%  	% calculate dsf spectra error
%  	data.val.svd.D_dsf = dUSr / dTr;

	data.val.svd.result = dsf_normalization(data.val.dE, data.val.svd.dsf(:,1)', data.val.svd.D_dsf(:,1)', data.val.svd.dsf(:,2)', data.val.svd.D_dsf(:,2)', expt.sample);
end
