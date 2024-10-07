function data = compute_dsf_linsolve(data, expt, tuples)
	% compute dynamic structure factor using linear system solution method

	% Elementary charge e in C
	[q_e, UM_q_e, D_q_e] = physical_constant('ELEMENTARY_CHARGE');
	% classical electron radius (r_e = e^2 / m_e / c^2) in m
	[r_e, UM_r_e, D_r_e] = physical_constant('CLASSICAL_ELECTRON_RADIUS');	
	
	beam = expt.beam;
	sample = expt.sample;
	monochr = expt.monochr;
	n_terms = expt.sample.n_terms;

	% other prefactors
	CI = r_e^2.0 / 4.0 / pi^2.0 / q_e^2.0;
	ea = 1; eb = 1; ec = 1;

	len_dE = length(data.val.dE);
	len_dth = length(data.val.dth);

	data.val.linsolve.tuples = tuples;
	% determine maximum number of good points
	data.val.linsolve.point.good = 1:len_dth;
	% determine point set
	data.val.linsolve.point.set = nchoosek(data.val.linsolve.point.good, n_terms);
	data.val.linsolve.point.set_len = length(data.val.linsolve.point.set(:,1));
	data.val.linsolve.dsf_raw = zeros(len_dE, 2 * n_terms, data.val.linsolve.point.set_len); data.val.linsolve.weight = zeros(len_dE, data.val.linsolve.point.set_len); data.val.linsolve.dsf = zeros(len_dE, n_terms); data.val.linsolve.D_dsf = zeros(len_dE, n_terms);
	for j=1:len_dE
		% experimental measures
		E_exp = data.val.dE(j) + expt.E_in;
		dth_exp = data.val.dth;
		ss_exp = data.val.ss(j,:);
		D_ss_exp = data.val.D_ss(j,:);
		rr_exp = data.val.rr(j,:);
		% theoretical quantities
		aa_th = data.val.sw.aa(j,:);
		bb_th = data.val.sw.bb(j,:);
		cc_th = data.val.sw.cc(j,:);
		dd_th = data.val.sw.dd(j,:);
		rr_th = data.val.sw.rr(j,:);
		% deviation of theoretical reflectivity from experimental
		rr_norm = norm(rr_exp);
		rr_exp = rr_exp ./ rr_norm;
		rr_th = normalize(rr_th);
		rr_dev = abs(rr_exp - rr_th) ./ rr_exp;
		%ref_rr_dev = 0.1;
		% deviation of derivatives
% 		diff_rr_th = diff(rr_th);
% 		drr_dev = [diff_rr_dev(1), (diff_rr_dev(1:end-1) + diff_rr_dev(2:end)) ./ 2.0, diff_rr_dev(end)];
% 		ref_drr_dev = 4;
		% numerical stability of theoretical values
		num_stab = aa_th.^2 .* bb_th.^2 .* cc_th.^2 .* dd_th.^2;
		% minimum angular separation for independent photons
		ref_th_dev = darwin_width(E_exp, monochr, monochr.hkl, beam.polarization);
		% dielectric tensor epsilon as a function of energy (n-point calculation)
		np=1:n_terms;
		for l=1:data.val.linsolve.point.set_len
			% forward-diffracted wave amplitude A^2
			aa = aa_th(data.val.linsolve.point.set(l,:));
			% Bragg-reflected wave amplitude B^2
			bb = bb_th(data.val.linsolve.point.set(l,:));
			% interference 2*A*B*cos(Delta)
			cc = cc_th(data.val.linsolve.point.set(l,:));
			% interference 2*A*B*sin(Delta)
			dd = dd_th(data.val.linsolve.point.set(l,:));
			% measured intensity and its interpolated background
			ss = ss_exp(data.val.linsolve.point.set(l,:));
			D_ss = D_ss_exp(data.val.linsolve.point.set(l,:));
			% diagonal term of dielectric tensor epsilon (proportional to forward-diffracted wave amplitude)
			switch n_terms
				case 2
					M = CI .* [ea .* aa' + eb .* bb', ec .* cc'];
				case 3
					M = CI .* [ea .* aa', ec .* cc', eb .* bb'];
				case 4
					M = CI .* [ea .* aa', ec .* cc', eb .* bb', ec .* dd'];
				otherwise
					error('Too many terms in equations for non-diagonal dynamic structure factor determination');
			end
			X = (M\ss')';
			D_X = (M\D_ss')';
			% collect result
			data.val.linsolve.dsf_raw(j,:,l) = [X, D_X];
			% numerical stability check
			%data.val.linsolve.chi2 = sqrt(sum(((M * X')' - ss).^2.0)) ./ (length(ss) - n_terms - 1.0);
			% terms dependency due to angular deviation
			th_dev = diff(dth_exp(data.val.linsolve.point.set(l,:)));
			% calcualte data.val.linsolve.weight(j,k,:)s
			switch tuples
				case 'all'
					data.val.linsolve.weight(j,l) = 1;
				case 'selected'
					data.val.linsolve.weight(j,l) = 1;
					% remove point sets with rr_exp below average
					m = 1;
					while ((m <= n_terms) && (data.val.linsolve.weight(j,l) == 1))
						if (rr_exp(data.val.linsolve.point.set(l,m)) < mean(rr_exp))
							data.val.linsolve.weight(j,l) = 0;
						end
						m = m + 1;
					end
					% remove point sets with th too close
					m = 1;
					while ((m <= length(th_dev)) && (data.val.linsolve.weight(j,l) == 1))
						if (th_dev(m) < ref_th_dev)
							data.val.linsolve.weight(j,l) = 0;
						end
						m = m + 1;
					end
				case 'weighted'
					data.val.linsolve.weight(j,l) = prod(num_stab(data.val.linsolve.point.set(l,:))) .* prod(1 ./ rr_dev(data.val.linsolve.point.set(l,:)).^2.0) .* prod(1 ./ drr_dev(data.val.linsolve.point.set(l,:)).^2.0) .* (exp((prod(th_dev.^2.0)/ref_th_dev)^n_terms));
			end
		end
		% calculate point average for all calculed coefficient
		v = zeros(1, data.val.linsolve.point.set_len);
		w = zeros(1, data.val.linsolve.point.set_len);
		for k=np
			v(:) = data.val.linsolve.dsf_raw(j,k,:);
			w(:) = data.val.linsolve.weight(j,:);
			[data.val.linsolve.dsf(j,k), data.val.linsolve.D_dsf(j,k)] = average(v, w.^2.0);
			
		end
	end

	data.val.linsolve.result = dsf_normalization(data.val.dE, data.val.linsolve.dsf(:,1)', data.val.linsolve.D_dsf(:,1)', data.val.linsolve.dsf(:,2)', data.val.linsolve.D_dsf(:,2)', expt.sample);
end
