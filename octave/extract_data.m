function data = extract_data(src, param, E_in, beam, monochr, sample)
	% extract data

	% assume:
	% - energy scans all have different energy
	% - angular oversampling does not introduce errors
	
	% spec data column identifiers
	spec.id.th = 1; % Bragg angle (schi2)
	spec.id.E = 11; % energy shift
	spec.id.mo = 16; % monitor
	spec.id.ss = 17; % signal intensity
	spec.id.rr = 6; % Bragg intensity

	spec.offset.align = 1;
	spec.offset.main = 2;
	spec.offset.tail = 3;
	spec.offset.len = length(fieldnames(spec.offset));

	raw = specread(src.file);

	% new data column identifiers
	id.th = 1; % Bragg angle (schi2)
	id.E = 2; % E shift
	id.mo = 3; % monitor
	id.ss = 4; % signal intensity
	id.rr = 5; % Bragg intensity

	%% extract signal
	sequence.main = src.scan(spec.offset.main:spec.offset.len:end);
	sequence.tail = src.scan(spec.offset.tail:spec.offset.len:end);
	raw_data.num.th_main = length(raw(sequence.main(1)).y(:,1));
	raw_data.num.th_tail = length(raw(sequence.tail(1)).y(:,1));
	raw_data.num.th = raw_data.num.th_main + raw_data.num.th_tail;
	raw_data.num.E = length(sequence.main);
	raw_data.num.pt = raw_data.num.th * raw_data.num.E;

	% construct data matrix
	data.shift = zeros(1, raw_data.num.E);
	raw_data.dth = zeros(raw_data.num.E, raw_data.num.th); raw_data.dE = zeros(raw_data.num.E, raw_data.num.th); raw_data.ss = zeros(raw_data.num.E, raw_data.num.th); raw_data.D_ss = zeros(raw_data.num.E, raw_data.num.th); raw_data.rr = zeros(raw_data.num.E, raw_data.num.th);
	raw_data.min.left_tail.dth = 1; raw_data.max.left_tail.dth = 0;
	raw_data.min.main.dth = 0; raw_data.max.main.dth = 0;
	raw_data.min.right_tail.dth = 0; raw_data.max.right_tail.dth = -1;
	for i=1:raw_data.num.E
		% extract data
		range_l_i = 1:raw_data.num.th_tail/2; % left input
		range_m_i = 1:raw_data.num.th_main; % middle input
		range_r_i = raw_data.num.th_tail/2+1:raw_data.num.th_tail;	% right input	
		range_l_o = range_l_i; % left output
		range_m_o = length(range_l_i) + range_m_i; % middle output
		range_r_o = length(range_l_i) + length(range_m_i) - length(range_r_i) + range_r_i; % right output
		% extract left tail part of data
		raw_data.val(range_l_o,i,id.E) = raw(sequence.tail(i)).y(range_l_i,spec.id.E);
		raw_data.val(range_l_o,i,id.th) = raw(sequence.tail(i)).y(range_l_i,spec.id.th);
		raw_data.val(range_l_o,i,id.mo) = raw(sequence.tail(i)).y(range_l_i,spec.id.mo);
		raw_data.val(range_l_o,i,id.ss) = raw(sequence.tail(i)).y(range_l_i,spec.id.ss);
		raw_data.val(range_l_o,i,id.rr) = raw(sequence.tail(i)).y(range_l_i,spec.id.rr);
		% extract central part of data
		raw_data.val(range_m_o,i,id.E)=raw(sequence.main(i)).y(range_m_i,spec.id.E);
		raw_data.val(range_m_o,i,id.th)=raw(sequence.main(i)).y(range_m_i,spec.id.th);
		raw_data.val(range_m_o,i,id.mo)=raw(sequence.main(i)).y(range_m_i,spec.id.mo);
		raw_data.val(range_m_o,i,id.ss)=raw(sequence.main(i)).y(range_m_i,spec.id.ss);
		raw_data.val(range_m_o,i,id.rr)=raw(sequence.main(i)).y(range_m_i,spec.id.rr);
		% extract right tail part of data
		raw_data.val(range_r_o,i,id.E)=raw(sequence.tail(i)).y(range_r_i,spec.id.E);
		raw_data.val(range_r_o,i,id.th)=raw(sequence.tail(i)).y(range_r_i,spec.id.th);
		raw_data.val(range_r_o,i,id.mo)=raw(sequence.tail(i)).y(range_r_i,spec.id.mo);
		raw_data.val(range_r_o,i,id.ss)=raw(sequence.tail(i)).y(range_r_i,spec.id.ss);
		raw_data.val(range_r_o,i,id.rr)=raw(sequence.tail(i)).y(range_r_i,spec.id.rr);

		% energy and energy loss
		dE_exp = fliplr(raw_data.val(:,i,id.E)');
		E_exp = dE_exp + E_in;
		% angle and Bragg angle
		th_B = bragg_angle(sample, E_exp, sample.hkl);
		th_exp = geometric_angle(fliplr(deg2rad(abs(raw_data.val(:,i,id.th)'))), mean(th_B), sample.phi, beam);
		dth_exp = angular_departure(th_exp, th_B);
		% signal, monitor, reflectivity and errors
		ss_exp = fliplr(raw_data.val(:,i,id.ss)'); D_ss_exp = sqrt(ss_exp);
		mo_exp = fliplr(raw_data.val(:,i,id.mo)'); D_mo_exp = sqrt(mo_exp);
		rr_exp = fliplr(raw_data.val(:,i,id.rr)') ./ mo_exp;
		% calculate error on signal
		diff_ss_var = [1./mo_exp; -ss_exp./mo_exp.^2.0]; 
		D_ss_var = [D_ss_exp; D_mo_exp];
		ss_exp = ss_exp ./ mo_exp;
		D_ss_exp = zeros(size(ss_exp));
		for j=1:length(ss_exp)
			for k=1:length(D_ss_var(:,1))
				D_ss_exp(j) = D_ss_exp(j) + (D_ss_var(k,j) * diff_ss_var(k,j))^2.0;
			end
		end
		D_ss_exp = sqrt(D_ss_exp);

		% calculate theoretical quantities for alignment
		if (param.precise_convolution)
			% calculate E sampling in eV
			E_min = E_exp - param.E.bound;
			E_max = E_exp + param.E.bound;
			E = sequenced(E_min, E_max, param.E.density);
			% calculate th sampling in rad
			th_min = min(bragg_angle(monochr, E_max, monochr.hkl), bragg_angle(sample, E_max, sample.hkl)) + param.dth.min;
			th_max = max(bragg_angle(monochr, E_min, monochr.hkl), bragg_angle(sample, E_min, sample.hkl)) + param.dth.max;
			th = sequenced(th_min, th_max, param.th.density);
			% WARNING: normalization issues?
			beam.bragg_angle = bragg_angle(monochr, mean(E_exp), monochr.hkl);
			ibt = zeros(length(th), length(E)); mmt = zeros(length(th), length(E)); rrt = zeros(length(th), length(E));
			for j=1:length(E)
				% ibt beam before monochromator
				ibt(:,j) = beam_profile(th, E(j), 'gauss', [beam.bragg_angle, beam.divergence]);
				for nu=1:2
					% monochromator transfer function
					mmt(:,j) = mmt(:,j)' + beam.polarization(nu) .* reflectivity(th, E(j), monochr, monochr.hkl, nu);
					% reflectivity of the sample
					rrt(:,j) = rrt(:,j)' + beam.polarization(nu) .* reflectivity(th, E(j), sample, sample.hkl, nu);
				end
			end
			% 2D convolution for (E,th) dispersion, integrated over Energy
			mmt = normalize(ibt .* mmt);
			rr_th = iconv2(rrt', mmt', 'same');
		else
			% calculate th sampling in rad
			th_min = bragg_angle(sample, E_exp, sample.hkl) + param.dth.min;
			th_max = bragg_angle(sample, E_exp, sample.hkl) + param.dth.max;
			th = sequenced(th_min, th_max, param.th.density);
			thmm_min = bragg_angle(monochr, E_exp, monochr.hkl) + param.dth.min;
			thmm_max = bragg_angle(monochr, E_exp, monochr.hkl) + param.dth.max;
			thmm = sequenced(thmm_min, thmm_max, param.th.density);
			% calculate theoretical rocking curve
			mmt = zeros(size(th)); rrt = zeros(size(th));
			for nu=1:2
				% monochromator transfer function
				mmt = mmt + beam.polarization(nu) .* reflectivity(thmm, E_in, monochr, monochr.hkl, nu);
				% reflectivity of the sample
				rrt = rrt + beam.polarization(nu) .* reflectivity(th, mean(E_exp), sample, sample.hkl, nu);
			end
			mmt = normalize(mmt); % should solve normalization issues
			rr_th = iconv2(rrt, mmt, 'same');
		end

		% width adjustment and alignment of experimental data using reflectivity
		dth = angular_departure(th, mean(th_B));
		rrn_exp = rr_exp ./ max(rr_exp) .* max(rr_th);
		align_exp_range = (rrn_exp > max(rrn_exp) / 2.0);
		% alignment
		x_shift = x_align(dth, rr_th, dth_exp(align_exp_range), rrn_exp(align_exp_range));
		dth_exp = dth_exp + x_shift;
		data.shift(i) = x_shift;
		% adjustment (?)
		fwhm_exp = fwhm(dth_exp, rr_exp);
		fwhm_th = fwhm(dth, rr_th);
		dth_exp = dth_exp .* fwhm_th ./ fwhm_exp;
		% re-alignment (?)
		x_shift = x_align(dth, rr_th, dth_exp(align_exp_range), rrn_exp(align_exp_range));
		dth_exp = dth_exp + x_shift;

		% calculate matrix of data
		data_raw.dth(i,:) = dth_exp;
		data_raw.dE(i,:) = dE_exp;
		data_raw.ss(i,:) = ss_exp;
		data_raw.D_ss(i,:) = D_ss_exp;
		data_raw.rr(i,:) = rrn_exp;

		% find limits
		if (min(dth_exp(range_l_o)) < raw_data.min.left_tail.dth)
			raw_data.min.left_tail.dth = min(dth_exp(range_l_o));
		elseif (max(dth_exp(range_l_o)) > raw_data.max.left_tail.dth)
			raw_data.max.left_tail.dth = max(dth_exp(range_l_o));
		end
		if (min(dth_exp(range_m_o)) < raw_data.min.main.dth)
			raw_data.min.main.dth = min(dth_exp(range_m_o));
		elseif (max(dth_exp(range_m_o)) > raw_data.max.main.dth)
			raw_data.max.main.dth = max(dth_exp(range_m_o));
		end
		if (min(dth_exp(range_r_o)) < raw_data.min.right_tail.dth)
			raw_data.min.right_tail.dth = min(dth_exp(range_r_o));
		elseif (max(dth_exp(range_r_o)) > raw_data.max.right_tail.dth)
			raw_data.max.right_tail.dth = max(dth_exp(range_r_o));
		end 
	end

	% generate interpolated and mesh grid for signal and reflectivity
	% NOTE: assume tail has only 1 point
	% NOTE: left and right got inverted while importing data (should be: left < right)
	left_tail_dth = mean(data_raw.dth(:,1));
	right_tail_dth = mean(data_raw.dth(:,end));
	main_dth = sequenced(raw_data.min.main.dth, raw_data.max.main.dth, round(raw_data.num.th_main * param.dth.sampling_ratio));
	data.dth_raw = [left_tail_dth, main_dth, right_tail_dth];
	data.dE_raw = mean(data_raw.dE');
	data.rr_raw = zeros(length(data.dE_raw), length(data.dth_raw)); data.ss_raw = zeros(length(data.dE_raw), length(data.dth_raw)); data.D_ss_raw = zeros(length(data.dE_raw), length(data.dth_raw));
	for i=1:raw_data.num.E
		% reflectivity
		data.rr_raw(i,:) = interp1(data_raw.dth(i,:), data_raw.rr(i,:), data.dth_raw);
		if (isnan(data.rr_raw(i,1)))
			data.rr_raw(i,1) = data_raw.rr(i,1);
		end
		if (isnan(data.rr_raw(i,end)))
			data.rr_raw(i,end) = data_raw.rr(i,end);
		end
		% signal
		data.ss_raw(i,:) = interp1(data_raw.dth(i,:), data_raw.ss(i,:), data.dth_raw);
		if (isnan(data.ss_raw(i,1)))
			data.ss_raw(i,1) = data_raw.ss(i,1);
		end
		if (isnan(data.ss_raw(i,end)))
			data.ss_raw(i,end) = data_raw.ss(i,end);
		end
		% error on signal
		data.D_ss_raw(i,:) = interp1(data_raw.dth(i,:), data_raw.D_ss(i,:), data.dth_raw);
		if (isnan(data.D_ss_raw(i,1)))
			data.D_ss_raw(i,1) = data_raw.D_ss(i,1);
		end
		if (isnan(data.D_ss_raw(i,end)))
			data.D_ss_raw(i,end) = data_raw.D_ss(i,end);
		end
	end
	[data.mesh_dth_raw, data.mesh_dE_raw] = meshgrid(data.dth_raw, data.dE_raw);
end
