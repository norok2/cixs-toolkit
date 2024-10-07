function data = calculate_energy_resolution(expt)
	% plot beam profile

	param = expt.param;
	monochr = expt.monochr;
	sample = expt.sample;
	beam = expt.beam;
	detector = expt.detector;

	E_exp = expt.E_in;
	data.th_B = bragg_angle(sample, E_exp, sample.hkl);
	% calculate E sampling in eV
	E_min = E_exp - param.E.bound;
	E_max = E_exp + param.E.bound;
	data.E = sequenced(E_min, E_max, param.E.density);
	len_E = length(data.E);
	% calculate th sampling in rad
	th_min = min(bragg_angle(monochr, E_max, monochr.hkl), bragg_angle(sample, E_max, sample.hkl)) - param.th.bound;
	th_max = max(bragg_angle(monochr, E_min, monochr.hkl), bragg_angle(sample, E_min, sample.hkl)) + param.th.bound;
	data.th = sequenced(th_min, th_max, param.th.density);
	data.dth = angular_departure(datat.th, data.th_B);
	len_th = length(data.th);

	data.ib = zeros(len_th, len_E); data.mm = zeros(len_th, len_E); data.rr = zeros(len_th, len_E); data.aa = zeros(len_th, len_E); data.bb = zeros(len_th, len_E); data.cc = zeros(len_th, len_E); data.dd = zeros(len_th, len_E);
	for k=1:len_E
		% Bragg angle in rad
		th_B = bragg_angle(monochr, data.E(k), monochr.hkl);
		% calculate theoretical quantities
		beam.bragg_angle = bragg_angle(monochr, E_exp, monochr.hkl);
		% polarization factors
		beam.diffracted = beam.direction; % NOTE: this may need reconsidering
		beam.reflected = [cos(2.0 .* th_B), 0, sin(2.0 .* th_B)];
		detector.direction = (rotation_matrix(beam.normal, detector.phi) * (rotation_matrix(beam.binormal, -detector.th) * beam.direction'))';
		proj_diffracted_detector = detector.direction * beam.diffracted';
		P0 = [sqrt(1 - tan(th_B).^2.0 .* proj_diffracted_detector), -proj_diffracted_detector];
		Ph = P0;
		% beam before monochromator
		data.ib(:,k) = beam_profile(data.th, data.E, 'gauss', [beam.bragg_angle, beam.divergence]);
		for nu=1:2
			% monochromator transfer function
			data.mm(:,k) = data.mm(:,k)' + beam.polarization(nu) * reflectivity(data.th, data.E(k), monochr, monochr.hkl, nu);
			% sample quantities (I: intensity; xi: amplitude ratio)
			[I, xi] = reflectivity(data.th, data.E(k), sample, sample.hkl, nu);
			% reflectivity 
			data.rr(:,k) = data.rr(:,k)' + beam.polarization(nu) .* I;
			% standing waves (zz: amplitude attenuation; aa, bb: diagonal; cc: non-diagonal real; dd: non-diagonal imag)
			zz = 1e5 .* amplitude_attenuation(data.th, data.E(k), sample, sample.hkl, nu, beam, detector);
			data.aa(:,k) = data.aa(:,k)' + beam.polarization(nu) .* (P0(nu).^2.0) .* ones(size(data.th)) .* zz;
			data.bb(:,k) = data.bb(:,k)' + beam.polarization(nu) .* (Ph(nu).^2.0) .* abs(xi).^2.0 .* zz;
			data.cc(:,k) = data.cc(:,k)' + beam.polarization(nu) .* (P0(nu).*Ph(nu)) .* 2.0 .* real(xi) .* zz; % sign is arbitrary?
			data.dd(:,k) = data.dd(:,k)' + beam.polarization(nu) .* (P0(nu).*Ph(nu)) .* 2.0 .* imag(xi) .* zz; % sign is arbitrary?
		end
	end
  	data.ibmm = data.ib .* data.mm;
	data.exact.aa = iconv2(data.aa', data.ibmm', 'same');
	data.exact.bb = iconv2(data.bb', data.ibmm', 'same');
	data.exact.cc = iconv2(data.cc', data.ibmm', 'same');
	data.exact.dd = iconv2(data.dd', data.ibmm', 'same');
	data.exact.rr = iconv2(data.rr', data.ibmm', 'same');
	data.exact.norm = norm(data.exact.rr);
	data.exact.aa = data.exact.aa ./ data.exact.norm;
	data.exact.bb = data.exact.bb ./ data.exact.norm;
	data.exact.cc = data.exact.cc ./ data.exact.norm;
	data.exact.dd = data.exact.dd ./ data.exact.norm;
	data.exact.rr = data.exact.rr ./ data.exact.norm;
	data.exact.th = data.th;
	data.exact.dth = data.dth;

	E_id = round(len_E / 2.0);
	data.approx.aa = iconv2(data.aa(:,E_id)', data.mm(:,E_id)', 'same');
	data.approx.bb = iconv2(data.bb(:,E_id)', data.mm(:,E_id)', 'same');
	data.approx.cc = iconv2(data.cc(:,E_id)', data.mm(:,E_id)', 'same');
	data.approx.dd = iconv2(data.dd(:,E_id)', data.mm(:,E_id)', 'same');
	data.approx.rr = iconv2(data.rr(:,E_id)', data.mm(:,E_id)', 'same');
	data.approx.norm = norm(data.approx.rr);
	data.approx.aa = data.approx.aa ./ data.approx.norm;
	data.approx.bb = data.approx.bb ./ data.approx.norm;
	data.approx.cc = data.approx.cc ./ data.approx.norm;
	data.approx.dd = data.approx.dd ./ data.approx.norm;
	data.approx.rr = data.approx.rr ./ data.approx.norm;
	data.approx.th = data.th;
	data.approx.dth = data.dth;
	
	% alignment
	data.shift = x_align(data.exact.th, data.exact.rr, data.approx.th, data.approx.rr);
	data.approx.th = data.approx.th + data.shift;
	data.approx.dth = data.approx.dth + data.shift;
end
