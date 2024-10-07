function data = compute_standing_wave_coeff(data, expt)
	% calculate theoretical values of standing wave parameters

	param = expt.param;
	beam = expt.beam;
	sample = expt.sample;
	monochr = expt.monochr;
	detector = expt.detector;

	len_dE = length(data.val.dE);
	len_dth = length(data.val.dth);

	% calculate theoretical standing wave parameters
	data.val.sw.aa = zeros(size(data.val.ss)); data.val.sw.bb = zeros(size(data.val.ss)); data.val.sw.cc = zeros(size(data.val.ss)); data.val.sw.dd = zeros(size(data.val.ss)); data.val.sw.rr = zeros(size(data.val.ss)); data.val.th_B = zeros(size(data.val.dE));
	for j=1:len_dE
		% energy in eV
		E_exp = data.val.dE(j) + expt.E_in;
		% Bragg angle in rad
		th_B_exp = bragg_angle(sample, E_exp, sample.hkl);
		% calculate theoretical quantities
		if (param.precise_convolution) % NOT WORKING UNLESS MONOCHROMATOR AND SAMPLE HELPER REFLECTION ARE THE SAME
			% calculate E sampling in eV
			E_min = E_exp - param.E.bound;
			E_max = E_exp + param.E.bound;
			E = sequenced(E_min, E_max, param.E.density);
			% calculate th sampling in rad
			th_min = min(bragg_angle(monochr, E_max, monochr.hkl), bragg_angle(sample, E_max, sample.hkl)) + param.dth.min;
			th_max = max(bragg_angle(monochr, E_min, monochr.hkl), bragg_angle(sample, E_min, sample.hkl)) + param.dth.max;
			th = sequenced(th_min, th_max, param.th.density);
			% WARNING: normalization issues?
			ibt = zeros(length(th), length(E)); mmt = zeros(length(th), length(E)); rrt = zeros(length(th), length(E)); aat = zeros(length(th), length(E)); bbt = zeros(length(th), length(E)); cct = zeros(length(th), length(E)); ddt = zeros(length(th), length(E));
			for k=1:length(E)
				% Bragg angle in rad
				th_B = bragg_angle(monochr, E(k), monochr.hkl);
				% polarization factors
				beam.diffracted = beam.direction; % NOTE: this may need reconsidering
				beam.reflected = [cos(2.0 .* th_B), 0, sin(2.0 .* th_B)];
				detector.direction = (rotation_matrix(beam.normal, detector.phi) * (rotation_matrix(beam.binormal, -detector.th) * beam.direction'))';
				proj_diffracted_detector = detector.direction * beam.diffracted';
				P0 = [sqrt(1 - tan(th_B).^2.0 .* proj_diffracted_detector), -proj_diffracted_detector];
				Ph = P0;
				% beam before monochromator
				ibt(:,k) = beam_profile(th, E, 'gauss', [th_B_exp, beam.divergence]);
				for nu=1:2
					% monochromator transfer function
					mmt(:,k) = mmt(:,k)' + beam.polarization(nu) * reflectivity(th, E(k), monochr, monochr.hkl, nu);
					% sample quantities (I: intensity; xi: amplitude ratio)
					[I, xi] = reflectivity(th, E(k), sample, sample.hkl, nu);
					% reflectivity 
					rrt(:,k) = rrt(:,k)' + beam.polarization(nu) .* I;
					% standing waves (zz: amplitude attenuation; aat, bbt: diagonal; cct: non-diagonal)
					zz = 1e5 .* amplitude_attenuation(th, E(k), sample, sample.hkl, nu, beam, detector);
					aat(:,k) = aat(:,k)' + beam.polarization(nu) .* (P0(nu).^2.0) .* ones(size(th)) .* zz;
					bbt(:,k) = bbt(:,k)' + beam.polarization(nu) .* (Ph(nu).^2.0) .* abs(xi).^2.0 .* zz;
					cct(:,k) = cct(:,k)' + beam.polarization(nu) .* (P0(nu).*Ph(nu)) .* 2.0 .* real(xi) .* zz; % sign is arbitrary?
					ddt(:,k) = ddt(:,k)' + beam.polarization(nu) .* (P0(nu).*Ph(nu)) .* 2.0 .* imag(xi) .* zz; % sign is arbitrary?
				end
			end
			% 2D convolution for (E,th) dispersion, integrated over Energy
			mmt = normalize(ibt .* mmt);
			sw.aa = iconv2(aat', mmt', 'same');
			sw.bb = iconv2(bbt', mmt', 'same');
			sw.cc = iconv2(cct', mmt', 'same');
			sw.dd = iconv2(ddt', mmt', 'same');
			sw.rr = iconv2(rrt', mmt', 'same');
		else
			th_min = bragg_angle(sample, E_exp, sample.hkl) + param.dth.min;
			th_max = bragg_angle(sample, E_exp, sample.hkl) + param.dth.max;
			th = sequenced(th_min, th_max, param.th.density);
			thmm_min = bragg_angle(monochr, E_exp, monochr.hkl) + param.dth.min;
			thmm_max = bragg_angle(monochr, E_exp, monochr.hkl) + param.dth.max;
			thmm = sequenced(thmm_min, thmm_max, param.th.density);
			% polarization factors
			beam.diffracted = beam.direction; % NOTE: this may need reconsidering
			beam.reflected = [cos(2.0 .* th_B_exp), 0, sin(2.0 .* th_B_exp)];
			detector.direction = (rotation_matrix(beam.normal, detector.phi) * (rotation_matrix(beam.binormal, -detector.th) * beam.direction'))';
			proj_diffracted_detector = detector.direction * beam.diffracted';
			P0 = [sqrt(1 - tan(th_B_exp).^2.0 .* proj_diffracted_detector), -proj_diffracted_detector];
			Ph = P0;
			mmt = zeros(size(thmm)); rrt = zeros(size(th)); aat = zeros(size(th)); bbt = zeros(size(th)); cct = zeros(size(th)); ddt = zeros(size(th));
			for nu=1:2
				% monochromator transfer function
				mmt = mmt + beam.polarization(nu) * reflectivity(thmm, E_exp, monochr, monochr.hkl, nu);
				% sample quantities (I: intensity; xi: amplitude ratio)
				[I, xi] = reflectivity(th, E_exp, sample, sample.hkl, nu);
				% reflectivity 
				rrt = rrt + beam.polarization(nu) .* I;
				% standing waves (zz: amplitude attenuation; aat, bbt: diagonal; cct: non-diagonal)
				zz = 1e5 .* amplitude_attenuation(th, E_exp, sample, sample.hkl, nu, beam, detector);
				aat = aat + beam.polarization(nu) .* (P0(nu).^2.0) .* ones(size(th)) .* zz;
				bbt = bbt + beam.polarization(nu) .* (Ph(nu).^2.0) .* abs(xi).^2.0 .* zz;
				cct = cct + beam.polarization(nu) .* (P0(nu).*Ph(nu)) .* 2.0 .* real(xi) .* zz; % sign is arbitrary?
				ddt = ddt + beam.polarization(nu) .* (P0(nu).*Ph(nu)) .* 2.0 .* imag(xi) .* zz; % sign is arbitrary?
			end
			% convolution for energy dispersion projected over th
			mmt = normalize(mmt); % should solve normalization issues
			sw.aa = iconv2(aat, mmt, 'same');
			sw.bb = iconv2(bbt, mmt, 'same');
			sw.cc = iconv2(cct, mmt, 'same');
			sw.dd = iconv2(ddt, mmt, 'same');
			sw.rr = iconv2(rrt, mmt, 'same');
		end

		% 3D grid construction
		data.val.th_B(j) = th_B_exp;
		data.val.sw.aa(j,:) = interp1(th, sw.aa, data.val.dth + data.val.th_B(j), 'linear');
		data.val.sw.bb(j,:) = interp1(th, sw.bb, data.val.dth + data.val.th_B(j), 'linear');
		data.val.sw.cc(j,:) = interp1(th, sw.cc, data.val.dth + data.val.th_B(j), 'linear');
		data.val.sw.dd(j,:) = interp1(th, sw.dd, data.val.dth + data.val.th_B(j), 'linear');
		data.val.sw.rr(j,:) = interp1(th, sw.rr, data.val.dth + data.val.th_B(j), 'linear');
	end
end
