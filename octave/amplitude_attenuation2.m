function [z, mu, mu_e] = amplitude_attenuation2(th, E, crystal, hkl, polarization, beam, detector)
	% depth integral of standing wave field damped by absorption (including abnormal correction)

	% Bragg angle in rad
	th_B = bragg_angle(crystal, E, hkl);
	% structure form factor in forward direction
	[F_0, chi_0] = structure_form_factor(crystal, E, [0 0 0]);
	% structure form factor in [hkl] direction
	[F_h, chi_h] = structure_form_factor(crystal, E, hkl);
	% structure form factor in -[hkl] direction
	[F_hb, chi_hb] = structure_form_factor(crystal, E, -hkl);
	
	mu_0 = 1 ./ attenuation_length(crystal, E);
	mu_1 = 1 ./ attenuation_length(crystal, E);
	surface_normal = (rotation_matrix([0, 1, 0], crystal.asymmetry) * rotation_matrix([0, 1, 0], -th_B) * beam.normal')';
	detector.direction = (rotation_matrix(beam.normal, detector.phi) * (rotation_matrix(beam.binormal, -detector.th) * beam.direction'))';
	zeta = acos(detector.direction * surface_normal'); % NOTE: this value should be positive

	C = polarization_factor(polarization, th_B);
	psi_0 = (pi/2 - th + crystal.asymmetry);
	psi_h = (pi/2 - th - crystal.asymmetry);
	psi_0B = (pi/2 - th_B + crystal.asymmetry);
	gamma_0 = cos(psi_0);
	gamma_h = cos(psi_h);
	beta_r = (2 .* (psi_0 - psi_0B) .* sin(2 .* th_B) + abs(real(chi_0)) .* (1 + abs(gamma_h)./gamma_0));
	beta_i = abs(imag(chi_0)) .* (1 + abs(gamma_h)./gamma_0);
	Psi_h = real(chi_h) .* imag(chi_hb) + real(chi_hb) .* imag(chi_h);
	W = sqrt(beta_r.^2 - beta_i.^2 - 4 .* C.^2.0 .* abs(real(chi_h)).^2 .* abs(gamma_h) ./ gamma_0 + 2i .* (beta_r .* beta_i - 2.0 .* C.^2.0 .* Psi_h .* abs(gamma_h) ./ gamma_0));
	u = sign(imag(W)); % NOTE: +- sign not specified in article
	sigma = (mu_0 ./ 2.0) .* ((1 ./ gamma_0) - (1 ./ abs(gamma_h)) + u .* (imag(W) ./ abs(imag(chi_0)) ./ abs(gamma_h)));
	mu_e = sigma .* gamma_0;
	mu = (sigma  + mu_1 ./ cos(zeta));
		
	z = (1.0 - exp(-mu .* crystal.thickness)) ./ mu;
end
