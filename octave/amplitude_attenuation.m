function [z, mu, mu_e] = amplitude_attenuation(th, E, crystal, hkl, polarization, beam, detector)
	% depth integral of standing wave field damped by absorption (including abnormal correction)

	% Bragg angle in rad
	th_B = bragg_angle(crystal, E, hkl);
	% structure form factor in forward direction
	[F_0, chi_0] = structure_form_factor(crystal, E, [0 0 0]);
	% structure form factor in [hkl] direction
	[F_h, chi_h] = structure_form_factor(crystal, E, hkl);
	% structure form factor in -[hkl] direction
	[F_hb, chi_hb] = structure_form_factor(crystal, E, -hkl);
	% wavelength in m
	lambda = conv_E('eV', 'm', E);	
	
	%mu_0 = 1 ./ attenuation_length(crystal, E);
	mu_0 = -2.0 * pi ./ lambda .* imag(chi_0);
	mu_1 = 1 ./ attenuation_length(crystal, E);
	surface_normal = (rotation_matrix([0, 1, 0], crystal.asymmetry) * rotation_matrix([0, 1, 0], -th_B) * beam.normal')';
	detector.direction = (rotation_matrix(beam.normal, detector.phi) * (rotation_matrix(beam.binormal, -detector.th) * beam.direction'))';
	zeta = acos(detector.direction * surface_normal'); % NOTE: this value should be positive and smaller than PI/2

	C = polarization_factor(polarization, th_B);
	gamma_0 = cos(pi/2 - th_B + crystal.asymmetry);
	gamma_h = cos(pi/2 + th_B + crystal.asymmetry);
	eta = norm_deviation(th, th_B, chi_0, chi_h, chi_hb, crystal.asymmetry, polarization);	
	Lambda_0 = lambda .* sqrt(gamma_0 .* abs(gamma_h)) ./ (abs(C) .* sqrt(chi_h .* chi_hb));
	u = -sign(real(eta)); % NOTE: using u = +1 gives non-continuos function z
	W = (eta + u .* sqrt(eta.^2.0 + sign(gamma_h))) ./ Lambda_0;
	mu_e = mu_0 + 2.0 .* pi .* gamma_0 .* imag(W);
	mu = (mu_e ./ gamma_0 + mu_1 ./ cos(zeta));
	
	% \int_0^L {\exp{-mu z}}{dz} = \frac{1 - \exp{-mu L}}{mu}
	z = (1.0 - exp(-mu .* crystal.thickness)) ./ mu;
end
