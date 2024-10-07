function [I, xi, eta] = reflectivity3(theta, E, crystal, hkl, polarization)
	% reflectivity of the Bragg case of diffraction in dynamical theory of diffraction (Matsushita)

	% Bragg angle in rad
	theta_B = bragg_angle(crystal, E, hkl);
	% asymmetry correction
	b = sin(theta_B + crystal.asymmetry) ./ sin(theta_B - crystal.asymmetry);
	% angle variation
	Delta_theta = (theta - theta_B);
	
	% sign of susceptivity to form factor conversion seems incorrect in Matsushita
	% structure form factor in forward direction
	[F_0, chi_0] = structure_form_factor(crystal, E, [0 0 0], +1, +1, -1);
	% structure form factor in [hkl] direction
	[F_h, chi_hh] = structure_form_factor(crystal, E, hkl, +1, +1, -1);
	% structure form factor in -[hkl] direction
	[F_hb, chi_hb] = structure_form_factor(crystal, E, -hkl, +1, +1, -1);
	% structure form factor in -[hkl] direction (centrosymmetric)
	chi_h = sqrt(chi_hh .* chi_hb);
	
	C = polarization_factor(polarization, theta_B);
	W = (((1 + abs(b)) .* real(chi_0)) + (2.0 .* abs(b) .* sin(2 .* theta_B) .* Delta_theta)) ./ (2 .* C * sqrt(b) .* abs(real(chi_h)));
	g = ((1 + abs(b)) .* imag(chi_0)) ./ (2.0 .* C * sqrt(b) .* abs(real(chi_h)));
	kappa = imag(chi_h) ./ real(chi_h);
	L = (W.^2.0 + g.^2.0 + sqrt((W.^2.0 - g.^2.0 - 1.0 + kappa.^2.0).^2.0 + 4.0 .* (g .* W - kappa).^2.0)) ./ (1.0 + kappa.^2.0);
	I = L - sqrt(L.^2.0 - 1);
	% NOTE: from Authier
	eta = W + 1i .* g;
	xi = (sqrt(abs(b))) .* (eta - sign(W) .* sqrt(eta.^2 - 1)) .* exp(1i .* (atan2(imag(chi_h), real(chi_h))));
end