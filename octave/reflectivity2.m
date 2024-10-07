function [I, xi, eta] = reflectivity2(theta, E, crystal, hkl, polarization)
	% reflectivity of the Bragg case of diffraction in dynamical theory of diffraction (Schuelke)

	% Bragg angle in rad
	theta_B = bragg_angle(crystal, E, hkl);
	
	% structure form factor in forward direction (fourier: -1, im_ff: +1, suscept: +1)
	[F_0, chi_0] = structure_form_factor(crystal, E, [0 0 0], -1, +1, -1);
	% structure form factor in [hkl] direction
	[F_h, chi_hh] = structure_form_factor(crystal, E, hkl, -1, +1, -1);
	% structure form factor in -[hkl] direction
	[F_hb, chi_hb] = structure_form_factor(crystal, E, -hkl, -1, +1, -1);
	% structure form factor in -[hkl] direction (centrosymmetric)
	chi_h = sqrt(chi_hh .* chi_hb);
	
	C = polarization_factor(polarization, theta_B);
	
	% NOTE: tau dependence dropped because it was used only to compute sign of C
 	psi_0 = (pi/2 - theta + crystal.asymmetry);
 	psi_h = (pi/2 - theta - crystal.asymmetry);
	psi_0B = (pi/2 - theta_B + crystal.asymmetry);
	gamma_0 = cos(psi_0);
	gamma_h = cos(psi_h);
	v = (sqrt(gamma_0./abs(gamma_h)) + sqrt(abs(gamma_h)./gamma_0))./2;
	%B = 1 ./ sqrt(abs(gamma_h) ./ gamma_0);
	beta_r = (2 .* (psi_0 - psi_0B) .* sin(2 .* theta_B) + abs(real(chi_0)) .* (1 + abs(gamma_h)./gamma_0));
	%beta_i = abs(imag(chi_0)) .* (1 + abs(gamma_h)./gamma_0);
	Phi_h = real(chi_h) .* real(chi_hb) - imag(chi_h) .* imag(chi_hb);
	Psi_h = real(chi_h) .* imag(chi_hb) + real(chi_hb) .* imag(chi_h);
	y = beta_r ./ (2 .* C .* sqrt(abs(gamma_h)./gamma_0) .* sqrt(Phi_h));
	d = abs(imag(chi_0)) .* Psi_h .* v ./ sqrt(Phi_h) ./ Phi_h ./ abs(C) ./ 2.0 + y;
	f = abs(imag(chi_0)) .* v ./ abs(C) ./ sqrt(Phi_h) - y .* Psi_h ./ Phi_h ./ 2.0;
	E = d.^2 + f.^2 + sqrt(1 + (d.^2 + f.^2).^2.0 - 2.*(d.^2 - f.^2));
	eta = y;
	xi = -sign(C) .* sqrt(gamma_0 ./ abs(gamma_h)) .* sqrt(chi_h ./ chi_hb) .* (d .* (1 - sqrt((E - 1) ./ (E + 1))) + 1i .* f .* (1 - sqrt((E + 1) ./ (E - 1))));
	I = abs(xi).^2.0 .* abs(gamma_h ./ gamma_0); % NOTE: this is borrowed from Authier	
end