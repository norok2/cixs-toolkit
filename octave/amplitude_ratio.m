function [xi, eta] = amplitude_ratio(th, E, crystal, hkl, polarization)
    % intensity and amplitude ratio between forward diffracted and Bragg reflected wavefields in Bragg geometry using dynamical theory of diffraction
	
	% Bragg angle in rad
	th_B = bragg_angle(crystal, E, hkl);
	% structure form factor in forward direction (fourier: -1, im_ff: +1, suscept: -1)
	[F_0, chi_0] = structure_form_factor(crystal, E, [0 0 0]);
	% structure form factor in [hkl] direction
	[F_h, chi_h] = structure_form_factor(crystal, E, hkl);
	% structure form factor in -[hkl] direction
	[F_hb, chi_hb] = structure_form_factor(crystal, E, -hkl);
	
	C = polarization_factor(polarization, th_B);
	gamma_0 = cos(pi/2 - th_B + crystal.asymmetry);
	gamma_h = cos(pi/2 + th_B + crystal.asymmetry);
	gamma = gamma_h ./ gamma_0;
	eta = norm_deviation(th, th_B, chi_0, chi_h, chi_hb, crystal.asymmetry, polarization);
	u = -sign(real(eta)); 
	xi = (sign(C) .* sign(gamma_h) ./ sqrt(abs(gamma))) .* (sqrt(chi_h .* chi_hb) ./ chi_hb) .* (eta + u .* sqrt(eta.^2.0 + sign(gamma_h)));
	xi = xi .* exp(1i .* pi); % phase shift alignment
end
