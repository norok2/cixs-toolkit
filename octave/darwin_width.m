function darwin = darwin_width(E, crystal, hkl, polarization)
	% Darwin width, corresponding to the width of the total reflection domain in Bragg geometry

	% Bragg angle in rad
	th_B = bragg_angle(crystal, E, hkl);
	% structure form factor in [hkl] direction (fourier: -1, im_ff: +1, suscept: -1)
	[F_h, chi_h] = structure_form_factor(crystal, E, hkl);
	% structure form factor in -[hkl] direction
	[F_hb, chi_hb] = structure_form_factor(crystal, E, -hkl);
	
	
	darwin = zeros(size(E));
	for nu=1:length(polarization)
		C = polarization_factor(nu, th_B);
		psi_0 = (pi/2 - th_B + crystal.asymmetry);
		psi_h = (pi/2 + th_B + crystal.asymmetry);
		gamma_0 = cos(psi_0);
		gamma_h = cos(psi_h);
		gamma = gamma_h ./ gamma_0;
		darwin = darwin + polarization(nu) .* 2.0 .* real((C .* sqrt(abs(gamma) .* chi_h .* chi_hb)) ./ sin(2 .* th_B));
	end
end