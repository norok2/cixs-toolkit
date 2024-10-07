function [I, xi, eta] = reflectivity4(theta, E, crystal, hkl)
	% reflectivity of the Bragg case of diffraction in dynamical theory of diffraction (Als-Nielsen McMorrow)

	% classical electron radius (r_e = e^2 / m_e / c^2) in m
	[r_e, UM_r_e, D_r_e] = physical_constant('CLASSICAL_ELECTRON_RADIUS');

	% Bragg angle in rad
	theta_B = bragg_angle(crystal, E, hkl);
	% departure parametre
	zeta = (theta - theta_B) .* cot(theta_B);
	% asymmetry correction
	b = sin(theta_B + crystal.asymmetry) ./ sin(theta_B - crystal.asymmetry);
	zeta = zeta ./ sqrt(b);
	
	% structure form factors
	% Als-Nielsen and McMorrow, Elements of Modern X-Ray Physics (John Wiley & Sons, 2001) p. 193
	%F_0 = 8 * (14 + 0.25 - 1i * 0.33);
	%F_h = 4 * abs(1 - 1i) * (10.54 + 0.25 - 1i * 0.33);
	% structure form factor in forward direction (fourier: -1, im_ff: +1, suscept: -1)
	F_0 = structure_form_factor(crystal, E, [0 0 0], +1, -1, +1);
	% structure form factor in [hkl] direction
	F_hh = structure_form_factor(crystal, E, hkl, +1, -1, +1);
	% structure form factor in -[hkl] direction
	F_hb = structure_form_factor(crystal, E, -hkl, +1, -1, +1);
	% structure form factor in [hkl] direction (centrosymmetric)
	% conjugation can fix +1 imaginary part atomic form factor correction
    F_h = sqrt(F_hh .* F_hb);
	
	V = unit_cell_volume(crystal.cell_size, crystal.cell_angle); % in m
	d = bragg_distance(crystal, hkl); % in m
	m = bragg_order(crystal, hkl);
	L = (m .* V) ./ (2.0 .* d.^2.0 .* r_e); 
	g_0 = F_0 ./ L;
	g_h = F_h ./ L;
	eta = (m .* pi .* zeta ./ g_h - g_0 ./ g_h);
	
	% formula of (S_0/T_0)
	xi = zeros(length(eta),1);
	for i=1:length(eta)
		if (real(eta(i)) > 1.0)
			xi(i) = eta(i) - sqrt((eta(i)^2.0 - 1.0));
		elseif (real(eta(i)) < -1.0)
			xi(i) = eta(i) + sqrt((eta(i)^2.0 - 1.0));
		else %if ((wr(i) < 1.0) && (wr(i) > -1.0))
			xi(i) = eta(i) - 1i * sqrt(1.0 -(eta(i)^2.0));
		end
	end
	
	% 1-line (almost) equivalent
	%r = eta - sign(real(eta)) .* sqrt(eta.^2.0 - 1.0);	
	
	I = abs(xi).^2.0;
end