function eta = norm_deviation(th, th_B, chi_0, chi_h, chi_hb, alpha, polarization)
	% angle to normalized deviation from Bragg angle in complex form (thus including absorption) with asymmetry   

	C = polarization_factor(polarization, th_B);
	gamma_0 = cos(pi/2 - th_B + alpha);
	gamma_h = cos(pi/2 + th_B + alpha);
	gamma = gamma_h ./ gamma_0;
	Delta_theta = angular_departure(th, th_B);
	eta = (2 .* Delta_theta .* sin(2 .* th_B) + chi_0 .* (1 - gamma)) ./ (2 .* sqrt(abs(gamma)) .* abs(C) .* sqrt(chi_h .* chi_hb));
end
