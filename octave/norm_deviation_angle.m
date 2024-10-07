function theta = norm_deviation_angle(eta, theta_B, chi_0, chi_h, chi_hb, alpha, polarization)
    % normalized deviation from Bragg angle to angle
    
	C = polarization_factor(nu, theta_B);
	gamma_0 = cos(pi/2 - theta_B + alpha);
	gamma_h = cos(pi/2 + theta_B + alpha);
	gamma = gamma_h ./ gamma_0;
	theta = real((2.0 .* sqrt(chi_h .* chi_hb) .* eta .* abs(C) .* sqrt(abs(gamma)) + chi_0 .* gamma + 2.0 .* theta_B .* sin(2.0 .* theta_B) - chi_0) / (2.0 .* sin(2.0 .* theta_B)));
end
