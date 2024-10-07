function [I, xi, eta] = reflectivity(th, E, crystal, hkl, polarization)
    % reflectivity of the Bragg case of diffraction in dynamical theory of diffraction using Laue approach.
	
	% Bragg angle in rad
	th_B = bragg_angle(crystal, E, hkl);
	
	gamma_0 = cos(pi/2 - th_B + crystal.asymmetry);
	gamma_h = cos(pi/2 + th_B + crystal.asymmetry);
	gamma = gamma_h ./ gamma_0;
	[xi, eta] = amplitude_ratio(th, E, crystal, hkl, polarization);
	I = abs(xi).^2.0 .* abs(gamma);
end
