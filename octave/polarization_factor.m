function C = polarization_factor(polarization, theta_B)
	switch polarization
		case {1, 'sigma', 's', 'horizonta', 'h', 'H', 'hor'}
			C = 1;
		case {2, 'pi', 'p', 'vertical', 'v', 'V', 'ver'}
			C = cos(2 * theta_B);
		otherwise
			error('unknown polarization state');
	end
end