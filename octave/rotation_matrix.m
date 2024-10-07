function R = rotation_matrix(v, theta)
	% generate rotation matrix of angle theta around a 3D vector v
	
	if (([1, 3] == size(v)))
		v = v';
	end
	
	I = eye(3);
	V = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
	VV = v * v';
	R = I .* cos(theta) + V .* sin(theta) + VV .* (1 - cos(theta));
end