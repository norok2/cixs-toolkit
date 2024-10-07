function V = unit_cell_volume(size, angle)
	% calculate volume of lattice unit cell (size and angle are both 3D vectors)
	
	V = size(1) * size(2) * size(3) * sqrt(1.0 - cos(angle(1))^2.0 - cos(angle(2))^2.0 - cos(angle(3))^2 + 2 * cos(angle(1)) * cos(angle(2)) * cos(angle(3)));
end