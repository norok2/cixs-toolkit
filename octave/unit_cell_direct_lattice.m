function [a, b, c] = unit_cell_direct_lattice(size, angle)
	% direct space lattice vectors [a, b, c] in cartesian coordinates (with a parallel to X and b in the XY plane)

	a = size(1) .* [1.0, 0.0, 0.0];
	b = size(2) .* [cos(angle(3)), sin(angle(3)), 0.0];
	c = size(3) .* [cos(angle(2)) * sin(angle(3)), cos(angle(1)) - cos(angle(2)) * cos(angle(3)), sqrt(1.0 - cos(angle(1))^2.0 - cos(angle(2))^2.0 - cos(angle(3))^2 + 2 * cos(angle(1)) * cos(angle(2)) * cos(angle(3)))] ./ sin(angle(3));
end
