function [aa, bb, cc] = unit_cell_reciprocal_lattice(size, angle)
	% real space lattice vectors [a, b, c] in cartesian coordinates (with a parallel to X and b in the XY plane)

	[a, b, c] = unit_cell_direct_lattice(size, angle);
	V = (a * cross(b, c)');
	aa = (2.0 * pi / V) .* cross(b, c);
	bb = (2.0 * pi / V) .* cross(c, a);
	cc = (2.0 * pi / V) .* cross(a, b);
end