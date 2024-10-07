function d = bragg_distance(crystal, hkl)
    % given lattice parameter, calculate Bragg plane spacing of selected reflection 

	if (hkl == [0 0 0])
		d = 0;
	else
		[aa, bb, cc] = unit_cell_reciprocal_lattice(crystal.cell_size, crystal.cell_angle);
		g = [aa; bb; cc;] ./ (2 * pi);
		d = 1 ./ sqrt(sum((hkl*g).^2.0));
	end
end
