function p = plane_projection(v, p1, p2)
	pn = cross(p1, p2) ./ norm(p1) ./ norm(p2);
	p = v - ((pn' * pn) * v')';
end