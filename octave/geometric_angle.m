function th_o = geometric_angle(th_i, th_B, phi, beam)
	% convert angle to account for sample rotation around beam normal
	
	th_o = th_i;
	for i=1:length(th_o)
		% Method 1 (WRONG?)
%  		th_o(i) = atan(tan(th_i(i) - th_B) .* cos(phi)) + th_B;
%  		th_o(i) = atan(tan(th_i(i)) .* cos(phi));

%  		% Method 2 (WRONG?)
%  		surface_normal = (rotation_matrix(beam.binormal, -th_i(i)) * beam.normal')';
%  		reflection.direction = (rotation_matrix(surface_normal, -th_i(i)) * (rotation_matrix(beam.normal, phi) * beam.direction'))';
%  		reflection.projection = plane_projection(reflection.direction, beam.direction, beam.normal);
%  		th_o(i) = acos(reflection.projection * reflection.direction') ./ norm(reflection.projection);

%  		% Method 3 (WRONG?)
%  		reflection.rotation_axis = (rotation_matrix(beam.normal, phi) * beam.binormal')';
%  		reflection.direction = (rotation_matrix(beam.binormal, -2*th_B) * beam.direction')';
%  		reflection.distorted = (rotation_matrix(reflection.rotation_axis, -th_i(i)) * (rotation_matrix(beam.normal, phi) * reflection.direction'))';
%  		reflection.projection = plane_projection(reflection.distorted, beam.direction, beam.normal);
%  		th_o(i) = acos(normalize(reflection.projection) * reflection.direction');

%  		% Method 4 (almost correct?)
%  		reflection.rotation_axis = (rotation_matrix(beam.normal, phi) * beam.binormal')';
%  		reflection.distorted = (rotation_matrix(reflection.rotation_axis, -2*th_i(i)) * (rotation_matrix(beam.normal, phi) * beam.direction'))';
%  		reflection.projection = plane_projection(reflection.distorted, beam.direction, beam.normal);
%  		th_o(i) = acos(normalize(reflection.projection) * beam.direction') / 2;		

		% Method 5 (correct?)
		reflection.rotation_axis = (rotation_matrix(beam.normal, phi) * beam.binormal')';
		reflection.direction = (rotation_matrix(beam.binormal, -2*th_B) * beam.direction')';
		reflection.distorted = (rotation_matrix(reflection.rotation_axis, -2*(th_i(i) - th_B)) * rotation_matrix(beam.normal, phi) * rotation_matrix(beam.binormal, -2*th_B) * beam.direction')';
		reflection.projection = normalize(plane_projection(reflection.distorted, beam.direction, beam.normal));
		th_o(i) = acos(reflection.projection * beam.direction') / 2;
	end
end