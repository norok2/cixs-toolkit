function a = color_list(method)
	% array of available colors for plotting

	if (nargin < 1)
		method = 'full';
	end

	switch (method)
		case {'full'}
%  			a = ['r', 'g', 'b', 'c', 'm', 'y', 'k']'; % Matlab
%  			a = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w']'; % Octave
			a = ['b', 'r', 'g', 'c', 'm', 'y', 'k']'; % colors in common to matlab and octave
		case {'reduced', 'red'}
			a = ['b', 'r', 'g', 'm']';
		case {'rgb'}
			a = ['b', 'r', 'g']';
		case {'cmyk'}
			a = ['c', 'm', 'y', 'k']';
		otherwise
			error('invalid color option.');
end
