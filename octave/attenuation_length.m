function mu = attenuation_length(crystal, E, data_directory)
    % extract experimental attenuation length of Silicon (Si) in m from data file (source: CXRO)
	% http://henke.lbl.gov/optical_constants/atten2.html

	if (nargin < 3)
		data_directory.attenuation_length = [pwd '/../attenuation_length/'];
	end
	
	M = load([data_directory.attenuation_length crystal.symbol '.txt']);
	mu = interp1(M(:,1), M(:,2), E) .* 1e-6 ./ (crystal.density ./ 1e3);
end

