function FF = atomic_form_factor(elem_id, E, Q, im_form_factor_sign, data_directory)
    % experimental atomic form factor of given element as a function of E and Q
	
	if (nargin < 4)
		im_form_factor_sign = +1;		
	end
	if (nargin < 5)
		data_directory.atomic_form_factor = [pwd '/../atomic_form_factor/'];
		data_directory.atomic_form_factor_momentum = [pwd '/../atomic_form_factor_momentum/'];
	end

	elem = chemical_properties(elem_id);
	
	% estimated Silicon (Si) crystal structure factor in function of incident energy and angular direction
	% Als-Nielsen and McMorrow, Elements of Modern X-Ray Physics (John Wiley & Sons, 2001) p. 115
	FF_Q_param = load([data_directory.atomic_form_factor_momentum elem.symbol '.txt']);
	FF_Q = 0;
	n = ((length(FF_Q_param) - 1) / 2);
	for i=1:n
		FF_Q = FF_Q + (FF_Q_param(i) .* exp(-FF_Q_param(i+n+1) .* (Q./4.0./pi./1e10).^2.0));
	end
	FF_Q = FF_Q + FF_Q_param(n+1);
	
    % extract experimental complex atomic form factor from data file (source: CXRO)
	% http://henke.lbl.gov/optical_constants/asf.html
	FF_E = load([data_directory.atomic_form_factor elem.symbol '.txt']);
	FFr = interp1(FF_E(:,1), FF_E(:,2), E);
	FFi = interp1(FF_E(:,1), FF_E(:,3), E);
	FF_a = (FFr + 1i .* im_form_factor_sign .* FFi) - elem.Z;
	
	FF = FF_Q + FF_a;
end

