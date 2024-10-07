function fdata = filter_spectra(cdata, n_point)
	% filter spectra

	fdata = cdata;
	for i = 1:length(cdata)
		x = -length(cdata(i).dE):length(cdata(i).dE);
		H = f_gauss(x, 0, n_point);
%  		fdata(i).ddsf = conv2(cdata(i).ddsf, H, 'same');
%  		fdata(i).D_ddsf = max(abs(fdata(i).ddsf - cdata(i).ddsf), cdata(i).D_ddsf);
		fdata(i).ndsf = conv2(cdata(i).ndsf, H, 'same');
		fdata(i).D_ndsf = max(abs(fdata(i).ndsf - cdata(i).ndsf), cdata(i).D_ndsf);
		fdata(i).rho_ratio = (trapz(fdata(i).dE, fdata(i).dE .* fdata(i).ndsf) / trapz(fdata(i).dE, fdata(i).dE .* fdata(i).ddsf)) * ((fdata(i).q0 * fdata(i).q0') / (fdata(i).q0 * fdata(i).qh'));
		fdata(i).D_rho_ratio = fdata(i).rho_ratio * abs(trapz(fdata(i).dE, fdata(i).dE .* fdata(i).D_ndsf) / trapz(fdata(i).dE, fdata(i).dE .* fdata(i).ndsf)) + abs(trapz(fdata(i).dE, fdata(i).dE .* fdata(i).D_ddsf) / trapz(fdata(i).dE, fdata(i).dE .* fdata(i).ddsf));	
	end
end
