function cdata = combine_spectra(expt, groups)
	% combine spectra

	cdata.name = '';
	cdata.q0 = [];
	cdata.qh = [];
	cdata.hkl = [];
	cdata.dE = [];
	cdata.ddsf = [];
	cdata.D_ddsf = [];
	cdata.ndsf = [];
	cdata.D_ndsf = [];
	cdata.rho_ratio = [];
	cdata.D_rho_ratio = [];
	cdata = repmat(cdata, 1, groups.n);

	for i=1:groups.n
		n = groups.dim(i);
		% find first expt of each group
		j = 1;
		while (expt(j).group ~= i)
			j = j + 1;
		end
		% initialize combined spectra
		cdata(i).name = expt(j).name(1:6);
		cdata(i).q0 = expt(j).sample.q;
		cdata(i).qh = expt(j).sample.q + expt(j).sample.hkl;
		cdata(i).hkl = expt(j).sample.hkl;
		% determine energy ranges
		dE_min = expt(j).result.dE(1);
		dE_max = expt(j).result.dE(end);
		dE_step = min(abs(diff(expt(j).result.dE)));
		for j=1:length(expt)
			if (expt(j).group == i)
				dE_min = max(expt(j).result.dE(1), dE_min);
				dE_max = min(expt(j).result.dE(end), dE_max);
				dE_step = min(min(abs(diff(expt(j).result.dE))), dE_step);
			end
		end
		dE = dE_min:(dE_step / 2):dE_max;
		% calculate average of spectra and error
		ddsf_M = zeros(n, length(dE));
		D_ddsf_M = zeros(n, length(dE));
		ndsf_M = zeros(n, length(dE));
		D_ndsf_M = zeros(n, length(dE));
		k = 1;
		for j=1:length(expt)
			if (expt(j).group == i)
				ddsf_M(k,:) = interp1(expt(j).result.dE, expt(j).result.ddsf, dE);
				D_ddsf_M(k,:) = interp1(expt(j).result.dE, expt(j).result.D_ddsf, dE);
				ndsf_M(k,:) = interp1(expt(j).result.dE, expt(j).result.ndsf, dE);
				D_ndsf_M(k,:) = interp1(expt(j).result.dE, expt(j).result.D_ndsf, dE);
				k = k + 1;
			end
		end
		ddsf = mean(ddsf_M);
		D_ddsf = std(ddsf_M) / n; % mean(D_ddsf_M);
		ndsf = mean(ndsf_M);
		D_ndsf = std(ndsf_M) / n; % mean(D_ndsf_M);
		% calculate rho ratio
		rho_ratio = (trapz(dE, dE .* ndsf) / trapz(dE, dE .* ddsf)) * ((cdata(i).q0 * cdata(i).q0') / (cdata(i).q0 * cdata(i).qh'));
		D_rho_ratio = rho_ratio * abs(trapz(dE, dE .* D_ndsf) / trapz(dE, dE .* ndsf)) + abs(trapz(dE, dE .* D_ddsf) / trapz(dE, dE .* ddsf));
		% store results in cdata
		cdata(i).dE = dE;
		cdata(i).ddsf = ddsf;
		cdata(i).D_ddsf = D_ddsf;
		cdata(i).ndsf = ndsf;
		cdata(i).D_ndsf = D_ndsf;
		cdata(i).rho_ratio = rho_ratio;
		cdata(i).D_rho_ratio = D_rho_ratio;
	end
end
