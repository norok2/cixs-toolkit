function data = remove_elastic_peak(data)
	% subtract elastic peak (using Gaussian profile) from val

	% fit elastic peak
	data.val.ep = zeros(size(data.val.ss_nobg));
	elastic_E = zeros(size(data.val.dth_raw));
	elastic_width = zeros(size(data.val.dth_raw));
	fit_range = (data.val.dE_raw < 2 * abs(min(data.val.dE_raw)));
	x = data.val.dE_raw(fit_range);
	for j=1:length(data.val.dth_raw)
		% calculate elastic peak
		p0 = [0.1, 0.5, 0.5]; % best
		y = data.val.ss_nobg(fit_range,j)';
		f = @(x,p) (p(3) .* f_gauss(x, p(1), p(2)));
		p = datafit(x, y, p0, f, 'fminsearch'); % NOTE: this may affect fitting
		elastic_E(j) = p(1);
		elastic_width(j) = p(2) * (sqrt(8.0*log(2.0)));
		data.val.ep(:,j) = f(data.val.dE_raw, p);
		% subtract elastic peak (assume errorless value of elastic peak)
		data.val.ss_noep(:,j) = data.val.ss_nobg(:,j) - data.val.ep(:,j);
		data.val.D_ss_noep(:,j) = data.val.D_ss_nobg(:,j);
	end
	data.val.elastic_E = mean(elastic_E);
	data.val.elastic_width = mean(elastic_width);
end
