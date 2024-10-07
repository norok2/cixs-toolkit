function plot_plasmons_expt(expt, groups, workdir)
	% plot plasmons for each experiment

	% plot plasmons for each experiment
	for i=1:groups.n
		figure; clf;
		hold on
		for j=1:length(expt)
			if (expt(j).group == i)
				errorbar(expt(j).result.dE, (expt(j).result.ddsf + expt(j).result.ndsf) ./ 2.0, (expt(j).result.D_ddsf + expt(j).result.D_ndsf) ./ 2.0, ['>' id2color(j)]);
				errorbar(expt(j).result.dE, (expt(j).result.ddsf - expt(j).result.ndsf) ./ 2.0, (expt(j).result.D_ddsf + expt(j).result.D_ndsf) ./ 2.0, ['>' id2color(j)]);
				group_name = expt(j).name(1:6);
			end
		end
		hold off
		xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
		print_plot([workdir group_name '_plasmons']);
	end
end
