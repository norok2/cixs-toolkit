function plot_compare(expt, groups, workdir)
	% plot results comparison

	% compare different experiments
	for i=1:groups.n
		k = 1;
		figure; clf;
		hold on
		for j=1:length(expt)
			if (expt(j).group == i)
				errorbar(expt(j).result.dE, expt(j).result.ddsf, expt(j).result.D_ddsf, ['>' id2color(k, 'red')]);
				errorbar(expt(j).result.dE, expt(j).result.ndsf, expt(j).result.D_ndsf, ['>' id2color(k+1, 'red')]);
				group_name = expt(j).name(1:6);
				plot_legend{k} = sprintf('%d (dDSF)', fix(k/2)+1);
				plot_legend{k+1} = sprintf('%d (nDSF)', fix(k/2)+1);
				k = k + 2;
			end
		end
		hold off
		legend(plot_legend, 'location', 'northeast');
		xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
		print_plot([workdir group_name '_comp']);
	end
end
