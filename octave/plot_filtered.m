function plot_filtered(cdata, workdir)
	% plot results comparison

	% compare different experiments
	for i=1:length(cdata)
		figure; clf;
		hold on
			errorbar(cdata(i).dE, cdata(i).ddsf, cdata(i).D_ddsf, ['>' id2color(1)]);
			errorbar(cdata(i).dE, cdata(i).ndsf, cdata(i).D_ndsf, ['>' id2color(2)]);
		hold off
		xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
		print_plot([workdir cdata(i).name '_filter']);
	end
end
