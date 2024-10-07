function plot_literature(cdata, ldata, workdir)
	% plot comparison with literature

	for i=1:length(ldata)
		if (ldata(i).available == true)
			figure; clf;
			hold on
				errorbar(cdata(i).dE, cdata(i).ddsf, cdata(i).D_ddsf, ['>' id2color(1, 'red')]);
				errorbar(cdata(i).dE, cdata(i).ndsf, cdata(i).D_ndsf, ['>' id2color(2, 'red')]);
				plot(ldata(i).data(:,1), ldata(i).data(:,2) .* ldata(i).norm, ['-o' id2color(3, 'red')]);
				plot(ldata(i).data(:,1), ldata(i).data(:,3) .* ldata(i).norm, ['-o' id2color(4, 'red')]);
			hold off
			legend('dDSF', 'nDSF', 'dDSF lit.', 'nDSF lit.');
			xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
			print_plot([workdir ldata(i).name '_lit_comp']);

			figure; clf;
			hold on
				errorbar(cdata(i).dE, cdata(i).ddsf, cdata(i).D_ddsf, ['>' id2color(1, 'red')]);
				errorbar(cdata(i).dE, cdata(i).ndsf, cdata(i).D_ndsf, ['>' id2color(2, 'red')]);
				plot(ldata(i).data(:,1) .* ldata(i).magnify, ldata(i).data(:,2) .* ldata(i).norm, ['-o' id2color(3, 'red')]);
				plot(ldata(i).data(:,1) .* ldata(i).magnify, ldata(i).data(:,3) .* ldata(i).norm, ['-o' id2color(4, 'red')]);
			hold off
			legend('dDSF', 'nDSF', 'dDSF lit.', 'nDSF lit.');
			xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
			print_plot([workdir ldata(i).name '_lit_comp_magn']);

			figure; clf;
			hold on
				errorbar(cdata(i).dE, cdata(i).ddsf, cdata(i).D_ddsf, ['>' id2color(1, 'red')]);
				errorbar(cdata(i).dE, cdata(i).ndsf, cdata(i).D_ndsf, ['>' id2color(2, 'red')]);
				plot(ldata(i).data(:,1) + ldata(i).shift, ldata(i).data(:,2) .* ldata(i).norm, ['-o' id2color(3, 'red')]);
				plot(ldata(i).data(:,1) + ldata(i).shift, ldata(i).data(:,3) .* ldata(i).norm, ['-o' id2color(4, 'red')]);
			hold off
			legend('dDSF', 'nDSF', 'dDSF lit.', 'nDSF lit.');
			xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
			print_plot([workdir ldata(i).name '_lit_comp_shift']);
		end
	end
end