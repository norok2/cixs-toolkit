function pdata = plot_plasmons(cdata, workdir)
	% plot plasmons

	% plot plasmons
	i = 1; % only case where it is sensible
	pdata.name = cdata(i).name;
	% calculate plasmon band
	pdata.dE = cdata(i).dE;
	pdata.u = (cdata(i).ddsf - cdata(i).ndsf) ./ 2.0;
	pdata.D_u = (cdata(i).D_ddsf + cdata(i).D_ndsf) ./ 2.0;
	pdata.l = (cdata(i).ddsf + cdata(i).ndsf) ./ 2.0;
	pdata.D_l = (cdata(i).D_ddsf + cdata(i).D_ndsf) ./ 2.0;
	% plot plasmon band
	figure; clf;
	hold on
		errorbar(cdata(i).dE, pdata.l, pdata.D_l, ['>' id2color(2)]);
		errorbar(cdata(i).dE, pdata.u, pdata.D_u, ['>' id2color(1)]);
	hold off
	legend('P_l', 'P_u');
	xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
	print_plot([workdir pdata.name '_plasmons']);
end
