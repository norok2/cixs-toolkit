function plot_linsolve_results(data, expt, workdir)
	% plot linsolve results

	% plot calculated S(q,w) e S(q0,qh,w) point by point before averaging
	figure; clf;
	hold on
	len_dE = length(data.val.dE);
	for l=1:data.val.linsolve.point.set_len
		for k=1:expt.sample.n_terms
			switch (mod(k, 4))
				case {1}
					symbol = 'x';
				case {2}
					symbol = 'o';
				case {3}
					symbol = '*';
				case {0}
					symbol = '>';
			end
			if (sum(data.val.linsolve.weight(:,l)) >= len_dE)
				plot(data.val.dE, data.val.linsolve.dsf_raw(:,k,l), [symbol id2color(k+2, 'red')]);
			end
		end
	end
	hold off
	legend('dDSF', 'nDSF');
	xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (arb. units)');
	print_plot([workdir expt.name '_linsolve_raw_dsf_' sprintf('%d', expt.sample.n_terms) '_' data.val.linsolve.tuples '']);

	% plot calculated S(q,w) e S(q0,qh,w)
	figure; clf;
	hold on
		errorbar(data.val.dE, data.val.linsolve.result.ddsf, data.val.linsolve.result.D_ddsf, ['>' id2color(1+2, 'red')]);
		errorbar(data.val.dE, data.val.linsolve.result.ndsf, data.val.linsolve.result.D_ndsf, ['>' id2color(2+2, 'red')]);
	hold off
	legend('dDSF', 'nDSF');
	xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
	print_plot([workdir expt.name '_linsolve_dsf']);
end
