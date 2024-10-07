function plot_svd_results(data, expt, workdir)
	% plot Singular Value Decomposition results

	% plot Singular Value Decomposition result to infer number of linearly independent components
	figure; clf;
	hold on;
	for k=1:length(data.val.dth)
		plot(data.val.dE, data.val.svd.US(:,k), ['-' id2color(k)]);
	end
	hold off
	xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (arb. units)');
	print_plot([workdir expt.name '_svd_raw']);
	
	% plot Singular Value Decomposition result to infer number of linearly independent components
	figure; clf;
	hold on;
	for k=1:expt.sample.n_terms
		plot(data.val.dE, data.val.svd.dsf(:,k), ['-' id2color(k)]);
	end
	hold off
	xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (arb. units)');
	print_plot([workdir expt.name '_svd_reduced']);

	% plot calculated S(q0,w) e S(q0,qh,w)
	figure; clf; xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
	hold on
		errorbar(data.val.dE, data.val.svd.result.ddsf, data.val.svd.result.D_ddsf, ['>' id2color(1)]);
		errorbar(data.val.dE, data.val.svd.result.ndsf, data.val.svd.result.D_ndsf, ['>' id2color(2)]);
	hold off
	legend('dDSF', 'nDSF');
	xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (arb. units)');
	print_plot([workdir expt.name '_svd_dsf']);
end
