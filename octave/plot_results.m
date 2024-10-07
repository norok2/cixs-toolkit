function plot_results(data, expt, workdir)
	% plot final results

	% plot comparison between linsolve and svd method for S(q,w) e S(q0,qh,w)
	figure; clf;
	hold on
		errorbar(data.val.dE, data.val.linsolve.result.ddsf, data.val.linsolve.result.D_ddsf, ['>' id2color(3, 'red')]);
		errorbar(data.val.dE, data.val.linsolve.result.ndsf, data.val.linsolve.result.D_ndsf, ['>' id2color(4, 'red')]);
		errorbar(data.val.dE, data.val.svd.result.ddsf, data.val.svd.result.D_ddsf, ['>' id2color(1, 'red')]);
		errorbar(data.val.dE, data.val.svd.result.ndsf, data.val.svd.result.D_ndsf, ['>' id2color(2, 'red')]);
	hold off
	legend('dDSF (IAS)', 'nDSF (IAS)', 'dDSF (SVD)', 'nDSF(SVD)');
	xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
	print_plot([workdir expt.name '_dsf_comparison']);

	% plot calculated S(q,w) e S(q0,qh,w)
	figure; clf;
	hold on
		errorbar(data.val.dE, data.val.result.ddsf, data.val.result.D_ddsf, ['>' id2color(1)]);
		errorbar(data.val.dE, data.val.result.ndsf, data.val.result.D_ndsf, ['>' id2color(2)]);
	hold off
	legend('dDSF', 'nDSF');
	xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
	print_plot([workdir expt.name '_norm_dsf']);

	% plot E*S(q,E) e E*S(q0,qh,E)
	figure; clf;
	hold on
		errorbar(data.val.dE, data.val.dE .* data.val.result.ddsf, data.val.dE .* data.val.result.D_ddsf, ['>' id2color(1)]);
		errorbar(data.val.dE, data.val.dE .* data.val.result.ndsf, data.val.dE .* data.val.result.D_ndsf, ['>' id2color(2)]);
	hold off
	legend('dDSF', 'nDSF');
	xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor S times Energy (pure number)');
	print_plot([workdir expt.name '_norm_E-x-dsf']);

%  	figure; clf;
%  	hold on
%  		loglog(data.val.dE, data.val.dE .* data.val.result.ddsf, ['-*' id2color(1)]);
%  		loglog(data.val.dE, data.val.dE .* data.val.result.ndsf, ['-*' id2color(2)]);
%  	hold off
%	legend('E dDSF', 'E nDSF');
%  	xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor S times Energy (pure number)');

%  	% plot calculated Im[epsilon^-1](q,w) e Im[epsilon^-1](q0,qh,w)
%  	figure; clf;
%  	hold on
%  		errorbar(data.val.dE, data.val.result.deps, data.val.result.D_deps, ['>' id2color(1)]);
%  		errorbar(data.val.dE, data.val.result.neps, data.val.result.D_neps, ['>' id2color(2)]);
%  	hold off
%  	xlabel('Energy (eV)'); ylabel('Im[-1/\epsilon] (pure number)');
%  	print_plot([workdir expt.name '_norm_epsilon']);

%  	% plot comparison with literature
%  	if (~isempty(data.lit))
%  		figure; clf;
%  		hold on
%  			errorbar(data.val.dE, data.val.result.ddsf, data.val.result.D_ddsf, ['>' id2color(1)]);
%  			errorbar(data.val.dE, data.val.result.ndsf, data.val.result.D_ndsf, ['>' id2color(2)]);
%  			plot(data.lit.data(:,1) + data.lit.shift, data.lit.data(:,2) .* data.lit.norm, ['--' id2color(1)]);
%  			plot(data.lit.data(:,1) + data.lit.shift, data.lit.data(:,3) .* data.lit.norm, ['--' id2color(2)]);
%  		hold off
%  		legend('dDSF (lit.)', 'nDSF (lit.)', 'dDSF', 'nDSF');
%  		xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
%  		print_plot([workdir expt.name '_lit_comp']);
%  	end
end
