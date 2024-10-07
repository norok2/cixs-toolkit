function plot_check(data, expt, workdir)
	% plot check on wavefield amplitudes and phase shift using reflectivity

	% check shift
	figure; clf;
	plot(data.val.dE_raw * 1e3, data.val.shift, ['-*', 'r']);
	xlabel('Angular Departure (rad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_rr_check_shift']);

	% check reproducibility of reflectivity
	figure; clf;
	hold on
		k = fix(mean([length(data.val.dE) length(data.val.dE_raw)]) / 2.0);
		errorbar(data.val.dth(:) * 1e3, data.val.rr(k,:), data.val.rr(k,:) * (5/100), ['>', 'r']);
		plot(data.val.dth(:) * 1e3, data.val.sw.rr(k,:), ['-', 'b']);
	hold off
	xlabel('Angular Departure (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_rr_check_resolution']);

	% check alignment of experimental reflectivity
	figure; clf;
	hold on
	for k=1:length(data.val.dE_raw)
		kc = k / length(data.val.dE_raw) * (length(color_list()) - 1);
		plot(data.val.dth_raw * 1e3, data.val.rr_raw(k,:), ['-', id2color(kc)]);
	end
	hold off
	xlabel('Angular Departure (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_rr_check_exp_alignment']);

	% check alignment of theoretical reflectivity
	figure; clf;
	hold on
	for k=1:length(data.val.dE)
		kc = k / length(data.val.dE_raw) * (length(color_list()) - 1);
		plot(data.val.dth * 1e3, data.val.sw.rr(k,:), ['-', id2color(kc)]);
	end
	hold off
	xlabel('Angular Departure (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_rr_check_theor_alignment']);
end
