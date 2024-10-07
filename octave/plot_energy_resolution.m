function plot_energy_resolution(data, expt, workdir)
	% plot beam profile

	data.dth = angular_departure(data.th, bragg_angle(expt.sample, expt.E_in, expt.sample.hkl)); %DEBUG
	data.exact.dth = data.dth;
	data.approx.dth = data.dth + data.shift;

	plot_dumond(data.dth, data.E - expt.E_in, data.ib ./ 100000, expt, 'beam', workdir);
	plot_dumond(data.dth, data.E - expt.E_in, data.mm, expt, 'monochr', workdir);
	plot_dumond(data.dth, data.E - expt.E_in, data.ibmm ./ 100000, expt, 'beam_monochr', workdir);
	plot_dumond(data.dth, data.E - expt.E_in, data.rr, expt, 'rr', workdir);
	plot_dumond(data.dth, data.E - expt.E_in, data.aa, expt, 'aa', workdir);
	plot_dumond(data.dth, data.E - expt.E_in, data.bb, expt, 'bb', workdir);
	plot_dumond(data.dth, data.E - expt.E_in, data.cc, expt, 'cc', workdir);
	plot_dumond(data.dth, data.E - expt.E_in, data.dd, expt, 'dd', workdir);


	len_raw = length(data.dth);
	range_raw = fix((len_raw*4/16):(len_raw*12/16));
	energy_raw = length(data.E) / 2;
	len_exact = length(data.exact.dth);
	range_exact = fix((len_exact*4/16):(len_exact*12/16));
	len_approx = length(data.approx.dth);
	range_approx = fix((len_approx*4/16):(len_approx*12/16));
	
	% aa coefficient comparison
	figure; clf;
	hold on
		plot(data.exact.dth(range_exact) .* 1e3, data.exact.aa(range_exact), ['-' 'k']);
		plot(data.approx.dth(range_approx) .* 1e3, data.approx.aa(range_approx), ['--' id2color(1)]);
	hold off
	legend('A''', 'A''''', 'location', 'northeast');
	xlabel('Angle (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_comp_conv_aa']);

	% bb coefficient comparison
	figure; clf;
	hold on
		plot(data.exact.dth(range_exact) .* 1e3, data.exact.bb(range_exact), ['-' 'k']);
		plot(data.approx.dth(range_approx) .* 1e3, data.approx.bb(range_approx), ['--' id2color(3)]);
	hold off
	legend('B''', 'B''''', 'location', 'northeast');
	xlabel('Angle (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_comp_conv_bb']);

	% cc coefficient comparison
	figure; clf;
	hold on
		plot(data.exact.dth(range_exact) .* 1e3, data.exact.cc(range_exact), ['-' 'k']);
		plot(data.approx.dth(range_approx) .* 1e3, data.approx.cc(range_approx), ['--' id2color(2)]);
	hold off
	legend('C''', 'C''''', 'location', 'northeast');
	xlabel('Angle (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_comp_conv_cc']);

	% dd coefficient comparison
	figure; clf;
	hold on
		plot(data.exact.dth(range_exact) .* 1e3, data.exact.dd(range_exact), ['-' 'k']);
		plot(data.approx.dth(range_approx) .* 1e3, data.approx.dd(range_approx), ['--' id2color(4)]);
	hold off
	legend('D''', 'D''''', 'location', 'northeast');
	xlabel('Angle (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_comp_conv_dd']);

	% rr coefficient comparison
	figure; clf;
	hold on
		plot(data.exact.dth(range_exact) .* 1e3, data.exact.rr(range_exact), ['-' 'k']);
		plot(data.approx.dth(range_approx) .* 1e3, data.approx.rr(range_approx), ['--' id2color(5)]);
	hold off
	legend('R''', 'R''''', 'location', 'northeast');
	xlabel('Angle (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_comp_conv_rr']);

	% raw wavefield coefficients comparison
	figure; clf;
	hold on
		plot(data.dth(range_raw) .* 1e3, data.aa(range_raw,energy_raw), ['-' id2color(1)]);
		plot(data.dth(range_raw) .* 1e3, data.bb(range_raw,energy_raw), ['-' id2color(3)]);
		plot(data.dth(range_raw) .* 1e3, data.cc(range_raw,energy_raw), ['-' id2color(2)]);
		plot(data.dth(range_raw) .* 1e3, data.dd(range_raw,energy_raw), ['-' id2color(4)]);
		plot(data.dth(range_raw) .* 1e3, data.rr(range_raw,energy_raw), ['-' id2color(5)]);
	hold off
	legend('A', 'B', 'C', 'D', 'R', 'location', 'northeast');
	xlabel('Angle (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_comp_raw']);

	% exact energy resolution comparison
	figure; clf;
	hold on
		plot(data.exact.dth(range_exact) .* 1e3, data.exact.aa(range_exact), ['-' id2color(1)]);
		plot(data.exact.dth(range_exact) .* 1e3, data.exact.bb(range_exact), ['-' id2color(3)]);
		plot(data.exact.dth(range_exact) .* 1e3, data.exact.cc(range_exact), ['-' id2color(2)]);
		plot(data.exact.dth(range_exact) .* 1e3, data.exact.dd(range_exact), ['-' id2color(4)]);
		plot(data.exact.dth(range_exact) .* 1e3, data.exact.rr(range_exact), ['-' id2color(5)]);
	hold off
	legend('A''', 'B''', 'C''', 'D''', 'R''', 'location', 'northeast');
	xlabel('Angle (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_comp_exact']);

	% approximate energy resolution comparison
	figure; clf;
	hold on
		plot(data.approx.dth(range_approx) .* 1e3, data.approx.aa(range_approx), ['-' id2color(1)]);
		plot(data.approx.dth(range_approx) .* 1e3, data.approx.bb(range_approx), ['-' id2color(3)]);
		plot(data.approx.dth(range_approx) .* 1e3, data.approx.cc(range_approx), ['-' id2color(2)]);
		plot(data.approx.dth(range_approx) .* 1e3, data.approx.dd(range_approx), ['-' id2color(4)]);
		plot(data.approx.dth(range_approx) .* 1e3, data.approx.rr(range_approx), ['-' id2color(5)]);
	hold off
	legend('A''''', 'B''''', 'C''''', 'D''''', 'R''''', 'location', 'northeast');
	xlabel('Angle (rad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_comp_approx']);
end
