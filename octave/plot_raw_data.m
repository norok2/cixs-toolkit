function plot_raw_data(data, expt, workdir)
	% plot raw data

	% plot raw data
	plot_3D_data(data.val.dth_raw, data.val.dE_raw, data.val.ss_raw, expt, 'raw_data', workdir);
	plot_3D_data(data.val.dth_raw, data.val.dE_raw, data.val.rr_raw, expt, 'reflectivity', workdir);
	plot_3D_data(data.bg.dth_raw, data.bg.dE_raw, data.bg.ss_raw, expt, 'background', workdir);
	plot_3D_data(data.val.dth_raw, data.val.dE_raw, data.val.ss_nobg, expt, 'raw_data_nobg', workdir);
	plot_3D_data(data.val.dth_raw, data.val.dE_raw, data.val.ss_noep, expt, 'raw_data_noep', workdir);
	plot_3D_data(data.val.dth_raw, data.val.dE_raw, data.val.ep, expt, 'elastic_peak', workdir);

	% plot standing wave coefficients
	plot_3D_data(data.val.dth, data.val.dE, data.val.sw.aa, expt, 'sw_aa', workdir);
	plot_3D_data(data.val.dth, data.val.dE, data.val.sw.bb, expt, 'sw_bb', workdir);
	plot_3D_data(data.val.dth, data.val.dE, data.val.sw.cc, expt, 'sw_cc', workdir);
	plot_3D_data(data.val.dth, data.val.dE, data.val.sw.dd, expt, 'sw_dd', workdir);
	plot_3D_data(data.val.dth, data.val.dE, data.val.sw.rr, expt, 'sw_rr', workdir);
end
