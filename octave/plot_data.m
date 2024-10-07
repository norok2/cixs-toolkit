function plot_data(data, expt, workdir)
	% plot final input data

	% plot definitive input data
	plot_3D_data(data.val.dth, data.val.dE, data.val.ss, expt, 'input_data', workdir);
end
