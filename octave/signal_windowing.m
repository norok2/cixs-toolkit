function data = signal_windowing(data, param)
	% narrow down (windowing) energy and angle values effectively used

	% compute angular range
	range_dth = ((data.val.dth_raw > param.dth.min) & (data.val.dth_raw < param.dth.max));
	data.val.dth = data.val.dth_raw(range_dth);

	% compute energy range
	param.dE.min = max(param.dE.min, 2*data.val.elastic_width + data.val.elastic_E);
	range_dE = ((data.val.dE_raw > param.dE.min) & (data.val.dE_raw < param.dE.max));
	data.val.dE = data.val.dE_raw(range_dE) - data.val.elastic_E;

	% generate mesh
	[data.val.mesh_dth, data.val.mesh_dE] = meshgrid(data.val.dth, data.val.dE);
	% windowed input signal
	data.val.ss = data.val.ss_noep(range_dE, range_dth);
	data.val.D_ss = data.val.D_ss_noep(range_dE, range_dth);
	% windowed reflectivity
	data.val.rr = data.val.rr_raw(range_dE, range_dth);
end
