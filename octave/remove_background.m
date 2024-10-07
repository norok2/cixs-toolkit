function data = remove_background(data)
	% subtract background from val

	% interpolate background data at val (dth,dE) coordinates
	data.val.ss_bg = interp2(data.bg.mesh_dth_raw, data.bg.mesh_dE_raw, data.bg.ss_raw, data.val.mesh_dth_raw, data.val.mesh_dE_raw);
	data.val.D_ss_bg = interp2(data.bg.mesh_dth_raw, data.bg.mesh_dE_raw, data.bg.D_ss_raw, data.val.mesh_dth_raw, data.val.mesh_dE_raw);

	% remove NaN at borders of background
	if (sum(isnan(data.val.ss_bg(1,:)) > 0)) % lowest energy
		data.val.ss_bg(1,:) = interp1(data.bg.dth_raw, data.bg.ss_raw(1,:), data.val.dth_raw);
	end
	if (sum(isnan(data.val.ss_bg(end,:)) > 0)) % highest energy
		data.val.ss_bg(end,:) = interp1(data.bg.dth_raw, data.bg.ss_raw(end,:), data.val.dth_raw);
	end
	if (sum(isnan(data.val.ss_bg(:,1)) > 0)) % lowest angle
		data.val.ss_bg(:,1) = interp1(data.bg.dE_raw, data.bg.ss_raw(:,1), data.val.dE_raw);
	end
	if (sum(isnan(data.val.ss_bg(:,end)) > 0)) % highest angle
		data.val.ss_bg(:,end) = interp1(data.bg.dE_raw, data.bg.ss_raw(:,end), data.val.dE_raw);
	end
	% remove NaN at borders of background error
	if (sum(isnan(data.val.D_ss_bg(1,:)) > 0)) % lowest energy
		data.val.D_ss_bg(1,:) = interp1(data.bg.dth_raw, data.bg.D_ss_raw(1,:), data.val.dth_raw);
	end
	if (sum(isnan(data.val.D_ss_bg(end,:)) > 0)) % highest energy
		data.val.D_ss_bg(end,:) = interp1(data.bg.dth_raw, data.bg.D_ss_raw(end,:), data.val.dth_raw);
	end
	if (sum(isnan(data.val.D_ss_bg(:,1)) > 0)) % lowest angle
		data.val.D_ss_bg(:,1) = interp1(data.bg.dE_raw, data.bg.D_ss_raw(:,1), data.val.dE_raw);
	end
	if (sum(isnan(data.val.D_ss_bg(:,end)) > 0)) % highest angle
		data.val.D_ss_bg(:,end) = interp1(data.bg.dE_raw, data.bg.D_ss_raw(:,end), data.val.dE_raw);
	end

	
	for j=1:length(data.val.dth_raw)
		% remove other NaN background values
		replacement = nanmean(data.val.ss_bg(:,j));
		nan_index = isnan(data.val.ss_bg(:,j));
		data.val.ss_bg(nan_index,j) = replacement;
		% remove other NaN background error values
		replacement = nanmean(data.val.D_ss_bg(:,j));
		nan_index = isnan(data.val.D_ss_bg(:,j));
		data.val.D_ss_bg(nan_index,j) = replacement;
	end

	% subtract background and calculate error
	data.val.ss_nobg = data.val.ss_raw - data.val.ss_bg;
	data.val.D_ss_nobg = sqrt(data.val.D_ss_raw.^2.0 + data.val.D_ss_bg.^2.0);
end
