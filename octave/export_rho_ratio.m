function export_rho_ratio(data, experiment, workdir)
	% export rho ratio data calculated with different methods

	% export result value (should be equal to one below)
	if (~isempty(data.val.result.rho_ratio))
		rho_ratio(1,:) = [data.val.result.rho_ratio, data.val.result.D_rho_ratio];
	else
		rho_ratio(1,:) = [0, 0];
	end

	% export linsolve value
	if (~isempty(data.val.linsolve.result.rho_ratio))
		rho_ratio(2,:) = [data.val.linsolve.result.rho_ratio, data.val.linsolve.result.D_rho_ratio];
	else
		rho_ratio(2,:) = [0, 0];
	end

	% export SVD value
	if (~isempty(data.val.svd.result.rho_ratio))
		rho_ratio(3,:) = [data.val.svd.result.rho_ratio, data.val.svd.result.D_rho_ratio];
	else
		rho_ratio(3,:) = [0, 0];
	end

	% export literature value
%  	if (~isempty(data.val.svd.rho_ratio))
%  		rho_ratio(4,:) = [data.val.svd.rho_ratio, data.val.svd.D_rho_ratio];
%  	else
%  		rho_ratio(4,:) = [0, 0];
%  	end

	csvwrite([workdir experiment.name '_rho_ratio.csv'], rho_ratio);
end
