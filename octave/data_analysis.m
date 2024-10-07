clear; tic;

configure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% execute
%  execute.standing_wave_calculation.all = false;
execute.energy_resolution.all = false;
execute.data_analysis.all = true;
debug = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% standing waves calculation (theory)
%  if (execute.standing_wave_calculation.all == true)
%  	calculate_standing_wave(expt);
%  	% produce figures for standing wave quantities in dynamical theory of X-ray diffraction
%  	plot_standing_wave(expt);
%  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% precise convolution
if (execute.energy_resolution.all == true)
	execute.energy_resolution.calculate = false;
	execute.energy_resolution.plot = false;

	experiment = expt(1);
	expt_data_file = [export_dir experiment.name '_energy_resolution.mat'];

	if ((execute.energy_resolution.calculate == true) || isempty(dir(expt_data_file)))
		data = calculate_energy_resolution(expt);
		save(expt_data_file, 'data');
	else
		load(expt_data_file, 'data');
	end

	if ((execute.energy_resolution.plot == true) && ~isempty(data))
		plot_energy_resolution(data, experiment, plot_dir);
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data analysis technique

if (execute.data_analysis.all == true)
	execute.data_analysis.extract = true;
	execute.data_analysis.preprocess = true;
	execute.data_analysis.dsf_analysis = true;
	execute.data_analysis.plot = true;
	execute.data_analysis.export = true;
	execute.data_analysis.all_expt = true;

%  	expt_dbg = expt; %DEBUG%
	if (execute.data_analysis.all_expt == true)
		first_expt = 1;
		last_expt = expt_num;
	else
		n = 1; % change this to analyze single run
		first_expt = n;
		last_expt = n;
	end


	% calculate all results for each experiment
	for i=first_expt:last_expt
		experiment = expt(i);
		expt_data_file = [export_dir experiment.name '_data.mat'];

		% extract experimental data file
		if ((execute.data_analysis.extract == true) || isempty(dir(expt_data_file)))
			% extract signal (allow multiple experimental runs)
			data.val = extract_data(experiment.signal, experiment.param, experiment.E_in, experiment.beam, experiment.monochr, experiment.sample);
			% extract background
			if (~isempty(experiment.background))
				data.bg = extract_data(experiment.background, experiment.param, experiment.E_in, experiment.beam, experiment.monochr, experiment.sample);
			end

			data.executed = 'extract';
			save(expt_data_file, 'data');
		else
			load(expt_data_file, 'data');
		end

		% prepare input data
		if (execute.data_analysis.preprocess == true || strcmp(data.executed, 'extract'))
			% remove background
			if (~isempty(data.bg))
				data = remove_background(data);
			else
				data.val.ss_nobg = data.val.ss_raw;
			end
			% remove elastic peak
			data = remove_elastic_peak(data);
			% narrow down range of used energy and angle
			data = signal_windowing(data, experiment.param);
			% calculate standing wave parameters
			data = compute_standing_wave_coeff(data, experiment);

			data.executed = 'preprocess';
			save(expt_data_file, 'data');
		else
			load(expt_data_file, 'data');
		end

		% plot input data and checks
		if (execute.data_analysis.plot == true)
			if ((i == first_expt) && (debug == true))
				plot_raw_data(data, experiment, plot_dir);
			end
			plot_check(data, experiment, plot_dir);
			plot_data(data, experiment, plot_dir);
		end

		% calculate diagonal and non diagonal part of dynamic structure factor
		if (execute.data_analysis.dsf_analysis == true || strcmp(data.executed, 'preprocess'))
			% use solution of linear approximation problem selecting input
  			data = compute_dsf_linsolve(data, experiment, 'selected');
			% use singular value decomposition
			data = compute_dsf_svd(data, experiment);

			% choose best extraction method (SVD)
 			data.val.result = data.val.svd.result;
			data.val.result.dE = data.val.dE;

			data.executed = 'dsf_analysis';
			save(expt_data_file, 'data');
		else
			load(expt_data_file, 'data');
		end

		% plot results
		if (execute.data_analysis.plot == true)
%  			plot_linsolve_results(data, experiment, plot_dir);
			plot_svd_results(data, experiment, plot_dir);
			plot_results(data, experiment, plot_dir);
		end

		% export results
		if (execute.data_analysis.export == true)
			export_rho_ratio(data, experiment, export_dir);
		end

		% save result to experiment
		expt(i).result = data.val.result;
%  		expt_dbg(i).data = data; %DEBUG%
	end


	% final calculation involving all experiments
	if (execute.data_analysis.all_expt == true)
		% combine spectra
		data.combined = combine_spectra(expt, groups);
		% filter spectra for high frequency noise
		data.filtered = filter_spectra(data.combined, 3);
		% import literature values
		data.literature = import_literature(data.combined, import_dir);

		% export all results
		if (execute.data_analysis.export == true)
			save([export_dir 'data_results.mat'], 'data');
			save([export_dir 'expt_results.mat'], 'expt');
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% evidence of plasmon band theory
		if (execute.data_analysis.plot == false)
			% plot comparison of results with each other
			plot_compare(expt, groups, plot_dir);
			% plot combined results
			plot_combined(data.combined, plot_dir);
			% plot filtered results
			plot_filtered(data.filtered, plot_dir);
			% plot literature values
			plot_literature(data.combined, data.literature, plot_dir);
			% plot plasmons
			pdata = plot_plasmons(data.combined, plot_dir); % to be modified
		end
	end
end

time.sec = toc;
time.min = time.sec / 60.0;
time.hor = time.min / 60.0;
time.day = time.hor / 24.0;
time
