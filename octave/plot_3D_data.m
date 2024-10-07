function plot_3D_data(dth, dE, val, expt, name, workdir)
	% plot all figures for given 3D data
	
	param = expt.param;
	len_dth = length(dth);
	len_dE = length(dE);

	val = val * 1.1e1;
	dth = dth * 1e3;
	param.plot.dth.min = param.plot.dth.min * 1e3;
	param.plot.dth.max = param.plot.dth.max * 1e3;

	% plot full data (contour)
	figure; clf;
	contour(dth, dE, val);
	xlabel('Angular Departure (mrad)'); ylabel('Energy (eV)');
%  	xlim([min(dth) max(dth)]); ylim([min(dE), max(dE)]);
	colorbar();
	print_plot([workdir expt.name '_' name '_contour']);

	% plot central part (surface and contour)
	range_dth = ((dth > param.plot.dth.min) & (dth < param.plot.dth.max));
	range_dE = ((dE > param.plot.dE.min) & (dE < param.plot.dE.max));
	[mesh_range_dth, mesh_range_dE] = meshgrid(dth(range_dth), dE(range_dE));
	figure; clf; %set(gca, 'fontsize', 12);
	surf(mesh_range_dth, mesh_range_dE, val(range_dE,range_dth)); view(210.0, 30.0); set(gca,'XDir','reverse');
	xlabel('\Delta\theta (mrad)'); %set(get(gca, 'XLabel'), 'FontSize', 20.0);
	ylabel('E (eV)'); %set(get(gca, 'YLabel'), 'FontSize', 20.0);
  	zlabel('I (arb. units)'); set(get(gca, 'ZLabel'), 'Rotation', 90.0); %set(get(gca, 'ZLabel'), 'FontSize', 20.0);
	colorbar();
	print_plot([workdir expt.name '_' name '_surfc']);

	% plot energy sequence
	id_E = fix(1:(len_dE / param.plot.n_seq):len_dE);
	shift_E = (max(max(val(id_E,:))) - min(min(val(id_E,:)))) * (16 / 16);
	figure; clf;
	hold on
	for i=1:length(id_E)
		plot_legend{i} = ['eV = ' sprintf('%f', dE(id_E(i)))];
%  		plot(dth .* 1e3, ones(size(dth)) * shift_E * (i - 1), ['--', 'y']);
		plot(dth(range_dth), val(id_E(i),range_dth) + shift_E * (i - 1), ['-', id2color(i, 'red')]);
	end
	hold off
%  	legend(plot_legend, 'location', 'outside');
	xlabel('Angular Departure (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_' name '_dE_seq']);

	% plot angular sequence
	id_th = fix(1:(len_dth / param.plot.n_seq):len_dth);
	shift_th = (max(max(val(:,id_th))) - min(min(val(:,id_th)))) * (16 / 16);
	figure; clf;
	hold on
	for i=1:length(id_th)
		plot_legend{i} = ['mrad = ' sprintf('%f', dth(id_th(i)) * 1e3)];
%  		plot(dE, ones(size(dE)) * shift_th * (i - 1), ['--', 'y']);
		plot(dE(range_dE), val(range_dE,id_th(i)) + shift_th * (i - 1), ['-', id2color(i, 'red')]);
	end
	hold off
%  	legend(plot_legend, 'location', 'outside');
	xlabel('Energy (eV)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_' name '_dth_seq']);
end
