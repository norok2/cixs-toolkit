function plot_dumond(th, E, ii, expt, name, workdir)
	[r, c] = size(ii);
	mr = round(r / 2.0);
	mc = round(c / 2.0);
	len_th = length(th);
	len_E = length(E);
	len_mesh = expt.param.plot.mesh * 0.3;

	% contour plot of DuMond diagram
	figure; clf;
	contour(th * 1e3, E, ii');
	xlabel('Angle (mrad)'); ylabel('Energy (eV)');
	colorbar();
	print_plot([workdir expt.name '_dumond_' name '_contour']);

	% surface contour plot of DuMond diagram
	range_limit = 0.1;
	range_shift = 0.5;
	range_E = fix(sequenced(1, len_E * range_limit, len_mesh)) + fix(len_E * (1 - range_limit) * (1 - range_shift));
	range_th = fix(sequenced(1, len_th * range_limit, len_mesh)) + fix(len_th * (1 - range_limit) * (1 - range_shift));
	[mesh_range_th, mesh_range_E] = meshgrid(th(range_th) .* 1e3, E(range_E));
	figure; clf;
	surfc(mesh_range_th, mesh_range_E, ii(range_E,range_th));
	view(30.0, 30.0); %set(gca,'XDir','reverse');
	xlabel('\Delta\theta (mrad)'); ylabel('E (eV)');
  	zlabel('I (arb. units)'); set(get(gca, 'ZLabel'), 'Rotation', 90.0);
	colorbar();
	print_plot([workdir expt.name '_dumond_' name '_surfc']);

	% Energy slice of DuMond diagram
	range_th = fix((len_th*5/16):(len_th*11/16));
	figure; clf;
	plot(th(range_th) .* 1e3, ii(range_th,mc), '-k');
	xlabel('Angle (mrad)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_dumond_' name '_E_slicing']);

	% Angular slice of DuMond diagram
	range_E = fix((len_E*5/16):(len_E*11/16));
	figure; clf; 
	plot(E(range_E), ii(mr,range_E), '-k');
	xlabel('Energy (eV)'); ylabel('Intensity (arb. units)');
	print_plot([workdir expt.name '_dumond_' name '_th_slicing']);
end
