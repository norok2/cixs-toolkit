clear; tic;

configure;

load([export_dir 'data_results.mat'], 'data');
load([export_dir 'expt_results.mat'], 'expt');

function R = calc_q(q)
	q_c = 1.0325; % q_critical in SI units
	q_au = 1.6250; % atomic units conversion factor
	qq = sqrt(q * q');
	R = [qq / q_c, qq / q_au];
end


for i=1:2:8  
	calc_q(expt(i).sample.q)
end

% get extra literature data
data.extralit(1).name = 'ndsf-0.3[3,1,1]';
data.extralit(1).q0 = 0.3*[3,1,1];
data.extralit(1).hkl = [1,1,1];
data.extralit(1).qh = data.extralit(1).q0 + data.extralit(1).hkl;
data.extralit(1).norm = 0.032990;
data.extralit(1).shift = 1.9249;
data.extralit(1).file = [import_dir 'literature_' data.extralit(1).name '.csv'];
data.extralit(1).data = dlmread(data.extralit(1).file, ',', 2, 0);

first_group = 1
last_group = 3  % 'groups.n' for all
num_expt = last_group - first_group + 1

disp(groups.n)

figure; clf;
plot_legend = {};
hold on
	for i = first_group:last_group
		errorbar(data.filtered(i).dE, data.filtered(i).ddsf, data.filtered(i).D_ddsf, ['>' id2color(i) ';' data.filtered(i).name ';']);
		plot_legend{i} = ['D-' data.filtered(i).name];
	end
	%kaprolat
	plot(data.literature(1).data(:,1) + data.literature(1).shift, data.literature(1).data(:,2) .* data.literature(1).norm, ['-o' id2color(1 + num_expt)]);
	plot_legend{1 + num_expt} = ['D-kaprolat-' data.literature(1).name];
hold off
legend(plot_legend);
xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
print_plot([plot_dir 'ddsf_trend']);

figure; clf;
plot_legend = {};
hold on
	for i = first_group:last_group
		errorbar(data.filtered(i).dE, data.filtered(i).ndsf, data.filtered(i).D_ndsf, ['>' id2color(i) ';' data.filtered(i).name ';']);
		plot_legend{i} = ['N-' data.filtered(i).name];
	end
	%kaprolat
	plot(data.literature(1).data(:,1) + data.literature(1).shift, data.literature(1).data(:,3) .* data.literature(1).norm, ['-o' id2color(1 + num_expt)]);
	plot_legend{1 + num_expt} = ['N-kaprolat-' data.literature(1).name];
	%ehrnsperger
	plot(data.extralit(1).data(:,1) + data.extralit(1).shift, data.extralit(1).data(:,2) .* data.extralit(1).norm, ['-o' id2color(2 + num_expt)]);
	plot_legend{2 + num_expt} = ['N-ehrnsperger-' data.extralit(1).name];
hold off
legend(plot_legend);
xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
print_plot([plot_dir 'ndsf_trend']);

figure; clf;
plot_legend = {};
hold on
	for i = first_group:last_group
		errorbar(data.filtered(i).dE, data.filtered(i).ddsf, data.filtered(i).D_ddsf, ['>' id2color(i) ';' data.filtered(i).name ';']);
		errorbar(data.filtered(i).dE, data.filtered(i).ndsf, data.filtered(i).D_ndsf, ['>' id2color(i) ';' data.filtered(i).name ';']);
		plot_legend{2*i-1} = ['D-' data.filtered(i).name];
		plot_legend{2*i-0} = ['N-' data.filtered(i).name];
	end
	%kaprolat
	plot(data.literature(1).data(:,1) + data.literature(1).shift, data.literature(1).data(:,2) .* data.literature(1).norm, ['-o' id2color(1 + num_expt)]);
	plot(data.literature(1).data(:,1) + data.literature(1).shift, data.literature(1).data(:,3) .* data.literature(1).norm, ['-o' id2color(1 + num_expt)]);
	plot_legend{1 + 2*num_expt} = ['D-kaprolat-' data.literature(1).name];
	plot_legend{2 + 2*num_expt} = ['N-kaprolat-' data.literature(1).name];
	%ehrnsperger
	plot(data.extralit(1).data(:,1) + data.extralit(1).shift, data.extralit(1).data(:,2) .* data.extralit(1).norm, ['-o' id2color(2 + num_expt)]);
	plot_legend{3 + 2*num_expt} = ['N-ehrnsperger-' data.extralit(1).name];
hold off
legend(plot_legend);
xlabel('Energy (eV)'); ylabel('Dynamic Structure Factor (eV^{-1})');
print_plot([plot_dir 'dsf_trend']);

time.sec = toc;
time.min = time.sec / 60.0;
time.hor = time.min / 60.0;
time.day = time.hor / 24.0;
time

