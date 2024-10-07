function [background, D_background] = background_interp(th_exp, E_exp, th_B_exp, data_exp, crystal, hkl)
	
	% create numerical function from data points
	data_3si_linear = TriScatteredInterp(data_exp(:,1), data_exp(:,2), data_exp(:,3), 'linear');
	D_data_3si_linear = TriScatteredInterp(data_exp(:,1), data_exp(:,2), data_exp(:,4), 'linear');
	data_3si_nearest = TriScatteredInterp(data_exp(:,1), data_exp(:,2), data_exp(:,3), 'nearest');
	D_data_3si_nearest = TriScatteredInterp(data_exp(:,1), data_exp(:,2), data_exp(:,4), 'nearest');
	
	background = zeros(size(th_exp)); D_background = zeros(size(th_exp));
	for i=1:length(th_exp)
		background(i) = data_3si_linear(th_exp(i), E_exp);
		D_background(i) = D_data_3si_linear(th_exp(i), E_exp);
		if (isnan(background(i)) || isnan(D_background(i)))
			background(i) = data_3si_nearest(th_exp(i), E_exp);
			D_background(i) = D_data_3si_nearest(th_exp(i), E_exp);
		end
	end
	
% 	% parameters for manual handling
% 	N = length(data_exp(:,1));
% 	N_th = 0;
% 	while ((N_th < N) && (data_exp(N_th+1,1) == data_exp(1,1)))
% 		N_th = N_th + 1;		
% 	end
% 	N_E = N / N_th;
% 	E = zeros(1, N_E);
% 	for i=1:N_E
% 		E(i) = mean(data_exp((i-1)*N_th+1:i*N_th,1));
% 	end
% 	th = zeros(1, N_th);
% 	for i=1:N_th
% 		th(i) = mean(data_exp(i:N_th:N,2));
% 	end
% 		
% 	% perform manual energy interpolation
% 	for i=1:length(th_exp)
% 		E_var = i:N_th:N;
% 		background(i) = interp1(data_exp(E_var,1), data_exp(E_var,3), E_exp, 'linear');
% 		D_background(i) = interp1(data_exp(E_var,1), data_exp(E_var,4), E_exp, 'linear');
% 		if (isnan(background(i)) || isnan(D_background(i)))
% 			background(i) = interp1(data_exp(E_var,1), data_exp(E_var,3), E_exp, 'nearest');
% 			D_background(i) = interp1(data_exp(E_var,1), data_exp(E_var,4), E_exp, 'nearest');
% 		end
% 	end

	% perform manual th interpolation averaged over energy
% 	th_B = bragg_angle(crystal, E, hkl);
% 	energy_bias = zeros(size(E));
% 	for i=1:N_E
% 		energy_bias(i) = mean(data_exp((i-1)*N_th+1:i*N_th,3));
% 	end
% 	
% 	for i=1:N_th
% 		th_mean(i) = mean(data_exp(i:N_th:N,3));
% 		th_D_mean(i) = mean(data_exp(i:N_th:N,4));
% 		th_stdev(i) = std(data_exp(i:N_th:N,3));
% 	end

% 	figure; clf;
% 	hold on
% 		p = polyfit(E, energy_bias, 2);
% 		x = E(1):0.0001:E(end);
% 		plot(E, energy_bias);
% 		plot(x, polyval(p, x));
% 	hold off
% 	
% 	figure; clf;
% 	hold on
% 		plot(th, th_mean, '-r');
% 		plot(th, th_D_mean, '-b');
% 		plot(th, th_stdev, '-g');	
% 	hold off
end