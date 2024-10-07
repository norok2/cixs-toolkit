function print_plot(filename)
	
	targets = {'eps', 'png', 'svg'};
	
	font = sprintf('-F%s:%d', '', 12);
	
	for i=1:numel(targets)
		ext = char(targets(i));
		plotfile = [filename '.' ext];
		if (strcmp(ext, 'svg'))
%  			base = 200;
%  			size = sprintf('-S%d,%d', 4 * base, 3 * base);
%  			font = sprintf('-F%s:%d', '', 13);
			args = {'-solid', '-color', '-r600', '-dsvg', font};
		elseif (strcmp(ext, 'png'))
			args = {'-solid', '-color', '-r150', '-dpng', font};
		elseif (strcmp(ext, 'eps'))
			args = {'-solid', '-color', '-r600', '-depsc2', font};
		end
		args = [args {plotfile}];
		print(args{:});
	end
end