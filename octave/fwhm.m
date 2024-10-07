function [width, middle_point] = fwhm(x, y, method, range_index_size)
	% calculate Full-Width Half Maximum of given numeric function using specified interpolation method and index range
	
	if (size(x) ~= size(y))
		y = y';
	end
	if (size(x) ~= size(y))
		error('Vectors must be of equal size');
	end
	% define default interpolation method
	if (nargin < 3)
		method = 'linear';
	end
	% define default index range size over which interpolate
	if (nargin < 4) || (range_index_size < 2)
		range_index_size = 2;
	end
	
	% sort function in ascending order
	[x, sort_index] = sort(x);
	y = y(sort_index);
	% find half-maximum
	half_max = max(y) ./ 2;
	% find the support of the peak
	peak = find(y > (half_max));
	% if there is any peak
	if (length(peak) > 0)
		% determine index ranges over which interpolate
		if (peak(1) > floor(range_index_size/2))
			range_up = floor(peak(1) - (range_index_size - 1)/2):floor(peak(1) + (range_index_size - 1)/2);
		else
			range_up = 1:range_index_size;
		end
		if (peak(end) < length(y) - ceil(range_index_size/2)) 
			range_down = ceil(peak(end) - (range_index_size - 1)/2):ceil(peak(end) + (range_index_size - 1)/2);
		else
			range_down = (length(y) - (range_index_size - 1)):length(y);
		end
		% interpolate for most accurate value
		step_up = interp1(y(range_up), x(range_up), half_max, method);
		step_down = interp1(y(range_down), x(range_down), half_max, method);
		% calculate width
		width = abs(step_up - step_down);
		middle_point = mean([step_up step_down]);
	else
		width = abs(x(end) - x(1));
		middle_point = mean([x(end) x(1)]);
	end
end
