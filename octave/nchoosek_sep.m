function [M] = nchoosek_sep(vector, k, sep)
	M = nchoosek(vector, k);
	
	remove = [];
	for i=1:length(M(:,1))
		removed = false;
		j = 1;
		while ((j <= k - 1) && ~removed)
			jj = j + 1;
			while ((jj <= k) && ~removed)
				if (abs(M(i,j) - M(i,jj)) < sep)
					remove = [remove, i];
					removed = true;
				end
				jj = jj + 1;
			end
			j = j + 1;
		end
	end
	M(remove,:) = [];
end
		