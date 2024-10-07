function unit_cell = standard_unit_cell(unit_cell)
	for k=1:length(unit_cell(1))
		for i=1:3
			if (unit_cell(k,i) < 0)
				unit_cell(k,i) = unit_cell(k,i) + 1;
			end
		end
	end
end