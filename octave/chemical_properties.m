function elem = chemical_properties(id)
	% get chemical properties of elements (name, symbol, Z, A)

	switch id
		case {'Si', 'Silicon', 14}
			elem.name = 'Silicon';
			elem.symbol = 'Si';
			elem.Z = 14;
			elem.A = 28.0855;
		otherwise
			error('unknown element');
	end
end