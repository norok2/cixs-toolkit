function C = iconv2(A, B, shape)
	% calculate convolution of matrix A with matrix B and integrate result in integration dimension(s) dimension

	if (nargin < 3)
		shape = 'full';
	end
	
	[rA, cA] = size(A);
	[rB, cB] = size(B);
	if (rA ~= rB)
		error('Error: matrix dimension must agree');
	end

	switch shape
		case 'full'
			CC = zeros(rA, cA + cB - 1);
			C = zeros(1, cA + cB - 1);
		case 'same'
			CC = zeros(rA, cA);
			C = zeros(1, cA);
		otherwise
			error('Error: invalid shape of resulting matrix');
	end

	for i=1:rA
		CC(i,:) = conv2(A(i,:), centroidize(B(i,:)), shape);
	end

	for j=1:length(C)
		C(j) = sum(CC(:,j));
	end
end