clear; close all;

x = -5:0.1:5;
n = 50;
Y = zeros(length(x), n); W = zeros(n, n);
Y(:,1) = f_gauss(x, -1, 1);
W(:,1) = f_gauss((1:n) - (n/2.0), 0, (1.0/8.0)*n);
Y(:,2) = f_gauss(x, 1, 0.5);
W(:,2) = f_gauss((1:n) - (n/2.0), -(1.0/4.0)*n, (1.0/8.0)*n) - f_gauss((1:n) - (n/2.0), +(1.0/4.0)*n, (1.0/8.0)*n);
r = 2;
for i=r+1:n
	Y(:,i) = 0.03/i .* rand(1, length(x));
	W(:,i) = 0.1 .* rand(1, n);
end

Yr = Y(:,1:r);
Wr = W(:,1:r);

D = Y * W';
Dr = Yr * Wr';

[U, S, V] = svd(D, 'econ');
US = U * S;

USr = US(:,1:r);
Vr = V(:,1:r);

TI = Wr' / Vr';
T = pinv(TI);

%  figure; clf; plot(D); figure; clf; plot(Dr);
%  figure; clf; plot(US * V'); figure; clf; plot(USr * Vr');

figure; clf; plot(Y); figure; clf; plot(Yr);
figure; clf; plot(US); figure; clf; plot(USr);
figure; clf; plot((USr * T));

figure; clf; plot(W); figure; clf; plot(Wr);
figure; clf; plot((TI * Vr')');



