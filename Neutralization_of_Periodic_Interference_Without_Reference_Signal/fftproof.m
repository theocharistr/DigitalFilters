% "Proof" by MATLAB
% A simple technique to develop and verify the steps of a proof 
% using random data input
%
% N P P
% Cornell U 
% Sept 1992
%

clear

n = 8; % any even
x =[-0.1241 + 0.4889i;
   1.4897 + 1.0347i;
   1.4090 + 0.7269i;
   1.4172 - 0.3034i;
   0.6715 + 0.2939i;
  -1.2075 - 0.7873i;
   0.7172 + 0.8884i;
   1.6302 - 1.1471i];
% input 
%x = randn(n,1) + 1i*randn(n,1);
% correct answer
ys = fft(x);

% root of unity
w = @(n,e) exp(-2*pi*1i.*e/n);

k = (0:n-1)';

% DFT proof steps
y = zeros(n,1);

for j = 0:n-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end

fprintf('DFT : %e\n', norm(y - ys))

% split output top bottom
y = zeros(n,1);

for j = 0:n/2-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end
for j = n/2:n-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end

fprintf('split output top bottom : %e\n', norm(y - ys))

% split input even odd
y = zeros(n,1);

k = (0:n/2-1)';
for j = 0:n/2-1
  y(j +1) = sum(w(n,j*2*k) .* x(2*k +1)) + sum(w(n,j*(2*k+1)) .* x(2*k+1 +1));
end
for j = n/2:n-1
  y(j +1) = sum(w(n,j*2*k) .* x(2*k +1)) + sum(w(n,j*(2*k+1)) .* x(2*k+1 +1));
end

fprintf('split input even odd : %e\n', norm(y - ys))

%%completing the proof with the last line
k = (0:n/2-1)';

for j=0:(n/2)-1
 fe(j+1) = sum(w(n/2,j*k) .* x(2*k+1));
 fo(j+1) = sum(w(n/2,j*k) .* x(2*k+1+1));
 
end

wfo = w(n,(0:n/2-1)') .* fo.'; 
y = [fe.' + wfo; fe.' - wfo];

fprintf('done : %e\n', norm(y - ys))
