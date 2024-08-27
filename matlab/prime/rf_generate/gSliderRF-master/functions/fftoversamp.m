function y = fftoversamp(x,n,os,tflag)

% n: length of original filter
% os: oversampling factor

if strcmp(tflag,'transp')
   % inverse fft + decimate
   y = n*os*ifft(x);
   y = y(1:n);
elseif strcmp(tflag, 'notransp')
   % forward oversampled fft
   y = fft(x,os*n);
end