function [ out ] = ft2( kdata )
%2DFT Summary of this function goes here
%   Detailed explanation goes here
[nx, ny, ~] = size(kdata);
out = fftshift(fft2(ifftshift(kdata)));% / sqrt(nx*ny);

end

