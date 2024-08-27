function [ out ] = ift2( data )
%2DFT Summary of this function goes here
%   Detailed explanation goes here
[nx, ny, ~] = size(data);

out = fftshift(ifft2(ifftshift(data)));% * sqrt(nx*ny);

end

