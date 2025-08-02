function Nic = filtfilt(ncc, opt, N1, N2, Nc, R)
% Fast computation of zero phase filtering (for alg=4)
fltlen = size(ncc,2)/Nc;    % filter length
numflt = size(ncc,1);       % number of filters

% LORAKS kernel is circular.
% Following indices account for circular elements in a square patch
[in1,in2] = meshgrid(-R:R,-R:R);
idx = find(in1.^2+in2.^2<=R^2);
in1 = in1(idx)';
in2 = in2(idx)';

ind = sub2ind([2*R+1, 2*R+1],R+1+in1,R+1+in2);

filtfilt = zeros([(2*R+1)*(2*R+1),Nc,numflt],'like',ncc);
filtfilt(ind,:,:) = reshape(permute(ncc,[2,1]),[fltlen,Nc,numflt]);
filtfilt = reshape(filtfilt,(2*R+1),(2*R+1),Nc,numflt);

cfilt = conj(filtfilt);

if opt == 'S'       % for S matrix
    ffilt = conj(filtfilt);
else                % for C matrix
    ffilt = flip(flip(filtfilt,1),2);
end

ccfilt = fft2(cfilt,4*R+1, 4*R+1);
fffilt = fft2(ffilt,4*R+1, 4*R+1);

patch = ifft2(sum(reshape(ccfilt,4*R+1,4*R+1,1,Nc,numflt) ...
    .* reshape(fffilt,4*R+1,4*R+1,Nc,1,numflt),5));

if opt == 'S'       % for S matrix
    Nic = fft2(circshift(padarray(patch, [N1-1-2*R N2-1-2*R],'post'),[-4*R-rem(N1,2) -4*R-rem(N2,2)]));
else                % for C matrix
    Nic = fft2(circshift(padarray(patch, [N1-1-2*R N2-1-2*R], 'post'),[-2*R -2*R]));
end
end
