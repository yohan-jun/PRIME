function ref_shift = ep2d_caipishift(ref,freq, phs)
% so far, caipishift, and caipicorr are in the same function, maybe changed
% later, caipicorr used when freq == 0
% caipishift, applied in every slice excited simultaneously
% caipicorr, applied in single slice every phase encoding line
%
%input:     ref, ksp data, coils * adc*phs*slcs
%           freq, size should be same as size(ref,4), for image shift
%           phs, size should be same as sms factor, for correction
%           off-center, then freq == 0, so far.
% todo, caipi corr for 4D data
if nargin<3
    ref_shift = zeros(size(ref));
    for i = 1:length(freq)
        ref_shift(:,:,:,i) = ref(:,:,:,i).*permute(exp(1j*freq(i)*size(ref,3)*(0:size(ref,3)-1)/size(ref,3)*2*pi),[3 1 2]);
    end
elseif freq == 0
    ref_shift = zeros(size(ref));
    
    tmp = repmat(phs(:)*pi,ceil(size(ref,3)/length(phs(:))),1);
    ref_shift = ref.*exp(-1j*permute(tmp(1:size(ref,3)),[3 2 1]));
    
    
end

