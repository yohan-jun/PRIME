function [out] = ift2_b0( data, Enc_ap, Enc_pa)
%IFT2_BO Summary of this function goes here
%   Detailed explanation goes here


[nx ny ns nc] = size(data);

out = zeros(nx, size(Enc_ap,3), ns, nc) ;
for k = 1:nx
%     disp(size(squeeze(Enc_ap(k,:,:))' * squeeze(data(k,:,1,:))))
    out(k,:,1,:) = squeeze(Enc_ap(k,:,:))' * squeeze(data(k,:,1,:));
    out(k,:,2,:) = squeeze(Enc_pa(k,:,:))' * squeeze(data(k,:,2,:));

end

end