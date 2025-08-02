function [imgmat,recongra]= rspsg(img2k,ref20k,phs1,slcs,refs1,refs2,weights)
% recon using spsg method
% call function rspsg(img2k(:,:,:,1),ref20k(:,:,1:2:end,1:6),[-2 -1 0 1 2 3]/2,[6 2])
% call function rspsg(img2k,ref20k,[-2 -1 0 1 2 3]/2,[8 1 3],
%       [0 -2 -4 -6 -1 -3 -5 -7],[0:-1:-7]+[-3 -2 -1 0 1 2 3 4]'*8);
% call function rspsg(img2k,ref20k,[-2 -1 0 1 2 3]/2,[8 1 3],
%       [0 -2 -4],[0:-1:-7]+[-3 -2 -1 0 1 2 3 4]'*8);
% input:    img2k, sms rawdata
%           ref20k, ref rawdata
%           phs1, caipi shift pattern
%           slcs, [sms grappa caipi] factors
%           refs1, for sms index, optional, could also be aligned with
%               ref. if refs1 exists, loop decided by length of refs1
%           refs2, for sms corresponding ref index, optional, sorted as
%               each column refs for each sms element
%               only used for calculating the slice gap
% output:   imgmat, coils*x*y*z, z in the order defined in refs2
%           recongra, recon temp result of one slice of the input

if nargin<5
% refs1 = [0 -2 -4 -6 -8 -10 -1 -3 -5 -7 -9];%sms slices
% refs2 = (0:-1:-10) + [-2 -1 0 1 2 3]'*11;
    refs1 = 0;
    refs2 = (1:slcs(1))';%now change to image index, should be positive
    refsidx = 1;
elseif nargin<6
    refs2 = reshape(1:slcs(1)*length(refs1),slcs(1),length(refs1));
    refsidx = 1:length(refs1);
end
if nargin>=6 
    refsidx = ceil(knnsearch(refs2(:),refs1(:))/slcs(1));% search in ref the sms slices
    refs2 = refs2 - min(refs2(:))+1;%now change to image index, should be positive
end
if nargin<=6
    weights.ker = [3 3];
end
    
if length(slcs)<3
    caipi_factor = 2;
else
    caipi_factor = slcs(3);
end

imgmat = zeros(size(ref20k));% to be recon
for i = 1:length(refs1)
    tmp = mrz.ep2d_caipishift(img2k(:,:,:,i),0,([0:caipi_factor-1])*refs1(i)/(diff(refs2([1 slcs(1)]))/(slcs(1)-1))*2/caipi_factor);%calculate the ratio of 2pi
    [reconspsg,w2] = mrz.slcgra(tmp,mrz.ep2d_caipishift(ref20k(:,:,1:slcs(2):end,1+slcs(1)*(refsidx(i)-1):slcs(1)*refsidx(i)),phs1),weights,'spsg');    
%     reconspsg = mrz.g1.sg.spsg(tmp,mrz.ep2d_caipishift(ref20k(:,:,1:slcs(2):end,1+slcs(1)*(refsidx(i)-1):slcs(1)*refsidx(i)),phs1),[3 3]);
if slcs(2)>1
    recongra = zeros(size(reconspsg).*[1 1 slcs(2) 1]);
    recongra(:,:,1:2:end,:) = reconspsg;
    for j = 1:size(reconspsg,4)
%         recongra(:,:,:,j) = mrz.gra(recongra(:,:,:,j),mrz.ep2d_caipishift(ref20k(:,:,37:96-36,j+(refsidx(i)-1)*slcs(1)),phs1(j)/slcs(2)),[1 slcs(2)],[5 4]);% 
        recongra(:,:,:,j) = mrz.gra(recongra(:,:,:,j),mrz.ep2d_caipishift(ref20k(:,:,:,j+(refsidx(i)-1)*slcs(1)),phs1(j)/slcs(2)),[1 slcs(2)],[5 4]);% 
    end
%     for j = 1:size(reconspsg,4)
%         recongra(:,:,:,j) = mrz.ep2d_grappa(reconspsg(:,:,:,j),[192 96 0],[1 slcs(2)],[5 4],mrz.ep2d_caipishift(ref20k(:,:,:,j),phs1(j)/slcs(2)));
%     end
    imgmat(:,:,:,refs2(:,refsidx(i))) = mrz.ep2d_caipishift(recongra,-phs1/slcs(2));
else
    imgmat(:,:,:,refs2(:,refsidx(i))) = mrz.ep2d_caipishift(reconspsg,-phs1);
    recongra = mrz.ep2d_caipishift(reconspsg,-phs1);
end
end