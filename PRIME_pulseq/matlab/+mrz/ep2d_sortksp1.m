function [ksp_resample kxx]= ep2d_sortksp1(rawdata0,seq10,dims,fov,os,traj_recon_delay,kx_input)
% ksp calculation from one traj epi sequence, but still uses the resample
% instead of griding, so the gradient needs to be carefully set. 
% if using interp, gradient is not neatly controlled, so using instead 
% mr.makeExtendTrapozoid is better, the k-space is exactly the same as
% define multiple blocks in one EPI train
%
% Call: ksp_resample = ep2d_sortksp1(rawdata0,seq10,dims,fov,os,traj_recon_delay)
% input: rawdata0, data from twx, supposed to be ADC_1traj * Coils * Reps
%        seq, sequence contains the One traj EPI, it may contains several
%           repetitions, so chop to one repetition using the info from dims
%        dims, [Nx Ny Nz], different from ep2d_sortksp be careful, Nz is
%           how many slices contains in the Pulseq file, not necessary to 
%           be the same as Pulseq file Kspace slices
%        fov, used for calculating the kmax, added from ep2d_sortksp, since
%           one traj contains kspace region larger than fov needed
%        os, oversampling rate, default 2
%        traj_recon_delay, delay between odd even. -3e-6 Trio, 1e-6 Prisma,
%           if omitted, then using the calculated delay
%           
% output:ksp, now sorted to coils * Nx * Ny * Nz
%
% rawdata readout samples length is set to divide the Pulseq kspace lines,
% not what dims defines, dims defines how many slices are reconstructed,
% the Pulseq file may contains only one slice or more than acquried slices
%
if nargin < 6
    traj_recon_delay = 0;
end
if nargin < 5 || isempty(os)
    os = 2;
end
Nx = dims(1);
Ny = dims(2);%used for reshape

[ktraj_adc_old, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq10.calculateKspacePP('trajectory_delay', traj_recon_delay);
if(length(dims)<3) dims(3) = size(ktraj_adc_old,2)/size(rawdata0,1);end
if mod(dims(3),1) && dims(3)>1, error('number of slices wrong, please check');end
if dims(3)<1, dims(3)=1;warning('kspace is supposed to be one slice');end
ktraj_adc_tmp = ktraj_adc_old(:,1:size(rawdata0,1));% ktraj_adc to one repetition to save some time
% ktraj_adc_tmp = ktraj_adc_old(:,1:size(ktraj_adc_old,2)/dims(3));%dims(3)
% now is how many slices are needed to be reconstructed
idx1 = find(ktraj_adc_tmp(1,:)>=-1/fov(1)*dims(1)/os/2);
idx2 = find(ktraj_adc_tmp(1,:)<= 1/fov(1)*dims(1)/os/2);
idx = intersect(idx1,idx2);
ktraj_adc = ktraj_adc_tmp(:,idx);
fprintf('resampling dataset %g to %s \n', [length(idx)], [num2str(length(idx)/Ny) '*' num2str(Ny)]);
ksp_resample = zeros(dims(1),size(rawdata0,2),dims(2),size(rawdata0,3));
for I = 1:round(size(rawdata0,3))
    rawdata = permute(reshape(rawdata0(idx,:,I),[],dims(2),size(rawdata0,2)),[1 3 2]);

%%% automatic detection of the measurement parameters (FOV, matrix size, etc)
% not appliable for grappa

nADC = size(rawdata, 1);
k_last=ktraj_adc(:,end);
k_2last=ktraj_adc(:,end-nADC);
delta_ky=k_last(2)-k_2last(2);

Ny_post=round(abs(k_last(2)/delta_ky));
if k_last(2)>0
    Ny_pre=round(abs(min(ktraj_adc(2,:))/delta_ky));
else
    Ny_pre=round(abs(max(ktraj_adc(2,:))/delta_ky));
end


%%% classical phase correction / trajectory delay calculation 
%  here we assume we are dealing with the calibration data
data_odd=ifftshift(ifft(ifftshift(rawdata(:,:,1:2:end),1)),1);
data_even=ifftshift(ifft(ifftshift(rawdata(end:-1:1,:,2:2:end),1)),1);
cmplx_diff=data_even.*conj(data_odd);
cmplx_slope=cmplx_diff(2:end,:,:).*conj(cmplx_diff(1:end-1,:,:));
mslope_phs=angle(sum(cmplx_slope(:)));
dwell_time=(t_adc(nADC)-t_adc(1))/(nADC-1);
measured_traj_delay=mslope_phs/2/2/pi*nADC*dwell_time;
fprintf('measured trajectory delay (assuming it is a calibration data set) is %g s\n', measured_traj_delay);

% we do not calculate the constant phase term here because it depends on
% the definitions of the center of k-space and image-space 
if nargin < 6
    fprintf('using this calculated delay now!\n');
    [ktraj_adc_old, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq10.calculateKspacePP('trajectory_delay', double(measured_traj_delay));
%     ktraj_adc_tmp = ktraj_adc_old(:,1:size(ktraj_adc_old,2)/dims(3));% ktraj_adc to one repetition to save some time
%     idx1 = find(ktraj_adc_tmp(1,1:26100)>=-1/fov(1)/2);
%     idx2 = find(ktraj_adc_tmp(1,1:26100)<= 1/fov(1)/2);
%     idx = intersect(idx1,idx2);
    ktraj_adc = ktraj_adc_old(:,idx);
end

%%% analyze the trajecotory, resample the data
% here we expect rawdata ktraj_adc loaded (and having the same dimensions)
nCoils = size(rawdata, 2); % the incoming data order is [kx coils acquisitions]
nAcq=size(rawdata,3);
nD=size(ktraj_adc, 1);

if nargin==7
    kxx = kx_input;
else
    kxmin=min(ktraj_adc(1,:));
    kxmax=max(ktraj_adc(1,:));
    kxmax1=kxmax/(Nx/2-1)*(Nx/2); % this compensates for the non-symmetric center definition in FFT
    kmaxabs=max(kxmax1, -kxmin);

    kxx= ((-Nx/2):(Nx/2-1))/(Nx/2)*kmaxabs; % kx-sample positions
end
ktraj_adc2=reshape(ktraj_adc,[size(ktraj_adc,1), nADC, size(ktraj_adc,2)/nADC]);
t_adc_chop = t_adc(idx);
t_adc2=reshape(t_adc_chop,[nADC, length(t_adc_chop)/nADC]);

data_resampled=zeros(length(kxx), nCoils, nAcq);
ktraj_resampled=zeros(nD, length(kxx), nAcq);
t_adc_resampled=zeros(length(kxx), nAcq);
for a=1:nAcq
    for c=1:nCoils
        data_resampled(:,c,a)=interp1(ktraj_adc2(1,:,a),rawdata(:,c,a),kxx,'spline',0);
    end
    ktraj_resampled(1,:,a)=kxx;
    for d=2:nD
        ktraj_resampled(d,:,a)=interp1(ktraj_adc2(1,:,a),ktraj_adc2(d,:,a),kxx,'linear',NaN);
    end
    t_adc_resampled(:,a)=interp1(ktraj_adc2(1,:,a),t_adc2(:,a),kxx,'linear',NaN);
end

% figure;
% imagesc(squeeze(abs(data_resampled(:,1,:)))');axis('square');
% title('EPI k-space data');

%%% in some cases (e.g. because of the incorrectly calculated trajectory) phase correction may be needed
%  one such case is the use of the frequency shift proportional to gradient
%  in combination with the gradient delay and FOV offset in the RO direction
%  this calculation is best done with the calibration data, but also seems
%  to work with the actual image data

% here we assume we are dealing with the calibration data
data_odd=ifftshift(ifft(ifftshift(data_resampled(:,:,1:2:end),1)),1);
data_even=ifftshift(ifft(ifftshift(data_resampled(:,:,2:2:end),1)),1);
cmplx_diff1=data_even.*conj(data_odd);
cmplx_diff2=data_even(:,:,1:end-1).*conj(data_odd(:,:,2:end));
mphase1=angle(sum(cmplx_diff1(:)));
mphase2=angle(sum(cmplx_diff2(:)));
mphase=angle(sum([cmplx_diff1(:); cmplx_diff2(:)]));

%%%
pc_coef=0;
%pc_coef=mphase1/2/pi;

data_pc=data_resampled;
for c=1:nCoils
    for i=1:size(data_resampled,1)
        data_pc(i,c,:)=squeeze(data_resampled(i,1,:)).*exp(1i*2*pi*pc_coef*mod((1:size(data_pc,3))',2));
    end
end
% figure;
imagesc(squeeze(angle(data_pc(:,1,:)))');axis('square');
title('phase of hybrid (x/ky) data');

ksp_resample(:,:,:,I) = data_resampled;
end
ksp_resample = permute(ksp_resample,[2 1 3 4]);
end
