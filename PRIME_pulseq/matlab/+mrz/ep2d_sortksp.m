function ksp_resample = ep2d_sortksp(rawdata0,seq,dims,autokspace,isgrappa,traj_recon_delay)
%
% output: data_resampled, maybe convolution will get better results? but
% now it is just good, permuted as coils * Nx * Ny * Nz
% call example:
%   refk = mrz.ep2d_sortksp(twx.image(),seq,[192 96],0, 0, 1e-6);
%   img2k = mrz.ep2d_sortksp(twx2.image(),seqm2,[192 48],0, 0, 1e-6);
%
% ksp_resample = ep2d_sortksp(rawdata0,seq,dims,autokspace,isgrappa,traj_recon_delay)
% input:    
%       rawdata0:   adc * coils * lins, irectly from mapVBVD
%       seq:        .seq file, to calculate the ksp sampling pattern
%       dims:       [Nx Ny [Nz]], Nx Ny the desired resampled points
%                   now Nz is how many slices contains in 
%                   Pulseq file, saveas ep2d_sortksp1
%       autokspace: 0;%control for auto interpretation of kspace or not,
%                   now always set as 0, just for compatibility
%       isgrappa:   0; also, no grappa pattern calculating here now.
%       traj_recon_delay:
%                   delay between odd and even lines

Nx = dims(1);
phslines = dims(2);Ny = dims(2);
if nargin < 6
    traj_recon_delay = 0;
end
%%% if necessary re-tune the trajectory delay to supress ghosting
% traj_recon_delay=-3e-6;%3.23e-6;%-1e-6;%3.90e-6;%-1.03e-6; % adjust this parameter to supress ghosting (negative allowed) (our trio -1.0e-6, prisma +3.9e-6; avanto +3.88)
[ktraj_adc_old, t_adc_old, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
if length(dims)>2
    ktraj_adc = ktraj_adc_old(:,1:size(ktraj_adc_old,2)/dims(3));
    t_adc = t_adc_old(:,1:size(t_adc_old,2)/dims(3));
else
    ktraj_adc = ktraj_adc_old;
    t_adc = t_adc_old;
end
nADC = size(rawdata0, 1);
if autokspace
    k_last=ktraj_adc(:,end);
    k_2last=ktraj_adc(:,end-nADC);
    delta_ky=k_last(2)-k_2last(2);
    fov=1/abs(delta_ky);
    Ny_post=round(abs(k_last(2)/delta_ky));
    if k_last(2)>0
        Ny_pre=round(abs(min(ktraj_adc(2,:))/delta_ky));
    else
        Ny_pre=round(abs(max(ktraj_adc(2,:))/delta_ky));
    end
    Nx=2*max([Ny_post,Ny_pre]);
    Ny=Nx;
    Ny_sampled=Ny_pre+Ny_post+1;
end

%%% loop for resampling
ksp_resample = zeros(size(rawdata0,2),dims(1),dims(2),size(rawdata0,3)/dims(2));
for I = 1:round(size(rawdata0,3)/phslines)
    rawdata = rawdata0(:,:,1+phslines*(I-1):phslines*I);

%%% automatic detection of the measurement parameters (FOV, matrix size, etc)
% not appliable for grappa



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
    [ktraj_adc_old, t_adc_old] = seq.calculateKspacePP('trajectory_delay', double(measured_traj_delay));
    if length(dims)>2
        ktraj_adc = ktraj_adc_old(:,1:size(ktraj_adc_old,2)/dims(3));
        t_adc = t_adc_old(:,1:size(t_adc_old,2)/dims(3));
    else
        ktraj_adc = ktraj_adc_old;
        t_adc = t_adc_old;
    end
end

%%% analyze the trajecotory, resample the data
% here we expect rawdata ktraj_adc loaded (and having the same dimensions)
nCoils = size(rawdata, 2); % the incoming data order is [kx coils acquisitions]
nAcq=size(rawdata,3);
nD=size(ktraj_adc, 1);

kxmin=min(ktraj_adc(1,:));
kxmax=max(ktraj_adc(1,:));
kxmax1=kxmax/(Nx/2-1)*(Nx/2); % this compensates for the non-symmetric center definition in FFT
kmaxabs=max(kxmax1, -kxmin);

kxx= ((-Nx/2):(Nx/2-1))/(Nx/2)*kmaxabs; % kx-sample positions
ktraj_adc2=reshape(ktraj_adc,[size(ktraj_adc,1), nADC, size(ktraj_adc,2)/nADC]);
t_adc2=reshape(t_adc,[nADC, length(t_adc)/nADC]);

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
imagesc(squeeze(abs(data_resampled(:,1,:)))');axis('square');
title('EPI k-space data');

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

ksp_resample(:,:,:,I) = permute(data_resampled,[2 1 3]);
end
