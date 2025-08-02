clear; close all; clc;

coil_comp   = 0; % coil compression
num_cc      = 32; % number of compressed coils
load_ktraj  = 1;
gre_3d      = 1;

diff_set    = 1; % from set1 to set10 for recon.
diff_idx    = 2; % select one of diffusion directions to recon
                 % diff_idx >= 2 (diff_idx = 1 is dummy data)


%% path

set(0,'DefaultFigureWindowStyle','docked')
addpath(genpath('utils/'))


%% load data

file_path       = '/path/to/twix/';
file_name       = 'main_twix_file'; % main twix file
file_sens_name  = 'sens_twix_file'; % low-res GRE for coil sensitivity

seq_path        = '/path/to/seq/';
seq_name        = 'P_PR_1mm_b1k_dir3_shared.seq';

diff_table      =   [0 0 1; 0 1 0; 1 0 0;];
diff_table      =   [0 0 0; 0 0 0; diff_table]; % dummy + b0 + dwi

param           =   load([seq_path, seq_name(1:end-4), '_param.mat']);
if param.adcSamples ~= param.adcSamples_low
    dat_max     = mapVBVD_MAX([file_path, file_name]);
    dat_min     = mapVBVD_MIN([file_path, file_name]);
else
    dat         = mapVBVD_MAX([file_path, file_name]);
end

prot                = [];
prot.alRegridMode   = 1;
tic
fprintf('loading main twix file... ')
if param.adcSamples ~= param.adcSamples_low
    kspace_tmp1 = dat_max{end}.image.unsorted();
    kspace_tmp2 = dat_min{end}.image.unsorted();
else
    kspace_tmp1 = dat{end}.image.unsorted();
end
toc

dat_sens        = mapVBVD_MAX([file_path, file_sens_name]);
tic
fprintf('loading gre twix file... ')
kspace_sens_tmp = dat_sens{end}.image();
toc

tic
fprintf('loading seq files... ')
if load_ktraj == 1
    load([seq_path, seq_name(1:end-4), '_ktraj_adc_acq.mat']);
else
    seq             = mr.Sequence();
    read(seq,[seq_path, seq_name]);
    traj_recon_delay = 0e-6;
    [ktraj_adc_acq, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
    if ~isfile([file_path, seq_name(1:end-4), '_ktraj_adc_acq.mat'])
        save([seq_path, seq_name(1:end-4), '_ktraj_adc_acq.mat'], 'ktraj_adc_acq');
    end
end
toc


%% scan parameters


fov         =   param.fov;
thickness   =   param.thickness;
esp         =   param.actual_esp;       % actual esp
esp_low     =   param.actual_low_esp;   % actual esp

MB          =   param.MB;
Nslices     =   param.Nslices;          % slices
% if MB ~= 1
%     Nslices =   Nslices*MB;
% end
Ngslider    =   length(param.rf);
Nseg        =   param.RSegment;         % shots
R           =   param.R;                % reduction factor (defined in .seq)
R_low       =   param.R_low;
Necho       =   2;                      % echoes
Ndir        =   param.diffusion_count;  % directions
Nc          =   size(kspace_tmp1,2);    % coils
if isfield(param,'Nx_org')
    Nx      =   param.Nx_org;
else
    Nx      =   param.Nx;
end
if isfield(param,'Nx_low_org')
    Nx_low  =   param.Nx_low_org;
else
    Nx_low  =   param.Nx_low;
end
Ny          =   param.Ny;
Ny_low      =   param.Ny_low;
Nnav        =   3;                      % nav lines
Nnav_low    =   3;                      % nav lines
Npe_acq     =   size(kspace_tmp1,3)/Nslices/Ngslider/Nseg/Ndir;
Nfe         =   size(kspace_tmp1,1);
if param.adcSamples ~= param.adcSamples_low
    Nfe_low =   size(kspace_tmp2,1);
else
    Nfe_low =   Nfe;
end
pF          =   param.partFourierFactor;
if isfield(param,'partFourierFactor_fe')
    pF_fe   =   param.partFourierFactor_fe;
else
    pF_fe   =   1;
end
if isfield(param,'partFourierFactor_low_fe')
    pF_low_fe = param.partFourierFactor_low_fe;
else
    pF_low_fe = 1;
end
ro_os       =   param.ro_os;
ro_os_low   =   param.ro_os_low;

Npe_pre         =   param.Ny_pre;
Npe_post        =   param.Ny_post;
Npe_def         =   param.Ny_meas;

Npe_low_pre     =   param.Ny_low_pre;
Npe_low_post    =   param.Ny_low_post;
Npe_low_def     =   param.Ny_low_meas;


%% coil sensitivity

if isfile([file_path, file_sens_name, '_sens.mat'])
    load([file_path, file_sens_name, '_sens.mat']);
    if coil_comp == 1
        load([file_path, file_sens_name, '_sens_cc.mat']);
    end
else
    
    N_tmp               = [Nx,Ny,Nc,Nslices*MB];

    if gre_3d == 1
        kspace_sens_tmp = permute(fft2call(crop(circshift(ifft3call(permute(kspace_sens_tmp,[1,3,4,2])),[0,0,-1,0]),...
            [size(kspace_sens_tmp,1),size(kspace_sens_tmp,3),Nslices*MB,size(kspace_sens_tmp,2)])),[1,4,2,3]);
    end

    kspace_sens_tmp_    = permute(sq(kspace_sens_tmp),[1,3,2,4]);
    N_sens              = size(kspace_sens_tmp_);
    sens_tmp            = crop(ifft2call(kspace_sens_tmp_),[N_sens(1)/2,N_sens(2),N_sens(3),N_sens(4)]);
    kspace_sens_tmp_    = fft2call(sens_tmp);

    kspace_sens_pad     = padarray(kspace_sens_tmp_,N_tmp/2-size(kspace_sens_tmp_)/2);
    if gre_3d ~= 1
        kspace_sens_pad = reshape(kspace_sens_pad,[N_tmp(1:3),1,1,1,1,1,1,N_tmp(4)]);
        kspace_sens_pad = sq(mrir_image_slice_deinterleave(kspace_sens_pad));
    end

    N_low_tmp           = [Nx_low,Ny_low,Nc,Nslices*MB];

    kspace_sens_tmp_    = permute(sq(kspace_sens_tmp),[1,3,2,4]);
    N_sens              = size(kspace_sens_tmp_);
    sens_tmp            = crop(ifft2call(kspace_sens_tmp_),[N_sens(1)/2,N_sens(2),N_sens(3),N_sens(4)]);
    kspace_sens_tmp_    = fft2call(sens_tmp);

    kspace_sens_low_pad = padarray(kspace_sens_tmp_,N_low_tmp/2-size(kspace_sens_tmp_)/2);
    if gre_3d ~= 1
        kspace_sens_low_pad = reshape(kspace_sens_low_pad,[N_low_tmp(1:3),1,1,1,1,1,1,N_low_tmp(4)]);
        kspace_sens_low_pad = sq(mrir_image_slice_deinterleave(kspace_sens_low_pad));
    end

    % coil compression
    if coil_comp == 1
        comp_mtx            = zeros([Nc,Nc,Nslices*MB],'single');
        kspace_sens_pad_cc  = zeros([Nx,Ny,num_cc,Nslices*MB],'single');
        for nn = 1:Nslices*MB
            comp_mtx(:,:,nn) = calcSCCMtx(reshape(kspace_sens_pad(:,:,:,nn),[Nx,Ny,1,Nc]));
            kspace_sens_pad_cc(:,:,:,nn)    = reshape((comp_mtx(:,1:num_cc,nn)'...
                *reshape(kspace_sens_pad(:,:,:,nn),Nx*Ny,Nc).').',Nx,Ny,num_cc);
        end
        
        comp_mtx_low                = zeros([Nc,Nc,Nslices*MB],'single');
        kspace_sens_low_pad_cc      = zeros([Nx_low,Ny_low,num_cc,Nslices*MB],'single');
        for nn = 1:Nslices*MB
            comp_mtx_low(:,:,nn)    = calcSCCMtx(reshape(kspace_sens_low_pad(:,:,:,nn),[Nx_low,Ny_low,1,Nc]));
            kspace_sens_low_pad_cc(:,:,:,nn) = reshape((comp_mtx_low(:,1:num_cc,nn)'...
                *reshape(kspace_sens_low_pad(:,:,:,nn),Nx_low*Ny_low,Nc).').',Nx_low,Ny_low,num_cc);
        end
    end

    num_acs         = 32;
    kernel_size     = [6,6];
    eigen_thresh    = 0.4;

    sens            = zeros(size(kspace_sens_pad));
    sens_low        = zeros(size(kspace_sens_low_pad));

    sens_adp        = zeros(size(kspace_sens_pad));
    sens_low_adp    = zeros(size(kspace_sens_low_pad));

    if coil_comp == 1
        sens_cc     = zeros(size(kspace_sens_pad_cc));
        sens_low_cc = zeros(size(kspace_sens_low_pad_cc));

        sens_adp_cc = zeros(size(kspace_sens_pad_cc));
        sens_low_adp_cc = zeros(size(kspace_sens_low_pad_cc));
    end

    % delete(gcp('nocreate'))
    % c = parcluster('local');
    % total_cores = c.NumWorkers;
    % parpool(min(Nslices*MB,ceil(total_cores/8)))

    % parfor ss = 1:Nslices*MB
    for ss = 1:Nslices*MB
        tic
        fprintf('coil sensitivity estimation... %d/%d ',ss,Nslices*MB);
        [maps, weights]     = ecalib_soft(kspace_sens_pad(:,:,:,ss), num_acs, kernel_size, eigen_thresh);
        sens(:,:,:,ss)      = dot_mult(maps, weights >= eigen_thresh);

        sens_adp_           = estimate_sensitivities(ifft2call(crop(kspace_sens_pad(:,:,:,ss),[128,128,Nc])), 0, 8);
        for cc = 1:Nc
            sens_adp(:,:,cc,ss) = imresize(sens_adp_(:,:,cc),[size(kspace_sens_pad,1) size(kspace_sens_pad,2)]).*rsos(sens(:,:,:,ss),3);
        end

        if coil_comp == 1
            [maps, weights]     = ecalib_soft(kspace_sens_pad_cc(:,:,:,ss), num_acs, kernel_size, eigen_thresh);
            sens_cc(:,:,:,ss)   = dot_mult(maps, weights >= eigen_thresh);

            sens_adp_cc_        = estimate_sensitivities(ifft2call(crop(kspace_sens_pad_cc(:,:,:,ss),[128,128,num_cc])), 0, 8);
            for cc = 1:num_cc
                sens_adp_cc(:,:,cc,ss) = imresize(sens_adp_cc_(:,:,cc),[size(kspace_sens_pad_cc,1) size(kspace_sens_pad_cc,2)]).*rsos(sens_cc(:,:,:,ss),3);
            end
        end
        toc
    end

    % parfor ss = 1:Nslices*MB
    for ss = 1:Nslices*MB
        tic
        fprintf('coil sensitivity estimation... %d/%d ',ss,Nslices*MB);
        [maps_low, weights_low] = ecalib_soft(kspace_sens_low_pad(:,:,:,ss), num_acs, kernel_size, eigen_thresh);
        sens_low(:,:,:,ss)      = dot_mult(maps_low, weights_low >= eigen_thresh);

        sens_low_adp(:,:,:,ss)  = abs(sens_low(:,:,:,ss)) .* exp(1i * angle(estimate_sensitivities(ifft2call(kspace_sens_low_pad(:,:,:,ss)), 0, 8).*rsos(sens_low(:,:,:,ss),3)));

        if coil_comp == 1
            [maps_low, weights_low] = ecalib_soft(kspace_sens_low_pad_cc(:,:,:,ss), num_acs, kernel_size, eigen_thresh);
            sens_low_cc(:,:,:,ss)   = dot_mult(maps_low, weights_low >= eigen_thresh);

            sens_low_adp_cc(:,:,:,ss) = estimate_sensitivities(ifft2call(kspace_sens_low_pad_cc(:,:,:,ss)), 0, 8).*rsos(sens_low_cc(:,:,:,ss),3);
        end
        toc
    end
    delete(gcp('nocreate'))

    sens_adp = abs(sens) .* exp(1i * angle(sens_adp));

    sens    = crop(sens,N_tmp);
    sens_adp = crop(sens_adp,N_tmp);
    sens    = flip(permute(sens,[1,2,4,3]),1);
    sens_adp = flip(permute(sens_adp,[1,2,4,3]),1);
    if coil_comp == 1
        sens_adp_cc = abs(sens_cc) .* exp(1i * angle(sens_adp_cc));

        N_tmp_cc = [Nx,Ny,num_cc,Nslices*MB];
        sens_cc = crop(sens_cc,N_tmp_cc);
        sens_adp_cc = crop(sens_adp_cc,N_tmp_cc);
        sens_cc = flip(permute(sens_cc,[1,2,4,3]),1);
        sens_adp_cc = flip(permute(sens_adp_cc,[1,2,4,3]),1);
    end

    sens_low_adp = abs(sens_low_adp) .* exp(1i * angle(sens_low_adp));

    sens_low    = crop(sens_low,N_low_tmp);
    sens_low_adp = crop(sens_low_adp,N_low_tmp);
    sens_low    = flip(permute(sens_low,[1,2,4,3]),1);
    sens_low_adp = flip(permute(sens_low_adp,[1,2,4,3]),1);
    if coil_comp == 1
        sens_low_adp_cc = abs(sens_low_cc) .* exp(1i * angle(sens_low_adp_cc));

        N_low_tmp_cc = [Nx_low,Ny_low,num_cc,Nslices*MB];
        sens_low_cc = crop(sens_low_cc,N_low_tmp_cc);
        sens_low_adp_cc = crop(sens_low_adp_cc,N_low_tmp_cc);
        sens_low_cc = flip(permute(sens_low_cc,[1,2,4,3]),1);
        sens_low_adp_cc = flip(permute(sens_low_adp_cc,[1,2,4,3]),1);
    end

    save([file_path, file_sens_name, '_sens.mat'], 'sens','sens_low','sens_adp','sens_low_adp','-v7.3');
    if coil_comp == 1
        save([file_path, file_sens_name, '_sens_cc.mat'], 'sens_cc','sens_low_cc','sens_adp_cc','sens_low_adp_cc','-v7.3');
    end
end
clear sens sens_cc sens_low sens_low_cc % memory opt. (if OPT2 is used later)
% clear sens_adp sens_adp_cc sens_low_adp sens_low_adp_cc % memory opt. (if OPT1 is used later)


%% Reconstruction

img_fieldmap            = zeros([Nx,Ny,Nslices*MB,Ndir,Ngslider]);
img_low_fieldmap        = zeros([Nx_low,Ny_low,Nslices*MB,Ndir,Ngslider]);
img_low_pad_fieldmap    = zeros([Nx,Ny,Nslices*MB,Ndir,Ngslider]);
sense_recon             = zeros([Nx,Ny,Nseg,Nslices*MB,Ndir,Ngslider]);
sense_low_recon         = zeros([Nx_low,Ny_low,Nseg,Nslices*MB,Ndir,Ngslider]);
hybrid_sense_1          = zeros([Nx,Ny,Nslices*MB,Ndir,Ngslider]);
hybrid_sense_2          = zeros([Nx,Ny,Nslices*MB,Ndir,Ngslider]);
hybrid_sense_low        = zeros([Nx_low,Ny_low,Nslices*MB,Ndir,Ngslider]);
sense_b0_1              = zeros([Nx,Ny,Nseg,Nslices*MB,Ndir,Ngslider]);
sense_b0_2              = zeros([Nx,Ny,Nseg,Nslices*MB,Ndir,Ngslider]);
sense_low_b0            = zeros([Nx_low,Ny_low,Nseg,Nslices*MB,Ndir,Ngslider]);
BUDA_1                  = zeros([Nx,Ny,Nseg,Nslices*MB,Ndir,Ngslider]);
BUDA_2                  = zeros([Nx,Ny,Nseg,Nslices*MB,Ndir,Ngslider]);
BUDA_low                = zeros([Nx_low,Ny_low,Nseg,Nslices*MB,Ndir,Ngslider]);
BUDA_LORAKS_1           = zeros([Nx,Ny,Nseg,Nslices*MB,Ndir,Ngslider]);
BUDA_LORAKS_2           = zeros([Nx,Ny,Nseg,Nslices*MB,Ndir,Ngslider]);
BUDA_LORAKS_low         = zeros([Nx_low,Ny_low,Nseg,Nslices*MB,Ndir,Ngslider]);


%%

for ndwi = diff_idx:diff_idx % (index 1 : dummy scan)

    fprintf('\n===== =============== =====\n');
    fprintf('        Diff Dir: %d/%d',ndwi,Ndir);
    fprintf('\n===== =============== =====\n');

    %% k-space indexing

    tstart = tic;
    fprintf('\n===== =============== =====\n');
    fprintf('        Diff Dir: %d/%d',ndwi,Ndir);
    fprintf('\n===== =============== =====\n');
    fprintf('k-space indexing... ');

    Nc = size(kspace_tmp1,2);

    kspace          = zeros([floor(Nx*pF_fe*ro_os),Nc,Ny*pF,Necho,Nseg,Ngslider,Nslices,2]);
    kspace_low      = zeros([floor(Nx_low*pF_low_fe*ro_os_low),Nc,floor(Ny_low*pF),Necho,Nseg,Ngslider,Nslices,2]);
    kspace_nav      = zeros([floor(Nx*pF_fe*ro_os),Nc,Nnav,Necho,Nseg,Ngslider,Nslices,2]);
    kspace_low_nav  = zeros([floor(Nx_low*pF_low_fe*ro_os_low),Nc,Nnav_low,Necho,Nseg,Ngslider,Nslices,2]);

    if param.adcSamples ~= param.adcSamples_low
        ktraj_adc1  = reshape(ktraj_adc_acq(1,1:Nfe*Npe_def),[1,Nfe,Npe_def]);
        ktraj_adc2  = reshape(ktraj_adc_acq(1,Nfe*Npe_def+1:Nfe*Npe_def+Nfe_low*Npe_low_def),[1,Nfe_low,Npe_low_def]);
        ktraj_adc2  = cat(3,zeros(size(ktraj_adc2,1),size(ktraj_adc2,2),Npe_def),ktraj_adc2);
    else
        ktraj_adc1  = reshape(ktraj_adc_acq,[size(ktraj_adc_acq,1),Nfe,size(ktraj_adc_acq,2)/Nfe]);
        ktraj_adc2  = ktraj_adc1;
    end

    cnt = 0;
    % order should be dependent on the .seq definition
    for dd = 1:Ndir
        for gg = 1:Nseg
            for ll = 1:Ngslider
                for ss = 1:Nslices
                    tic
                    for ee = 1:Necho
                        if dd == 1 || dd == ndwi
                            if param.adcSamples ~= param.adcSamples_low
                                kspace_tmp_org1 = kspace_tmp1(:,:,Npe_acq*cnt+1:Npe_acq*(cnt+1));
                                kspace_tmp_org2 = kspace_tmp2(:,:,Npe_acq*cnt+1:Npe_acq*(cnt+1));
                            else
                                kspace_tmp_org1 = kspace_tmp1(:,:,Npe_acq*cnt+1:Npe_acq*(cnt+1));
                                kspace_tmp_org2 = kspace_tmp_org1;
                            end

                            ktraj_tmp1 = ktraj_adc1(1,:,1:Npe_def);
                            ktraj_tmp2 = ktraj_adc2(1,:,1:Npe_acq);

                            for cc = 1:Nc
                                for aa = 1:Npe_acq
                                    if aa <= Npe_def
                                        kxmin   = min(min(sq(ktraj_tmp1(1,:,1:Npe_def))))*(Nx*pF_fe*ro_os/2-1)/(Nx*pF_fe*ro_os/2);
                                        kxmax   = max(max(sq(ktraj_tmp1(1,:,1:Npe_def))))/(Nx*pF_fe*ro_os/2-1)*(Nx*pF_fe*ro_os/2);
                                        kxx     = ((-Nx*pF_fe*ro_os/2):(Nx*pF_fe*ro_os/2-1))/(Nx*pF_fe*ro_os/2);
                                        if floor(Nx*pF_fe*ro_os/2) == Nx*pF_fe*ro_os/2
                                            kxx     = (kxx+1)/2*(kxmax-kxmin)+kxmin;
                                        else
                                            kxx     = (kxx-min(kxx))/(max(kxx)-min(kxx))*(kxmax-kxmin)+kxmin;
                                        end
                                        if mod(aa,2) == 0
                                            kxx = flip(kxx);
                                        end
                                        kspace_tmp_(:,cc,aa) = crop(interp1(ktraj_tmp1(1,:,aa),kspace_tmp_org1(:,cc,aa),kxx,'spline',0),[1,size(kspace,1)]);
                                    else
                                        if R_low == 1
                                            kxmin   = min(min(sq(crop(ktraj_tmp2(1,:,Npe_def+1:Npe_def+Npe_low_def),[1,Nfe_low*0.91,Npe_low_def]))))*(Nx_low*pF_low_fe*ro_os_low/2-1)/(Nx_low*pF_low_fe*ro_os_low/2);
                                            kxmax   = max(max(sq(crop(ktraj_tmp2(1,:,Npe_def+1:Npe_def+Npe_low_def),[1,Nfe_low*0.91,Npe_low_def]))))/(Nx_low*pF_low_fe*ro_os_low/2-1)*(Nx_low*pF_low_fe*ro_os_low/2);
                                        elseif R_low == 1.5
                                            kxmin   = min(min(sq(crop(ktraj_tmp2(1,:,Npe_def+1:Npe_def+Npe_low_def),[1,Nfe_low*0.95,Npe_low_def]))))*(Nx_low*pF_low_fe*ro_os_low/2-1)/(Nx_low*pF_low_fe*ro_os_low/2);
                                            kxmax   = max(max(sq(crop(ktraj_tmp2(1,:,Npe_def+1:Npe_def+Npe_low_def),[1,Nfe_low*0.95,Npe_low_def]))))/(Nx_low*pF_low_fe*ro_os_low/2-1)*(Nx_low*pF_low_fe*ro_os_low/2);
                                        elseif R_low == 2
                                            kxmin   = min(min(sq(crop(ktraj_tmp2(1,:,Npe_def+1:Npe_def+Npe_low_def),[1,Nfe_low*0.95,Npe_low_def]))))*(Nx_low*pF_low_fe*ro_os_low/2-1)/(Nx_low*pF_low_fe*ro_os_low/2);
                                            kxmax   = max(max(sq(crop(ktraj_tmp2(1,:,Npe_def+1:Npe_def+Npe_low_def),[1,Nfe_low*0.95,Npe_low_def]))))/(Nx_low*pF_low_fe*ro_os_low/2-1)*(Nx_low*pF_low_fe*ro_os_low/2);
                                        end
                                        
                                        kxx     = ((-Nx_low*pF_low_fe*ro_os_low/2):(Nx_low*pF_low_fe*ro_os_low/2-1))/(Nx_low*pF_low_fe*ro_os_low/2);
                                        if floor(Nx_low*pF_low_fe*ro_os_low/2) == Nx_low*pF_low_fe*ro_os_low/2
                                            kxx     = (kxx+1)/2*(kxmax-kxmin)+kxmin;
                                        else
                                            kxx     = (kxx-min(kxx))/(max(kxx)-min(kxx))*(kxmax-kxmin)+kxmin;
                                        end
                                        if mod(aa,2) == mod(Npe_def,2)
                                            kxx = flip(kxx);
                                        end
                                        kspace_low_tmp_(:,cc,aa) = interp1(ktraj_tmp2(1,:,aa),kspace_tmp_org2(:,cc,aa),kxx,'spline',0);
                                    end
                                end
                            end

                            if Nnav ~= 0 && Nnav_low ~= 0
                                kspace_nav(:,:,1,ee,gg,ll,ss,1) = kspace_tmp_(:,:,1);
                                kspace_nav(:,:,2,ee,gg,ll,ss,2) = flip(kspace_tmp_(:,:,2),1);
                                kspace_nav(:,:,3,ee,gg,ll,ss,1) = kspace_tmp_(:,:,3);
                                kspace_low_nav(:,:,1,ee,gg,ll,ss,1) = kspace_low_tmp_(:,:,Npe_def+1);
                                kspace_low_nav(:,:,2,ee,gg,ll,ss,2) = flip(kspace_low_tmp_(:,:,Npe_def+2),1);
                                kspace_low_nav(:,:,3,ee,gg,ll,ss,1) = kspace_low_tmp_(:,:,Npe_def+3);
                            end

                            if gg == 1
                                if ee == 1
                                    k_idx_1 = (1+2:R*Nseg*2:pF*Ny);
                                    k_idx_2 = (1+2+R*Nseg:R*Nseg*2:pF*Ny);
                                    kspace(:,:,k_idx_1,ee,gg,ll,ss,1) = kspace_tmp_(:,:,1:2:Npe_def);
                                    kspace(:,:,k_idx_2,ee,gg,ll,ss,2) = flip(kspace_tmp_(:,:,2:2:Npe_def),1);
                                else
                                    k_idx_1 = (1:R_low*Nseg*2:pF*Ny_low);
                                    k_idx_2 = (1+R_low*Nseg:R_low*Nseg*2:pF*Ny_low);
                                    kspace_low(:,:,k_idx_1,ee,gg,ll,ss,1) = kspace_low_tmp_(:,:,Npe_def+1:2:Npe_def+Npe_low_def-1);
                                    kspace_low(:,:,k_idx_2,ee,gg,ll,ss,2) = flip(kspace_low_tmp_(:,:,Npe_def+2:2:Npe_def+Npe_low_def-2),1);
                                end
                            elseif gg == 2
                                if ee == 1
                                    k_idx_1 = pF*Ny-R*Nseg:-R*Nseg*2:1;
                                    k_idx_2 = pF*Ny:-R*Nseg*2:1;
                                    kspace(:,:,k_idx_2,ee,gg,ll,ss,1) = kspace_tmp_(:,:,1:2:Npe_def);
                                    kspace(:,:,k_idx_1,ee,gg,ll,ss,2) = flip(kspace_tmp_(:,:,2:2:Npe_def-1),1);
                                else
                                    k_idx_1 = Npe_low_def*R_low*Nseg-R_low*Nseg-R_low*2-3:-R_low*Nseg*2:1;
                                    k_idx_2 = Npe_low_def*R_low*Nseg-R_low*2-3:-R_low*Nseg*2:1;
                                    kspace_low(:,:,k_idx_2(1:length(Npe_def+1:2:Npe_def+Npe_low_def-1)),ee,gg,ll,ss,1) = kspace_low_tmp_(:,:,Npe_def+1:2:Npe_def+Npe_low_def-1);
                                    kspace_low(:,:,k_idx_1(1:length(Npe_def+2:2:Npe_def+Npe_low_def-2)),ee,gg,ll,ss,2) = flip(kspace_low_tmp_(:,:,Npe_def+2:2:Npe_def+Npe_low_def-2),1);
                                end
                            end
                        end
                    end
                    cnt = cnt + 1;
                    fprintf('k-space indexing... %d/%d ',cnt,Ndir*Nseg*Ngslider*Nslices);
                    toc
                end
            end
        end
        if dd == 1
            kspace_nav_ = kspace_nav;
            kspace_low_nav_ = kspace_low_nav;
        elseif dd == ndwi
            kspace_nav = kspace_nav_;
            kspace_low_nav = kspace_low_nav_;
        end
    end
    
    % crop k-space
    kspace_size_tmp         = size(kspace);
    kspace_low_size_tmp     = size(kspace_low);
    kspace_nav_size_tmp     = size(kspace_nav);
    kspace_low_nav_size_tmp = size(kspace_low_nav);

    clear kspace_tmp1 kspace_tmp2 % for single direction only (to reduce memory)
    kspace = ifftc(kspace,1); kspace = ifftc(kspace,3); % memory opt.
    kspace = crop(kspace,[Nx*pF_fe,kspace_size_tmp(2:end)]); % memory opt.
    kspace = fftc(kspace,1);  kspace = fftc(kspace,3); % memory opt.

    kspace_low = ifftc(kspace_low,1); kspace_low = ifftc(kspace_low,3); % memory opt.
    kspace_low = crop(kspace_low,[Nx_low*pF_low_fe,kspace_low_size_tmp(2:end)]); % memory opt.
    kspace_low = fftc(kspace_low,1);  kspace_low = fftc(kspace_low,3); % memory opt.

    kspace_nav_acq          = fftc(crop(ifftc(kspace_nav,1),[Nx*pF_fe,kspace_nav_size_tmp(2:end)]),1);
    kspace_low_nav_acq      = fftc(crop(ifftc(kspace_low_nav,1),[Nx_low*pF_low_fe,kspace_low_nav_size_tmp(2:end)]),1);
    tend = toc(tstart)


    %% k-space order

    tic
    fprintf('\n===== =============== =====\n');
    fprintf('        Diff Dir: %d/%d',ndwi,Ndir);
    fprintf('\n===== =============== =====\n');
    fprintf('k-space ordering... ');
    % twix order
    % ['Col','Lin','Cha','Set','Eco','Phs','Rep','Seg','Par','Sli','Ave'];
    % Rep : directions, Seg: 2 (even/odd)
    % Unknown: 'Set','Phs','Par'

    kspace_nav      = kspace_nav_acq;
    kspace_low_nav  = kspace_low_nav_acq;
    clear kspace_acq kspace_low_acq kspace_nav_acq kspace_low_nav_acq

    kspace = permute(kspace,[1,3,2,9,4,10,5,8,6,7]); % memory opt.
    clear kspace_

    kspace_low_ = zeros([Nx_low*pF_low_fe,floor(Ny_low*pF),Nc,1,Necho,1,Nseg,2,Ngslider,Nslices]);
    
    kspace_low_(:,:,:,1,:,1,:,:,:,:) = permute(kspace_low,[1,3,2,4,5,8,6,7]);
    kspace_low  = kspace_low_;
    clear kspace_low_

    kspace_nav_ = zeros([Nx*pF_fe,Nnav,Nc,1,Necho,1,Nseg,2,Ngslider,Nslices]);
    kspace_nav_(:,:,:,1,:,1,:,:,:,:) = permute(kspace_nav,[1,3,2,4,5,8,6,7]);
    kspace_nav_ = kspace_nav_(:,:,:,:,1,:,:,:,:,:);
    kspace_nav  = kspace_nav_;
    clear kspace_nav_

    kspace_low_nav_ = zeros([Nx_low*pF_low_fe,Nnav,Nc,1,Necho,1,Nseg,2,Ngslider,Nslices]);
    kspace_low_nav_(:,:,:,1,:,1,:,:,:,:) = permute(kspace_low_nav,[1,3,2,4,5,8,6,7]);
    kspace_low_nav_ = kspace_low_nav_(:,:,:,:,1,:,:,:,:,:);
    kspace_low_nav  = kspace_low_nav_;
    clear kspace_low_nav_
    toc


    %% deinterleave slices

    tic
    fprintf('\n===== =============== =====\n');
    fprintf('        Diff Dir: %d/%d',ndwi,Ndir);
    fprintf('\n===== =============== =====\n');
    fprintf('deinterleave slices... ');
    AccY            = R*Nseg;
    AccY_low        = R_low*Nseg;
    AccZ            = 2;
    PhaseShiftBase  = pi; % pi

    kspace          = mrir_image_slice_deinterleave(kspace); % memory opt.
    if mod(Nslices,2) == 1
        kspace      = circshift(kspace,1,10); % memory opt.
    end
    clear kspace_

    kspace_low_     = mrir_image_slice_deinterleave(kspace_low);
    if mod(Nslices,2) == 1
        kspace_low_ = circshift(kspace_low_,1,10);
    end
    kspace_low      = kspace_low_;
    clear kspace_low_

    kspace_nav_     = mrir_image_slice_deinterleave(kspace_nav);
    kspace_nav      = kspace_nav_;
    clear kspace_nav_

    kspace_low_nav_ = mrir_image_slice_deinterleave(kspace_low_nav);
    kspace_low_nav  = kspace_low_nav_;
    clear kspace_nav_
    toc


    %% ghost correction

    tic
    fprintf('\n===== =============== =====\n');
    fprintf('        Diff Dir: %d/%d',ndwi,Ndir);
    fprintf('\n===== =============== =====\n');
    fprintf('ghost correction... ');

    for rr = 1:Nseg % AP/PA

        % ghost correction
        kspace_cor      = zeross(size(kspace(:,:,:,:,:,:,rr,1,:,:)));
        kspace_low_cor  = zeross(size(kspace_low(:,:,:,:,:,:,rr,1,:,:)));

        for ll = 1:Ngslider
            for ee = 1:Necho
                if ee == 1
                    [kspace_cor(:,:,:,:,ee,:,:,:,ll,:)]      = ghost_correct_pat_ref_v1_STD_bb(prot, kspace(:,:,:,:,ee,:,rr,:,ll,:), kspace_nav(:,:,:,:,:,:,rr,:,ll,:), AccY);
                else
                    [kspace_low_cor(:,:,:,:,ee,:,:,:,ll,:)]  = ghost_correct_pat_ref_v1_STD_bb(prot, kspace_low(:,:,:,:,ee,:,rr,:,ll,:), kspace_low_nav(:,:,:,:,:,:,rr,:,ll,:), AccY_low);
                end
            end
        end

        % deblur due to peshift (for sms only)
        if PhaseShiftBase
            kspace_cor_deblur       = caipi_deblur_CLv1_YJ(kspace_cor, prot, AccZ, PhaseShiftBase, 1);
            kspace_low_cor_deblur   = caipi_deblur_CLv1_YJ(kspace_low_cor, prot, AccZ, PhaseShiftBase, 2);
        else
            kspace_cor_deblur       = kspace_cor;
            kspace_low_cor_deblur   = kspace_low_cor;
        end

        % zero pad due to partial fourier
        pad_factor          = Ny*(1-pF);
        pad_factor_low      = round(Ny_low*(1-pF));

        pad_factor_fe       = Nx*(1-pF_fe);
        pad_factor_low_fe   = Nx_low*(1-pF_low_fe);

        if rr == 1
            kspace_cor_full_ap      = sq(padarray(kspace_cor_deblur,[0,pad_factor],'pre'));
            kspace_cor_full_ap      = padarray(kspace_cor_full_ap,[pad_factor_fe,0],'post');
            kspace_low_cor_full_ap  = sq(padarray(kspace_low_cor_deblur,[0,pad_factor_low],'post'));
            kspace_low_cor_full_ap  = padarray(kspace_low_cor_full_ap,[pad_factor_low_fe,0],'post');
        else
            kspace_cor_full_pa      = sq(padarray(kspace_cor_deblur,[0,pad_factor],'post'));
            kspace_cor_full_pa      = padarray(kspace_cor_full_pa,[pad_factor_fe,0],'post');
            kspace_low_cor_full_pa  = sq(padarray(kspace_low_cor_deblur,[0,pad_factor_low],'pre'));
            kspace_low_cor_full_pa  = padarray(kspace_low_cor_full_pa,[pad_factor_low_fe,0],'post');
        end
    end
    clear kspace_cor kspace_low_cor kspace_ kspace_low kspace_cor_deblur kspace_low_cor_deblur kspace_nav kspace_low_nav
    toc


    %% apply delta_ky shift

    tic
    fprintf('\n===== =============== =====\n');
    fprintf('        Diff Dir: %d/%d',ndwi,Ndir);
    fprintf('\n===== =============== =====\n');
    fprintf('delta ky shift... ');

    if AccY == 3
        del_ky_shift    = [1,1];
    elseif AccY == 4
        del_ky_shift    = [3,0];
    elseif AccY == 5
        del_ky_shift        = [1,-3];
    end

    if AccY_low == 2
        del_ky_low_shift    = [0,1];
    elseif AccY_low == 3
        del_ky_low_shift    = [2,1];
    elseif AccY_low == 4
        del_ky_low_shift    = [3,-1];
    end

    % correct k-space positioning:
    del_ky          = mod(AccY-del_ky_shift,AccY); % delta_ky to use in recon after circshift

    kspace_ap       = circshift(kspace_cor_full_ap,[0,-del_ky_shift(1)]);
    clear kspace_cor_full_ap
    kspace_pa       = circshift(kspace_cor_full_pa,[0,-del_ky_shift(2)]);
    clear kspace_cor_full_pa

    Img_epi_ap      = ifft2call(kspace_ap);
    kspace_ap_tmp   = kspace_ap(:,:,1,1,1,1);
    clear kspace_ap
    Img_epi_pa      = ifft2call(kspace_pa);
    kspace_pa_tmp   = kspace_pa(:,:,1,1,1,1);
    clear kspace_pa

    kspace_low_ap   = circshift(kspace_low_cor_full_ap,[0,-del_ky_low_shift(1)]);
    clear kspace_low_cor_full_ap
    kspace_low_pa   = circshift(kspace_low_cor_full_pa,[0,-del_ky_low_shift(2)]);
    clear kspace_low_cor_full_pa

    Img_epi_low_ap  = ifft2call(kspace_low_ap);
    kspace_low_ap_tmp = kspace_low_ap(:,:,1,2,1,1);
    clear kspace_low_ap
    Img_epi_low_pa  = ifft2call(kspace_low_pa);
    kspace_low_pa_tmp = kspace_low_pa(:,:,1,2,1,1);
    clear kspace_low_pa

    toc


    %% forward model for subsampled data

    tic
    fprintf('\n===== =============== =====\n');
    fprintf('        Diff Dir: %d/%d',ndwi,Ndir);
    fprintf('\n===== =============== =====\n');
    fprintf('forward model for subsampled data... ');

    clear N

    [N(1), N(2), num_chan, num_echo, num_gslider, num_slice] = size(Img_epi_ap);

    num_dwi     = 1;
    Ry          = AccY;

    if AccY == 3
        del_ky_a    = [2,1];
    elseif AccY == 4
        del_ky_a    = [3,1];
    elseif AccY == 5
        del_ky_a    = [4,2];
    end

    A_1 = fftc(eye(N(2)),1);

    A_ap = A_1(1+del_ky_a(1):Ry:end,:);
    A_pa = A_1(1+del_ky_a(2):Ry:end,:);

    % (x,y,chan,ap/pa,echo,gslider,slc)
    img_patref = cat(4,permute(Img_epi_ap,[1,2,3,7,4,5,6]),permute(Img_epi_pa,[1,2,3,7,4,5,6]));
    clear Img_epi_ap Img_epi_pa % memory opt.

    % (x, ky, chan) subsampled data
    sgnl_ap = zeross([round(N ./ [1,Ry]), num_chan, num_echo, num_gslider, num_slice]);
    sgnl_pa = zeross([round(N ./ [1,Ry]), num_chan, num_echo, num_gslider, num_slice]);
    
    for nslc = 1:num_slice
        for ng = 1:num_gslider
            for ne = 1:num_echo
                for xn = 1:N(1)
                    m_ap = sq( img_patref(xn,:,:,1,ne,ng,nslc) );
                    m_pa = sq( img_patref(xn,:,:,2,ne,ng,nslc) );

                    % signal in x,ky,chan
                    sgnl_ap(xn, :, :, ne, ng, nslc) = A_ap  * sq(m_ap);
                    sgnl_pa(xn, :, :, ne, ng, nslc) = A_pa  * sq(m_pa);
                end
            end
        end
    end
    clear img_patref % memory opt.


    [N_low(1), N_low(2), num_chan, num_echo, num_gslider, num_slice] = size(Img_epi_low_ap);

    Ry_low          = AccY_low;
    del_ky_low_a    = [0,1];
    
    A_low = fftc(eye(N_low(2)),1);

    A_low_ap = A_low(1+del_ky_low_a(1):Ry_low:end,:);
    A_low_pa = A_low(1+del_ky_low_a(2):Ry_low:end,:);

    % (x,y,chan,ap/pa,echo,slc)
    img_low_patref = cat(4,permute(Img_epi_low_ap,[1,2,3,7,4,5,6]),permute(Img_epi_low_pa,[1,2,3,7,4,5,6]));
    clear Img_epi_low_ap Img_epi_low_pa % memory opt.

    % (x, ky, chan) subsampled data
    sgnl_low_ap = zeross([round(N_low ./ [1,Ry_low]), num_chan, num_echo, num_gslider, num_slice]);
    sgnl_low_pa = zeross([round(N_low ./ [1,Ry_low]), num_chan, num_echo, num_gslider, num_slice]);

    for nslc = 1:num_slice
        for ng = 1:num_gslider
            for ne = 1:num_echo
                for xn = 1:N_low(1)
                    m_low_ap = sq( img_low_patref(xn,:,:,1,ne,ng,nslc) );
                    m_low_pa = sq( img_low_patref(xn,:,:,2,ne,ng,nslc) );

                    % signal in x,ky,chan
                    sgnl_low_ap(xn, :, :, ne, ng, nslc) = A_low_ap  * sq(m_low_ap);
                    sgnl_low_pa(xn, :, :, ne, ng, nslc) = A_low_pa  * sq(m_low_pa);
                end
            end
        end
    end

    clear img_low_patref % memory opt.
    toc


    %% Sense for AP and PA separately without B0 model

    tic
    fprintf('\n===== =============== =====\n');
    fprintf('        Diff Dir: %d/%d',ndwi,Ndir);
    fprintf('\n===== =============== =====\n');
    fprintf('SENSE without B0... ');

    lambda_tik = 1e-5;

    img_sense_ap = zeross([N(1), N(2)*MB, num_echo, num_gslider, num_slice]);
    img_sense_pa = zeross([N(1), N(2)*MB, num_echo, num_gslider, num_slice]);

    num_read = N(1);

    if MB == 2

    shift_amount = (ceil(N(2)/(AccY*MB))*[0,1])*1;

    % parfor nslc = 1:num_slice
    for nslc = 1:num_slice
        disp(['slice: ', num2str(nslc), ' / ', num2str(num_slice)])

        sens_select = zeros(N(1),N(2),MB,num_chan);
        nslc_select = Nslices*MB-nslc+1:-Nslices:1;

        for zz = 1:MB
            % sens_select(:,:,zz,:) = circshift(sens(:,:,nslc_select(zz),:), shift_amount(zz), 2); % OPT 1 org
            sens_select(:,:,zz,:) = circshift(sens_adp(:,:,nslc_select(zz),:), shift_amount(zz), 2); % OPT 2 update: coil sens
        end

        for xn = 1:num_read
            AC_ap = zeross([num_chan * round(N(2) / Ry), N(2)*MB]);
            AC_pa = zeross([num_chan * floor(N(2) / Ry), N(2)*MB]);

            for c = 1:num_chan
                AC_ap(1 + (c-1)*round(N(2)/Ry) : c*round(N(2)/Ry), :) = repmat(A_ap,1,MB) * diag( reshape(permute(sens_select(xn,:,:,c),[2,3,1]),[N(2)*MB,1]) );
                AC_pa(1 + (c-1)*floor(N(2)/Ry) : c*floor(N(2)/Ry), :) = repmat(A_pa,1,MB) * diag( reshape(permute(sens_select(xn,:,:,c),[2,3,1]),[N(2)*MB,1]) );
            end

            [U_ap, S_ap, V_ap] = svd(AC_ap, 'econ');
            [U_pa, S_pa, V_pa] = svd(AC_pa, 'econ');

            Einv_ap = V_ap * diag(diag(S_ap) ./ (diag(S_ap).^2 + lambda_tik)) * U_ap';
            Einv_pa = V_pa * diag(diag(S_pa) ./ (diag(S_pa).^2 + lambda_tik)) * U_pa';

            for ng = 1:num_gslider
                for ne = 1:num_echo
                    rhs_ap = sgnl_ap(xn, :, :, ne, ng, nslc);
                    rhs_ap = rhs_ap(:);
                    img_sense_ap(xn, :, ne, ng, nslc) = Einv_ap * rhs_ap;

                    rhs_pa = sgnl_pa(xn, :, :, ne, ng, nslc);
                    rhs_pa = rhs_pa(:);
                    img_sense_pa(xn, :, ne, ng, nslc) = Einv_pa * rhs_pa;
                end
            end
        end
    end

    img_sense_ap_ = zeross([N, num_echo, num_gslider, num_slice]);
    img_sense_pa_ = zeross([N, num_echo, num_gslider, num_slice]);
    for nslc = 1:num_slice
        nslc_select = Nslices*MB-nslc+1:-Nslices:1;
        for zz = 1:MB
            img_sense_ap_(:,:,:,:,nslc_select(zz)) = circshift(img_sense_ap(:,N(2)*(zz-1)+1:N(2)*zz,:,:,nslc), -shift_amount(zz), 2);
            img_sense_pa_(:,:,:,:,nslc_select(zz)) = circshift(img_sense_pa(:,N(2)*(zz-1)+1:N(2)*zz,:,:,nslc), -shift_amount(zz), 2);
        end
    end
    img_sense_ap = img_sense_ap_;
    img_sense_pa = img_sense_pa_;
    clear img_sense_ap_ img_sense_pa_
    clear sgnl_ap sgnl_pa

    elseif MB == 1
        % parfor nslc = 1:num_slice
        for nslc = 1:num_slice
            disp(['slice: ', num2str(nslc), ' / ', num2str(num_slice)])

            for xn = 1:num_read
                AC_ap = zeross([num_chan * floor(N(2) / Ry), N(2)]);
                AC_pa = zeross([num_chan * floor(N(2) / Ry), N(2)]);

                for c = 1:num_chan
                    AC_ap(1 + (c-1)*floor(N(2)/Ry) : c*floor(N(2)/Ry), :) = A_ap * diag( sens_adp(xn,:,nslc,c) );
                    AC_pa(1 + (c-1)*floor(N(2)/Ry) : c*floor(N(2)/Ry), :) = A_pa * diag( sens_adp(xn,:,nslc,c) );
                end

                [U_ap, S_ap, V_ap] = svd(AC_ap, 'econ');
                [U_pa, S_pa, V_pa] = svd(AC_pa, 'econ');

                Einv_ap = V_ap * diag(diag(S_ap) ./ (diag(S_ap).^2 + lambda_tik)) * U_ap';
                Einv_pa = V_pa * diag(diag(S_pa) ./ (diag(S_pa).^2 + lambda_tik)) * U_pa';

                for ng = 1:num_gslider
                    for ne = 1:num_echo
                        rhs_ap = sgnl_ap(xn, :, :, ne, ng, nslc);
                        rhs_ap = rhs_ap(:);
                        img_sense_ap(xn, :, ne, ng, nslc) = Einv_ap * rhs_ap;

                        rhs_pa = sgnl_pa(xn, :, :, ne, ng, nslc);
                        rhs_pa = rhs_pa(:);
                        img_sense_pa(xn, :, ne, ng, nslc) = Einv_pa * rhs_pa;
                    end
                end
            end
        end
    end


    img_low_sense_ap = zeross([N_low(1), N_low(2)*MB, num_echo, num_gslider, num_slice]);
    img_low_sense_pa = zeross([N_low(1), N_low(2)*MB, num_echo, num_gslider, num_slice]);

    num_read = N_low(1);

    if MB == 2

    shift_amount_low = (ceil(N_low(2)/(AccY_low*MB))*[0,1])*1;

    % parfor nslc = 1:num_slice
    for nslc = 1:num_slice
        disp(['slice: ', num2str(nslc), ' / ', num2str(num_slice)])

        sens_low_select = zeros(N_low(1),N_low(2),MB,num_chan);
        nslc_select = Nslices*MB-nslc+1:-Nslices:1;

        for zz = 1:MB
            % sens_low_select(:,:,zz,:) = circshift(sens_low(:,:,nslc_select(zz),:), shift_amount_low(zz), 2); % OPT 1 org
            sens_low_select(:,:,zz,:) = circshift(sens_low_adp(:,:,nslc_select(zz),:), shift_amount_low(zz), 2); % OPT 2 update: coil sens
        end

        for xn = 1:num_read
            AC_low_ap = zeross([num_chan * round(N_low(2) / Ry_low), N_low(2)*MB]);
            AC_low_pa = zeross([num_chan * round(N_low(2) / Ry_low), N_low(2)*MB]);

            for c = 1:num_chan
                AC_low_ap(1 + (c-1)*round(N_low(2)/Ry_low) : c*round(N_low(2)/Ry_low), :) = repmat(A_low_ap,1,MB) * diag( reshape(permute(sens_low_select(xn,:,:,c),[2,3,1]),[N_low(2)*MB,1]) );
                AC_low_pa(1 + (c-1)*round(N_low(2)/Ry_low) : c*round(N_low(2)/Ry_low), :) = repmat(A_low_pa,1,MB) * diag( reshape(permute(sens_low_select(xn,:,:,c),[2,3,1]),[N_low(2)*MB,1]) );
            end

            [U_low_ap, S_low_ap, V_low_ap] = svd(AC_low_ap, 'econ');
            [U_low_pa, S_low_pa, V_low_pa] = svd(AC_low_pa, 'econ');

            Einv_low_ap = V_low_ap * diag(diag(S_low_ap) ./ (diag(S_low_ap).^2 + lambda_tik)) * U_low_ap';
            Einv_low_pa = V_low_pa * diag(diag(S_low_pa) ./ (diag(S_low_pa).^2 + lambda_tik)) * U_low_pa';

            for ng = 1:num_gslider
                for ne = 1:num_echo
                    rhs_low_ap = sgnl_low_ap(xn, :, :, ne, ng, nslc);
                    rhs_low_ap = rhs_low_ap(:);
                    img_low_sense_ap(xn,:, ne, ng, nslc) = Einv_low_ap * rhs_low_ap;

                    rhs_low_pa = sgnl_low_pa(xn, :, :, ne, ng, nslc);
                    rhs_low_pa = rhs_low_pa(:);
                    img_low_sense_pa(xn,:, ne, ng, nslc) = Einv_low_pa * rhs_low_pa;
                end
            end
        end
    end

    img_low_sense_ap_ = zeross([N_low, num_echo, num_gslider, num_slice]);
    img_low_sense_pa_ = zeross([N_low, num_echo, num_gslider, num_slice]);
    for nslc = 1:num_slice
        nslc_select = Nslices*MB-nslc+1:-Nslices:1;
        for zz = 1:MB
            img_low_sense_ap_(:,:,:,:,nslc_select(zz)) = circshift(img_low_sense_ap(:,N_low(2)*(zz-1)+1:N_low(2)*zz,:,:,nslc), -shift_amount_low(zz), 2);
            img_low_sense_pa_(:,:,:,:,nslc_select(zz)) = circshift(img_low_sense_pa(:,N_low(2)*(zz-1)+1:N_low(2)*zz,:,:,nslc), -shift_amount_low(zz), 2);
        end
    end
    img_low_sense_ap = img_low_sense_ap_;
    img_low_sense_pa = img_low_sense_pa_;
    clear img_low_sense_ap_ img_low_sense_pa_
    clear sgnl_low_ap sgnl_low_pa

    elseif MB == 1
        % parfor nslc = 1:num_slice
        for nslc = 1:num_slice
            disp(['slice: ', num2str(nslc), ' / ', num2str(num_slice)])

            for xn = 1:num_read
                AC_low_ap = zeross([num_chan * floor(N_low(2) / Ry_low), N_low(2)]);
                AC_low_pa = zeross([num_chan * floor(N_low(2) / Ry_low), N_low(2)]);

                for c = 1:num_chan
                    AC_low_pa(1 + (c-1)*floor(N_low(2)/Ry_low) : c*floor(N_low(2)/Ry_low), :) = A_low_ap * diag( sens_low_adp(xn,:,nslc,c) );
                    AC_low_pa(1 + (c-1)*floor(N_low(2)/Ry_low) : c*floor(N_low(2)/Ry_low), :) = A_low_pa * diag( sens_low_adp(xn,:,nslc,c) );
                end

                [U_low_ap, S_low_ap, V_low_ap] = svd(AC_low_pa, 'econ');
                [U_low_pa, S_low_pa, V_low_pa] = svd(AC_low_pa, 'econ');

                Einv_low_ap = V_low_ap * diag(diag(S_low_ap) ./ (diag(S_low_ap).^2 + lambda_tik)) * U_low_ap';
                Einv_low_pa = V_low_pa * diag(diag(S_low_pa) ./ (diag(S_low_pa).^2 + lambda_tik)) * U_low_pa';

                for ng = 1:num_gslider
                    for ne = 1:num_echo
                        rhs_low_ap = sgnl_low_ap(xn, :, :, ne, ng, nslc);
                        rhs_low_ap = rhs_low_ap(:);
                        img_low_sense_ap(xn, :, ne, ng, nslc) = Einv_low_ap * rhs_low_ap;

                        rhs_low_pa = sgnl_low_pa(xn, :, :, ne, ng, nslc);
                        rhs_low_pa = rhs_low_pa(:);
                        img_low_sense_pa(xn, :, ne, ng, nslc) = Einv_low_pa * rhs_low_pa;
                    end
                end
            end
        end
    end

    % delete(gcp('nocreate'))

    img_sense_ap        = apodize(img_sense_ap);
    img_sense_pa        = apodize(img_sense_pa);
    img_low_sense_ap    = apodize(img_low_sense_ap);
    img_low_sense_pa    = apodize(img_low_sense_pa);
    img_sense           = cat(6, img_sense_ap, img_sense_pa); %x,y,te,slc,ap/pa
    img_low_sense       = cat(6, img_low_sense_ap, img_low_sense_pa); %x,y,te,slc,ap/pa
    img_low_sense_pad   = ifft2call(padarray(fft2call(img_low_sense),(size(img_sense)-size(img_low_sense))/2));

    delete(gcp('nocreate'))
    c = parcluster('local');
    total_cores = c.NumWorkers;
    parpool(min(num_slice,ceil(total_cores/8)))

    img_low_sense_pad_reg = zeros(size(img_low_sense_pad));
    fprintf('image registration b/w first and second echoes\n');
    for gs_index = 1:Ngslider
        parfor slc_index = 1:Nslices*MB
        % for slc_index = 1:Nslices*MB
            fprintf('gSlider idx: %d / %d, slice idx: %d / %d... \n', gs_index,Ngslider,slc_index,Nslices*MB);
            [img_fieldmap_low2high_ap(:,:,slc_index,ndwi,gs_index,:),~] = imregdeform(sq(abs(img_low_sense_pad(:,:,2,gs_index,slc_index,1))).*1e5,...
                                                                            sq(abs(img_sense(:,:,1,gs_index,slc_index,1))).*1e5,...
                                                                            NumPyramidLevels=6,...
                                                                            GridRegularization=0.6,...
                                                                            DisplayProgress=0);
            [img_fieldmap_low2high_pa(:,:,slc_index,ndwi,gs_index,:),~] = imregdeform(sq(abs(img_low_sense_pad(:,:,2,gs_index,slc_index,2))).*1e5,...
                                                                            sq(abs(img_sense(:,:,1,gs_index,slc_index,2))).*1e5,...
                                                                            NumPyramidLevels=6,...
                                                                            GridRegularization=0.6,...
                                                                            DisplayProgress=0);

            img_low_sense_pad_reg_real_ap = imwarp(sq(real(img_low_sense_pad(:,:,2,gs_index,slc_index,1))),sq(img_fieldmap_low2high_ap(:,:,slc_index,ndwi,gs_index,:)));
            img_low_sense_pad_reg_imag_ap = imwarp(sq(imag(img_low_sense_pad(:,:,2,gs_index,slc_index,1))),sq(img_fieldmap_low2high_ap(:,:,slc_index,ndwi,gs_index,:)));
            img_low_sense_pad_reg_real_pa = imwarp(sq(real(img_low_sense_pad(:,:,2,gs_index,slc_index,2))),sq(img_fieldmap_low2high_pa(:,:,slc_index,ndwi,gs_index,:)));
            img_low_sense_pad_reg_imag_pa = imwarp(sq(imag(img_low_sense_pad(:,:,2,gs_index,slc_index,2))),sq(img_fieldmap_low2high_pa(:,:,slc_index,ndwi,gs_index,:)));
            img_low_sense_pad_reg(:,:,2,gs_index,slc_index,:) = cat(6,img_low_sense_pad_reg_real_ap + 1i*img_low_sense_pad_reg_imag_ap,...
                                                                      img_low_sense_pad_reg_real_pa + 1i*img_low_sense_pad_reg_imag_pa);
        end
    end
    delete(gcp('nocreate'))

    save(['diff_',num2str(ndwi),'_img_sense.mat'], 'img_sense','ndwi','-v7.3');
    save(['diff_',num2str(ndwi),'_img_low_sense.mat'], 'img_low_sense','ndwi','-v7.3');
    save(['diff_',num2str(ndwi),'_img_low_sense_pad.mat'], 'img_low_sense_pad','ndwi','-v7.3');
    save(['diff_',num2str(ndwi),'_img_low_sense_pad_reg.mat'], 'img_low_sense_pad_reg','ndwi','-v7.3');

    clear img_sense_ap img_sense_pa img_low_sense_ap img_low_sense_pa
    toc


    %% gSlider index

    for gs_index = 1:Ngslider

        %% save as nifti

        tic
        fprintf('\n===== =============== =====\n');
        fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d',ndwi,Ndir,gs_index,Ngslider);
        fprintf('\n===== =============== =====\n');
        fprintf('saving data as nifti... ');

        save_path               = strcat(file_path,'topup/',seq_name(1:end-4),'/');

        img_recon_pad           = permute(sq(img_sense(:,:,:,gs_index,:,:)), [1,2,4,3,5]);
        img_low_recon_pad       = permute(sq(img_low_sense(:,:,:,gs_index,:,:)), [1,2,4,3,5]);
        img_low_pad_recon_pad   = permute(sq(img_low_sense_pad_reg(:,:,:,gs_index,:,:)), [1,2,4,3,5]);

        num_slc_pad             = ceil(num_slice*MB/4)*4;

        if mod(num_slc_pad-num_slice*MB,2) == 0
            img_recon_pad           = padarray(img_recon_pad, [0,0, num_slc_pad - num_slice*MB]/2);
            img_low_recon_pad       = padarray(img_low_recon_pad, [0,0, num_slc_pad - num_slice*MB]/2);
            img_low_pad_recon_pad   = padarray(img_low_pad_recon_pad, [0,0, num_slc_pad - num_slice*MB]/2);
        else
            img_recon_pad           = padarray(img_recon_pad, floor([0,0, num_slc_pad - num_slice*MB]/2),'pre');
            img_recon_pad           = padarray(img_recon_pad, ceil([0,0, num_slc_pad - num_slice*MB]/2),'post');
            img_low_recon_pad       = padarray(img_low_recon_pad, floor([0,0, num_slc_pad - num_slice*MB]/2),'pre');
            img_low_recon_pad       = padarray(img_low_recon_pad, ceil([0,0, num_slc_pad - num_slice*MB]/2),'post');
            img_low_pad_recon_pad   = padarray(img_low_pad_recon_pad, floor([0,0, num_slc_pad - num_slice*MB]/2),'pre');
            img_low_pad_recon_pad   = padarray(img_low_pad_recon_pad, ceil([0,0, num_slc_pad - num_slice*MB]/2),'post');
        end

        voxel_size          = [fov fov thickness]./[Nx,Ny,1].*1e3; % mm
        voxel_size_low      = [fov fov thickness]./[Nx_low,Ny_low,1].*1e3; % mm

        fpB0Topup           = [save_path, 'img_ap_pa_diff_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'];
        fpB0Topup_low       = [save_path, 'img_low_ap_pa_diff_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'];
        fpB0Topup_low_pad   = [save_path, 'img_low_pad_ap_pa_diff_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'];

        genNii(single(abs(sq(img_recon_pad(:,:,:,1,:)))), voxel_size, fpB0Topup)
        genNii(single(abs(sq(img_low_recon_pad(:,:,:,2,:)))), voxel_size_low, fpB0Topup_low)
        genNii(single(abs(sq(img_low_pad_recon_pad(:,:,:,2,:)))), voxel_size, fpB0Topup_low_pad)
        toc


        %% topup

        tic
        fprintf('\n===== =============== =====\n');
        fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d',ndwi,Ndir,gs_index,Ngslider);
        fprintf('\n===== =============== =====\n');
        fprintf('topup... ');

        config_path = seq_path;

        % create topup text file
        fpAcqp                  = [save_path, 'acq_param.txt'];
        fpAcqp_low              = [save_path, 'acq_low_param.txt'];
        fpAcqp_low_pad          = [save_path, 'acq_low_pad_param.txt'];

        readout_duration        = (Ny - 1)*esp/Ry;
        readout_duration_low    = (Ny_low - 1)*esp_low/Ry_low;
        readout_duration_low_pad= (Ny - 1)*esp_low/Ry_low;

        acqp0                   = [0 -1 0 readout_duration];
        acqp1                   = [0 1 0 readout_duration];
        acqp0_low               = [0 -1 0 readout_duration_low];
        acqp1_low               = [0 1 0 readout_duration_low];
        acqp0_low_pad           = [0 -1 0 readout_duration_low_pad];
        acqp1_low_pad           = [0 1 0 readout_duration_low_pad];
        acq_par                 = [acqp0; acqp1];
        acq_par_low             = [acqp0_low; acqp1_low];
        acq_par_low_pad         = [acqp0_low_pad; acqp1_low_pad];

        write_acq_params(fpAcqp, acq_par)
        write_acq_params(fpAcqp_low, acq_par_low)
        write_acq_params(fpAcqp_low_pad, acq_par_low_pad)

        fpB0Topup               = [save_path, 'img_ap_pa_diff_',num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'];
        fpB0Topup_low           = [save_path, 'img_low_ap_pa_diff_',num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'];
        fpB0Topup_low_pad       = [save_path, 'img_low_pad_ap_pa_diff_',num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'];

        fpOut                   = [save_path, 'out_ap_pa_diff_',num2str(ndwi),'_gslr_',num2str(gs_index)];
        fpOut_low               = [save_path, 'out_low_ap_pa_diff_',num2str(ndwi),'_gslr_',num2str(gs_index)];
        fpOut_low_pad           = [save_path, 'out_low_pad_ap_pa_diff_',num2str(ndwi),'_gslr_',num2str(gs_index)];

        fpIout                  = [save_path, 'epi_unwarp_diff_',num2str(ndwi),'_gslr_',num2str(gs_index)];
        fpIout_low              = [save_path, 'epi_low_unwarp_diff_',num2str(ndwi),'_gslr_',num2str(gs_index)];
        fpIout_low_pad          = [save_path, 'epi_low_pad_unwarp_diff_',num2str(ndwi),'_gslr_',num2str(gs_index)];

        fpField                 = [save_path, 'fieldmap_diff_',num2str(ndwi),'_gslr_',num2str(gs_index)];
        fpField_low             = [save_path, 'fieldmap_low_diff_',num2str(ndwi),'_gslr_',num2str(gs_index)];
        fpField_low_pad         = [save_path, 'fieldmap_low_pad_diff_',num2str(ndwi),'_gslr_',num2str(gs_index)];

        cmd                     = ['topup --imain=' fpB0Topup ' --config=', config_path 'b02b0_topup_1mm.cnf --datain=' fpAcqp ' --out=' fpOut ' --iout=' fpIout ' --fout=' fpField];
        cmd_low                 = ['topup --imain=' fpB0Topup_low ' --config=', config_path 'b02b0_topup_low.cnf --datain=' fpAcqp_low ' --out=' fpOut_low ' --iout=' fpIout_low ' --fout=' fpField_low];
        cmd_low_pad             = ['topup --imain=' fpB0Topup_low_pad ' --config=', config_path 'b02b0_topup_1mm.cnf --datain=' fpAcqp_low_pad ' --out=' fpOut_low_pad ' --iout=' fpIout_low_pad ' --fout=' fpField_low_pad];

        delete(gcp('nocreate'))
        c = parcluster('local');
        total_cores = c.NumWorkers;
        parpool(3)

       
        parfor tt = 1:3
        % for tt = 1:3
            if tt == 1
                [status,result] = system(cmd, '-echo');
            elseif tt == 2
                [status,result] = system(cmd_low, '-echo');
            else
                [status,result] = system(cmd_low_pad, '-echo');
            end
        end
        delete(gcp('nocreate'))
        toc


        %% load topup results

        tic
        fprintf('\n===== =============== =====\n');
        fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d',ndwi,Ndir,gs_index,Ngslider);
        fprintf('\n===== =============== =====\n');
        fprintf('load topup results... ');

        system(['gunzip -f ', save_path, '*.gz']);

        dt1                 = load_untouch_nii([fpField,'.nii']); % field map
        dt2                 = load_untouch_nii([fpIout,'.nii']);  % distortion corrected volumes
        dt1_low             = load_untouch_nii([fpField_low,'.nii']); % field map
        dt2_low             = load_untouch_nii([fpIout_low,'.nii']);  % distortion corrected volumes
        dt1_low_pad         = load_untouch_nii([fpField_low_pad,'.nii']); % field map
        dt2_low_pad         = load_untouch_nii([fpIout_low_pad,'.nii']);  % distortion corrected volumes

        img_fieldmap_           = flip(dt1.img,1);
        img_low_fieldmap_       = flip(dt1_low.img,1);
        img_low_pad_fieldmap_   = flip(dt1_low_pad.img,1);

        if mod(num_slc_pad-num_slice*MB,2) == 0
            pad_size_l      = (num_slc_pad - num_slice*MB)/2;
            pad_size_r      = (num_slc_pad - num_slice*MB)/2;
        else
            pad_size_l      = floor((num_slc_pad - num_slice*MB)/2);
            pad_size_r      = ceil((num_slc_pad - num_slice*MB)/2);
        end

        img_fieldmap(:,:,:,ndwi,gs_index)            = img_fieldmap_(:,:,1+pad_size_l:end-pad_size_r);
        img_low_fieldmap(:,:,:,ndwi,gs_index)        = img_low_fieldmap_(:,:,1+pad_size_l:end-pad_size_r);
        img_low_pad_fieldmap(:,:,:,ndwi,gs_index)    = img_low_pad_fieldmap_(:,:,1+pad_size_l:end-pad_size_r);
        save(['diff_',num2str(ndwi),'_img_fieldmap.mat'], 'img_fieldmap','ndwi','gs_index','-v7.3');
        save(['diff_',num2str(ndwi),'_img_low_fieldmap.mat'], 'img_low_fieldmap','ndwi','gs_index','-v7.3');
        save(['diff_',num2str(ndwi),'_img_low_pad_fieldmap.mat'], 'img_low_pad_fieldmap','ndwi','gs_index','-v7.3');
        
        dt2.img             = dt2.img(:,:,1+pad_size_l:end-pad_size_r,:);
        dt2_low.img         = dt2_low.img(:,:,1+pad_size_l:end-pad_size_r,:);
        toc


        %% apply topup on sense recon

        tic
        fprintf('\n===== =============== =====\n');
        fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d',ndwi,Ndir,gs_index,Ngslider);
        fprintf('\n===== =============== =====\n');
        fprintf('apply topup... ');

        img_sense_topup         = [];
        img_sense_low_pad_topup = [];
        img_low_sense_topup     = [];

        cd(save_path)

        genNii( real(img_recon_pad(:,:,:,1,1)), voxel_size, [save_path, 'bf_topup_blipup_re_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])
        genNii( imag(img_recon_pad(:,:,:,1,1)), voxel_size, [save_path, 'bf_topup_blipup_im_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])
        genNii( real(img_recon_pad(:,:,:,1,2)), voxel_size, [save_path, 'bf_topup_blipdown_re_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])
        genNii( imag(img_recon_pad(:,:,:,1,2)), voxel_size, [save_path, 'bf_topup_blipdown_im_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])

        genNii( real(img_low_recon_pad(:,:,:,2,1)), voxel_size_low, [save_path, 'bf_topup_low_blipup_re_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])
        genNii( imag(img_low_recon_pad(:,:,:,2,1)), voxel_size_low, [save_path, 'bf_topup_low_blipup_im_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])
        genNii( real(img_low_recon_pad(:,:,:,2,2)), voxel_size_low, [save_path, 'bf_topup_low_blipdown_re_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])
        genNii( imag(img_low_recon_pad(:,:,:,2,2)), voxel_size_low, [save_path, 'bf_topup_low_blipdown_im_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])

        genNii( real(img_low_pad_recon_pad(:,:,:,2,1)), voxel_size, [save_path, 'bf_topup_low_pad_blipup_re_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])
        genNii( imag(img_low_pad_recon_pad(:,:,:,2,1)), voxel_size, [save_path, 'bf_topup_low_pad_blipup_im_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])
        genNii( real(img_low_pad_recon_pad(:,:,:,2,2)), voxel_size, [save_path, 'bf_topup_low_pad_blipdown_re_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])
        genNii( imag(img_low_pad_recon_pad(:,:,:,2,2)), voxel_size, [save_path, 'bf_topup_low_pad_blipdown_im_', num2str(ndwi),'_gslr_',num2str(gs_index),'.nii'])

        cmd = ['applytopup --imain=bf_topup_blipup_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp, ' --inindex=1 --topup=', fpOut, ' --method=jac --out=af_topup_blipup_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        cmd = ['applytopup --imain=bf_topup_blipup_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp, ' --inindex=1 --topup=', fpOut, ' --method=jac --out=af_topup_blipup_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        cmd = ['applytopup --imain=bf_topup_blipdown_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp, ' --inindex=2 --topup=', fpOut, ' --method=jac --out=af_topup_blipdown_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        cmd = ['applytopup --imain=bf_topup_blipdown_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp, ' --inindex=2 --topup=', fpOut, ' --method=jac --out=af_topup_blipdown_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        cmd = ['applytopup --imain=bf_topup_low_pad_blipup_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp, ' --inindex=1 --topup=', fpOut_low_pad, ' --method=jac --out=af_topup_low_pad_blipup_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        cmd = ['applytopup --imain=bf_topup_low_pad_blipup_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp, ' --inindex=1 --topup=', fpOut_low_pad, ' --method=jac --out=af_topup_low_pad_blipup_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        cmd = ['applytopup --imain=bf_topup_low_pad_blipdown_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp, ' --inindex=2 --topup=', fpOut_low_pad, ' --method=jac --out=af_topup_low_pad_blipdown_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        cmd = ['applytopup --imain=bf_topup_low_pad_blipdown_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp, ' --inindex=2 --topup=', fpOut_low_pad, ' --method=jac --out=af_topup_low_pad_blipdown_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        cmd = ['applytopup --imain=bf_topup_low_blipup_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp_low, ' --inindex=1 --topup=', fpOut_low, ' --method=jac --out=af_topup_low_blipup_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        cmd = ['applytopup --imain=bf_topup_low_blipup_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp_low, ' --inindex=1 --topup=', fpOut_low, ' --method=jac --out=af_topup_low_blipup_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        cmd = ['applytopup --imain=bf_topup_low_blipdown_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp_low, ' --inindex=2 --topup=', fpOut_low, ' --method=jac --out=af_topup_low_blipdown_re_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        cmd = ['applytopup --imain=bf_topup_low_blipdown_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --datain=', fpAcqp_low, ' --inindex=2 --topup=', fpOut_low, ' --method=jac --out=af_topup_low_blipdown_im_', num2str(ndwi), '_gslr_', num2str(gs_index), ' --interp=spline'];
        [status, result] = system(cmd, '-echo');

        system(['gunzip -f ', save_path, '*.gz']);

        dt_re           = load_untouch_nii([save_path, 'af_topup_blipup_re_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        dt_im           = load_untouch_nii([save_path, 'af_topup_blipup_im_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        img_1           = dt_re.img + 1i * dt_im.img;

        dt_re           = load_untouch_nii([save_path, 'af_topup_blipdown_re_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        dt_im           = load_untouch_nii([save_path, 'af_topup_blipdown_im_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        img_2           = dt_re.img + 1i * dt_im.img;

        dt_low_pad_re   = load_untouch_nii([save_path, 'af_topup_low_pad_blipup_re_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        dt_low_pad_im   = load_untouch_nii([save_path, 'af_topup_low_pad_blipup_im_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        img_low_pad_1   = dt_low_pad_re.img + 1i * dt_low_pad_im.img;

        dt_low_pad_re   = load_untouch_nii([save_path, 'af_topup_low_pad_blipdown_re_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        dt_low_pad_im   = load_untouch_nii([save_path, 'af_topup_low_pad_blipdown_im_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        img_low_pad_2   = dt_low_pad_re.img + 1i * dt_low_pad_im.img;

        dt_low_re       = load_untouch_nii([save_path, 'af_topup_low_blipup_re_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        dt_low_im       = load_untouch_nii([save_path, 'af_topup_low_blipup_im_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        img_low_1       = dt_low_re.img + 1i * dt_low_im.img;

        dt_low_re       = load_untouch_nii([save_path, 'af_topup_low_blipdown_re_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        dt_low_im       = load_untouch_nii([save_path, 'af_topup_low_blipdown_im_', num2str(ndwi), '_gslr_', num2str(gs_index), '.nii']);
        img_low_2       = dt_low_re.img + 1i * dt_low_im.img;

        img_sense_topup         = permute(cat(4, img_1, img_2), [1,2,4,3]);
        img_sense_low_pad_topup = permute(cat(4, img_low_pad_1, img_low_pad_2), [1,2,4,3]);
        img_low_sense_topup     = permute(cat(4, img_low_1, img_low_2), [1,2,4,3]);

        img_sense_topup         = img_sense_topup(:,:,:,1+pad_size_l:end-pad_size_r);
        img_sense_low_pad_topup = img_sense_low_pad_topup(:,:,:,1+pad_size_l:end-pad_size_r);
        img_low_sense_topup     = img_low_sense_topup(:,:,:,1+pad_size_l:end-pad_size_r);
        toc


        %% image reconstruction

        for slc_index = 1:Nslices*MB

            fprintf('\n===== =============== =====\n');
            fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d, Slice: %d/%d',ndwi,Ndir,gs_index,Ngslider,slc_index,Nslices*MB);
            fprintf('\n===== =============== =====\n');

            %% Data initialization

            tic
            fprintf('\n===== =============== =====\n');
            fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d, Slice: %d/%d',ndwi,Ndir,gs_index,Ngslider,slc_index,Nslices*MB);
            fprintf('\n===== =============== =====\n');
            fprintf('data init... ');

            % sens_use = sens; % OPT 1 org
            % sens_low_use = sens_low; % OPT 1 org
            sens_use = sens_adp; % OPT 2 update: coil sens
            sens_low_use = sens_low_adp; % OPT 2 update: coil sens
            if coil_comp == 1
                Nc = num_cc;
                % sens_use = sens_cc; % OPT 1 org
                % sens_low_use = sens_low_cc; % OPT 1 org
                sens_use = sens_adp_cc; % OPT 2 update: coil sens
                sens_low_use = sens_low_adp_cc; % OPT 2 update: coil sens
            end

            mask_ap         = zeros([Nx,Ny]);
            mask_pa         = zeros([Nx,Ny]);
            mask_ap_tmp     = sq(abs(kspace_ap_tmp)>1e-10);
            mask_pa_tmp     = sq(abs(kspace_pa_tmp)>1e-10);
            mask_ap_tmp_idx = single(find(sum(mask_ap_tmp,1)~=0));
            mask_pa_tmp_idx = single(find(sum(mask_pa_tmp,1)~=0));
            if length(mask_ap_tmp_idx) > length(mask_pa_tmp_idx)
                mask_ap_tmp_idx = mask_ap_tmp_idx(1:end-1);
            end
            mask_ap(1:Nx,mask_ap_tmp_idx) = 1;
            mask_pa(1:Nx,mask_pa_tmp_idx) = 1;
            mask_ap         = logical(mask_ap);
            mask_pa         = logical(mask_pa);
            kMask           = cat(3,mask_ap,mask_pa);

            data_ap_tmp     = fft2call(permute(img_sense(:,:,1,gs_index,slc_index,1),[1,2,6,3,4,5]).*permute(sens_use(:,:,slc_index,:),[1,2,4,5,6,3])).*mask_ap;
            data_pa_tmp     = fft2call(permute(img_sense(:,:,1,gs_index,slc_index,2),[1,2,6,3,4,5]).*permute(sens_use(:,:,slc_index,:),[1,2,4,5,6,3])).*mask_pa;
            data            = cat(3,permute(data_ap_tmp,[1,2,4,3]),permute(data_pa_tmp,[1,2,4,3]));
            coil_sens       = cat(3,sens_use(:,:,slc_index,:),sens_use(:,:,slc_index,:));

            mask_low_ap         = zeros([Nx_low,Ny_low]);
            mask_low_pa         = zeros([Nx_low,Ny_low]);
            mask_low_ap_tmp     = sq(abs(kspace_low_ap_tmp)>1e-10);
            mask_low_pa_tmp     = sq(abs(kspace_low_pa_tmp)>1e-10);
            mask_low_ap_tmp_idx = single(find(sum(mask_low_ap_tmp,1)~=0));
            mask_low_pa_tmp_idx = single(find(sum(mask_low_pa_tmp,1)~=0));
            if length(mask_low_ap_tmp_idx) > length(mask_low_pa_tmp_idx)
                mask_low_ap_tmp_idx = mask_low_ap_tmp_idx(1:end-1);
            end
            mask_low_ap(:,mask_low_ap_tmp_idx) = 1;
            mask_low_pa(:,mask_low_pa_tmp_idx) = 1;
            mask_low_ap         = logical(mask_low_ap);
            mask_low_pa         = logical(mask_low_pa);
            kMask_low           = cat(3,mask_low_ap,mask_low_pa);

            data_low_ap_tmp = fft2call(permute(img_low_sense(:,:,2,gs_index,slc_index,1),[1,2,6,3,4,5]).*permute(sens_low_use(:,:,slc_index,:),[1,2,4,5,6,3])).*mask_low_ap;
            data_low_pa_tmp = fft2call(permute(img_low_sense(:,:,2,gs_index,slc_index,2),[1,2,6,3,4,5]).*permute(sens_low_use(:,:,slc_index,:),[1,2,4,5,6,3])).*mask_low_pa;
            data_low        = cat(3,permute(data_low_ap_tmp,[1,2,4,3]),permute(data_low_pa_tmp,[1,2,4,3]));
            coil_sens_low   = cat(3,sens_low_use(:,:,slc_index,:),sens_low_use(:,:,slc_index,:));
            toc


            %% SENSE without B0

            if (1)
            tic
            fprintf('\n===== =============== =====\n');
            fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d, Slice: %d/%d',ndwi,Ndir,gs_index,Ngslider,slc_index,Nslices*MB);
            fprintf('\n===== =============== =====\n');
            fprintf('SENSE without B0... ');

            Ah_1  = @(x) vect(sum(conj(coil_sens) .* ift2(kMask .* reshape(x, [Nx Ny 2 Nc])),4));
            AhA_1 = @(x) vect(sum(conj(coil_sens) .* ift2(kMask .* ft2(coil_sens .* repmat(reshape(x, [Nx Ny 2]),[1 1 1 Nc]))),4));
            Ah_low  = @(x) vect(sum(conj(coil_sens_low) .* ift2(kMask_low .* reshape(x, [Nx_low Ny_low 2 Nc])),4));
            AhA_low = @(x) vect(sum(conj(coil_sens_low) .* ift2(kMask_low .* ft2(coil_sens_low .* repmat(reshape(x, [Nx_low Ny_low 2]),[1 1 1 Nc]))),4));

            % SENSE reconstruction (least-squares)
            [recon,~,~,~] = pcg(AhA_1, Ah_1(data));
            recon = reshape(recon, [Nx Ny 2]);
            sense_recon(:,:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_sense_recon.mat'], 'sense_recon','ndwi','gs_index','slc_index','-v7.3');
            end

            [recon,~,~,~] = pcg(AhA_low, Ah_low(data_low));
            recon = reshape(recon, [Nx_low Ny_low 2]);
            sense_low_recon(:,:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_sense_low_recon.mat'], 'sense_low_recon','ndwi','gs_index','slc_index','-v7.3');
            end
            toc
            end


            %% hybrid-space SENSE for EPI

            tic
            fprintf('\n===== =============== =====\n');
            fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d, Slice: %d/%d',ndwi,Ndir,gs_index,Ngslider,slc_index,Nslices*MB);
            fprintf('\n===== =============== =====\n');
            fprintf('hybrid-SENSE without B0... ');

            t_axis_ap = [0:Ry:Ny-1] * esp/Ry;
            t_axis_pa = t_axis_ap(end:-1:1);

            % Remove Partial Fourier factor from time shift

            num_phase_enc_ap = sum(mask_ap(1,:));
            num_phase_enc_pa = sum(mask_pa(1,:));
            t_axis_ap = t_axis_ap(end-num_phase_enc_ap+1:end);
            t_axis_pa = t_axis_pa(1:num_phase_enc_pa);

            % SENSE recon results
            t1 = img_sense_topup(:,:,1,:);
            t2 = img_sense_topup(:,:,2,:);

            img_sense_low_pad_topup_pn(:,:,:,slc_index) = img_sense_low_pad_topup(:,:,:,slc_index);
            if slc_index < num_slice + 1
                img_sense_low_pad_topup_k(:,:,:,slc_index) = fft2call(img_sense_low_pad_topup(:,:,:,slc_index));
                
                img_sense_low_pad_topup_k(:,:,1,slc_index) = abs(img_sense_low_pad_topup_k(:,:,1,slc_index))...
                                                                .*exp(1i*(angle(img_sense_low_pad_topup_k(:,:,1,slc_index))-pi));
                img_sense_low_pad_topup_k(:,:,2,slc_index) = abs(img_sense_low_pad_topup_k(:,:,2,slc_index))...
                                                                .*exp(1i*(angle(img_sense_low_pad_topup_k(:,:,2,slc_index))));
                img_sense_low_pad_topup_pn(:,:,:,slc_index) = ifft2call(img_sense_low_pad_topup_k(:,:,:,slc_index));
            end

            t1_low = img_sense_low_pad_topup_pn(:,:,1,:);
            t2_low = img_sense_low_pad_topup_pn(:,:,2,:);

            tmp_k       = fft2call(img_sense_low_pad_topup_pn);
            tmp_mask    = zeros(size(tmp_k));

            crop_size   = [32,32];
            tmp_mask(Nx/2-crop_size(1)/2:Nx/2+crop_size(1)/2-1,Ny/2-crop_size(2)/2:Ny/2+crop_size(2)/2-1,:) = 1;
            tmp_ku = tmp_k.*tmp_mask;
            tmp_u = ifft2call(tmp_ku);
            t1_low = tmp_u(:,:,1,:);
            t2_low = tmp_u(:,:,2,:);

            % Estimated phase
            phs1 = exp(1i * angle(t1(:,:,:,slc_index)));
            phs2 = exp(1i * angle(t2(:,:,:,slc_index)));
            if diff_table(ndwi-1,1) == 0
                phs_diff = exp(1i * angle(t2(:,:,:,slc_index)./t1(:,:,:,slc_index)));
            else
                phs_diff = exp(1i * medfilt2(angle(t1_low(:,:,:,slc_index)./t2_low(:,:,:,slc_index)),[15,15],'symmetric'));
            end
            phs_diff(isnan(phs_diff)) = 0;
            phs_diff(isinf(phs_diff)) = 0;

            % Phase embedded coil sensitivity
            coil_sens_p             = coil_sens;
            coil_sens_p(:,:,1,:)    = coil_sens_p(:,:,1,:);
            coil_sens_p(:,:,2,:)    = coil_sens_p(:,:,2,:) .* phs_diff;

            ft_mtx      = fftc(eye(Ny),1);

            uft_mtx_ap  = ft_mtx(mask_ap(1,:),:); % matrix for (undersampled) Fourier transform
            uft_mtx_pa  = ft_mtx(mask_pa(1,:),:);

            fieldmap_1  = img_fieldmap(:,:,slc_index,ndwi,gs_index); % OPT1
            fieldmap_2  = medfilt2(img_low_pad_fieldmap(:,:,slc_index,ndwi,gs_index),[7,7]); % OPT2

            % EPI Encoding matrix for each readout line
            Enc_ap_1    = zeros(Nx, size(uft_mtx_ap,1), size(uft_mtx_ap,2));
            Enc_ap_2    = zeros(Nx, size(uft_mtx_ap,1), size(uft_mtx_ap,2));
            Enc_pa_1    = zeros(Nx, size(uft_mtx_pa,1), size(uft_mtx_pa,2));
            Enc_pa_2    = zeros(Nx, size(uft_mtx_pa,1), size(uft_mtx_pa,2));

            for k = 1:Nx
                b0_1            = fieldmap_1(k,:) * 2 * pi;
                b0_2            = fieldmap_2(k,:) * 2 * pi;
                W_ap_1          = exp(1i * t_axis_ap.' *b0_1);
                W_ap_2          = exp(1i * t_axis_ap.' *b0_2);
                W_pa_1          = exp(1i * t_axis_pa.' *b0_1);
                W_pa_2          = exp(1i * t_axis_pa.' *b0_2);
                Enc_ap_1(k,:,:) = uft_mtx_ap .* W_ap_1;
                Enc_ap_2(k,:,:) = uft_mtx_ap .* W_ap_2;
                Enc_pa_1(k,:,:) = uft_mtx_pa .* W_pa_1;
                Enc_pa_2(k,:,:) = uft_mtx_pa .* W_pa_2;
            end

           
            t_axis_low_ap = [0:Ry_low:Ny_low-1] * esp_low/Ry_low;
            t_axis_low_pa = t_axis_low_ap(end:-1:1);

            % Remove Partial Fourier factor from time shift

            num_phase_enc_low_ap = sum(mask_low_ap(1,:));
            num_phase_enc_low_pa = sum(mask_low_pa(1,:));
            t_axis_low_ap = t_axis_low_ap(end-num_phase_enc_low_ap+1:end);
            t_axis_low_pa = t_axis_low_pa(1:num_phase_enc_low_pa);

            % SENSE recon results
            t1_low = img_low_sense_topup(:,:,1,:);
            t2_low = img_low_sense_topup(:,:,2,:);

            tmp_k       = fft2call(img_low_sense_topup);
            tmp_mask    = zeros(size(tmp_k));

            crop_size   = [32,32];
            tmp_mask(Nx_low/2-crop_size(1)/2:Nx_low/2+crop_size(1)/2-1,Ny_low/2-crop_size(2)/2:Ny_low/2+crop_size(2)/2-1,:) = 1;
            tmp_ku = tmp_k.*tmp_mask;
            tmp_u = ifft2call(tmp_ku);
            t1_low = tmp_u(:,:,1,:);
            t2_low = tmp_u(:,:,2,:);

            % Estimated phase
            phs1_low = exp(1i * angle(t1_low(:,:,:,slc_index)));
            phs2_low = exp(1i * angle(t2_low(:,:,:,slc_index)));
            if diff_table(ndwi-1,1) == 0
                phs_diff = exp(1i * angle(img_low_sense_topup(:,:,2,slc_index)./img_low_sense_topup(:,:,1,slc_index)));
            else
                phs_diff = exp(1i * medfilt2(angle(t2_low(:,:,:,slc_index)./t1_low(:,:,:,slc_index)),[15,15],'symmetric'));
            end
            phs_diff(isnan(phs_diff)) = 0;
            phs_diff(isinf(phs_diff)) = 0;

            % Phase embedded coil sensitivity
            coil_sens_low_p             = coil_sens_low;
            coil_sens_low_p(:,:,1,:)    = coil_sens_low_p(:,:,1,:);
            coil_sens_low_p(:,:,2,:)    = coil_sens_low_p(:,:,2,:) .* phs_diff;

            ft_mtx_low  = fftc(eye(Ny_low),1);

            uft_mtx_low_ap  = ft_mtx_low(mask_low_ap(1,:),:); % matrix for (undersampled) Fourier transform
            uft_mtx_low_pa  = ft_mtx_low(mask_low_pa(1,:),:);

            fieldmap_low    = img_low_fieldmap(:,:,slc_index,ndwi,gs_index);

            % EPI Encoding matrix for each readout line
            Enc_low_ap = zeros(Nx_low, size(uft_mtx_low_ap,1), size(uft_mtx_low_ap,2));
            Enc_low_pa = zeros(Nx_low, size(uft_mtx_low_pa,1), size(uft_mtx_low_pa,2));

            for k = 1:Nx_low
                b0_low          = fieldmap_low(k,:) * 2 * pi;
                W_low_ap        = exp(1i * t_axis_low_ap.' *b0_low);
                W_low_pa        = exp(1i * t_axis_low_pa.' *b0_low);
                Enc_low_ap(k,:,:)   = uft_mtx_low_ap .* W_low_ap;
                Enc_low_pa(k,:,:)   = uft_mtx_low_pa .* W_low_pa;
            end
            %
            toc


            %% hybrid-space SENSE with B0

            tic
            fprintf('\n===== =============== =====\n');
            fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d, Slice: %d/%d',ndwi,Ndir,gs_index,Ngslider,slc_index,Nslices*MB);
            fprintf('\n===== =============== =====\n');
            fprintf('hybrid-SENSE with B0... ');

            A_1     = @(x) ft2_b0(coil_sens_p.*reshape(x,[Nx Ny]), Enc_ap_1, Enc_pa_1);
            Ah_1    = @(x) vect(sum(conj(coil_sens_p).*ift2_b0(x, Enc_ap_1, Enc_pa_1),[3 4]));
            AhA_1   = @(x) Ah_1(A_1(x));
            A_2     = @(x) ft2_b0(coil_sens_p.*reshape(x,[Nx Ny]), Enc_ap_2, Enc_pa_2);
            Ah_2    = @(x) vect(sum(conj(coil_sens_p).*ift2_b0(x, Enc_ap_2, Enc_pa_2),[3 4]));
            AhA_2   = @(x) Ah_2(A_2(x));
            A_low   = @(x) ft2_b0(coil_sens_low_p.*reshape(x,[Nx_low Ny_low]), Enc_low_ap, Enc_low_pa);
            Ah_low  = @(x) vect(sum(conj(coil_sens_low_p).*ift2_b0(x, Enc_low_ap, Enc_low_pa),[3 4]));
            AhA_low = @(x) Ah_low(A_low(x));

            data_m  = fftshift(ifft(ifftshift(data),[],1)) / sqrt(Nx);
            data_ap = data_m(:,mask_ap(1,:),1,:);
            data_pa = data_m(:,mask_pa(1,:),2,:);
            data_m  = cat(3, data_ap, data_pa);

            if (1)
            [recon,~,~,~] = pcg(AhA_1 , Ah_1(data_m), 1e-6, 10);
            recon = reshape(recon, [Nx Ny]);

            hybrid_sense_1(:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_hybrid_sense_1.mat'], 'hybrid_sense_1','ndwi','gs_index','slc_index','-v7.3');
            end
            end
            
            if (1)
            [recon,~,~,~] = pcg(AhA_2 , Ah_2(data_m), 1e-6, 3);
            recon = reshape(recon, [Nx Ny]);

            hybrid_sense_2(:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_hybrid_sense_2.mat'], 'hybrid_sense_2','ndwi','gs_index','slc_index','-v7.3');
            end
            end

            if (1)
            data_low_m  = fftshift(ifft(ifftshift(data_low),[],1)) / sqrt(Nx_low);
            data_low_ap = data_low_m(:,mask_low_ap(1,:),1,:);
            data_low_pa = data_low_m(:,mask_low_pa(1,:),2,:);
            data_low_m  = cat(3, data_low_ap, data_low_pa);

            [recon,~,~,~] = pcg(AhA_low , Ah_low(data_low_m), 1e-6, 3);
            recon = reshape(recon, [Nx_low Ny_low]);

            hybrid_sense_low(:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_hybrid_sense_low.mat'], 'hybrid_sense_low','ndwi','gs_index','slc_index','-v7.3');
            end
            end
            toc


            %% SENSE with B0 (Forward model for later reconstructions)

            if (1)
            tic
            fprintf('\n===== =============== =====\n');
            fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d, Slice: %d/%d',ndwi,Ndir,gs_index,Ngslider,slc_index,Nslices*MB);
            fprintf('\n===== =============== =====\n');
            fprintf('SENSE with B0... ');

            % SENSE with B0
            A_1   = @(x) ft2_b0(coil_sens.*reshape(x,[Nx Ny 2]), Enc_ap_1, Enc_pa_1);
            Ah_1  = @(x) vect(sum(conj(coil_sens).*ift2_b0(x, Enc_ap_1, Enc_pa_1),4));
            AhA_1 = @(x) Ah_1(A_1(x));
            A_2   = @(x) ft2_b0(coil_sens.*reshape(x,[Nx Ny 2]), Enc_ap_2, Enc_pa_2);
            Ah_2  = @(x) vect(sum(conj(coil_sens).*ift2_b0(x, Enc_ap_2, Enc_pa_2),4));
            AhA_2 = @(x) Ah_2(A_2(x));
            A_low   = @(x) ft2_b0(coil_sens_low.*reshape(x,[Nx_low Ny_low 2]), Enc_low_ap, Enc_low_pa);
            Ah_low  = @(x) vect(sum(conj(coil_sens_low).*ift2_b0(x, Enc_low_ap, Enc_low_pa),4));
            AhA_low = @(x) Ah_low(A_low(x));

            if (1)
            Ahd_1 = Ah_1(data_m);

            [recon,~,~,~] = pcg(AhA_1 , Ahd_1);
            recon = reshape(recon, [Nx Ny 2]);

            sense_b0_1(:,:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_sense_b0_1.mat'], 'sense_b0_1','ndwi','gs_index','slc_index','-v7.3');
            end
            end

            Ahd_2 = Ah_2(data_m);

            [recon,~,~,~] = pcg(AhA_2 , Ahd_2);
            recon = reshape(recon, [Nx Ny 2]);

            sense_b0_2(:,:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_sense_b0_2.mat'], 'sense_b0_2','ndwi','gs_index','slc_index','-v7.3');
            end

            if (1)
            Ahd_low = Ah_low(data_low_m);

            [recon,~,~,~] = pcg(AhA_low , Ahd_low);
            recon = reshape(recon, [Nx_low Ny_low 2]);

            sense_low_b0(:,:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_sense_low_b0.mat'], 'sense_low_b0','ndwi','gs_index','slc_index','-v7.3');
            end
            end
            toc
            end


            %% BUDA-LORAKS (single-channel matrix)

            if (1)
            tic
            fprintf('\n===== =============== =====\n');
            fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d, Slice: %d/%d',ndwi,Ndir,gs_index,Ngslider,slc_index,Nslices*MB);
            fprintf('\n===== =============== =====\n');
            fprintf('BUDA-LORAKS Reconstruction (Single-channel matrix)\n');

            LORAKS_type = 1; % Type 1: S-matrix (Support + Phase + parallel imaging)
            % Type 2: C matrix (Support + Parallel imaging)

            z_1 = vect(sense_b0_1(:,:,:,slc_index,ndwi,gs_index)); %initialization
            z_2 = vect(sense_b0_2(:,:,:,slc_index,ndwi,gs_index)); %initialization
            z_low = vect(sense_low_b0(:,:,:,slc_index,ndwi,gs_index)); %initialization

            % parameters for tuning
            Rad     = 5;     % kernel radius (3~5 is recommended)
            rank    = 150;   % rank for the matrix
            lambda  = 1e-3;  % lambda

            %
            tol = 1e-4;      % tolerance
            max_iter = 20;   % change this based on necessary # iterations

            [in1,in2] = meshgrid(-Rad:Rad,-Rad:Rad);
            idx = find(in1.^2+in2.^2<=Rad^2);
            patchSize = numel(idx);

            B = @(x) vect(ft2(reshape(x, [Nx Ny Nseg])));           % Fourier Transform
            Bh = @(x) vect(sum(ift2(reshape(x,[Nx, Ny, Nseg])),4));         % Adjoint
            B_low = @(x) vect(ft2(reshape(x, [Nx_low Ny_low Nseg])));           % Fourier Transform
            Bh_low = @(x) vect(sum(ift2(reshape(x,[Nx_low, Ny_low, Nseg])),4));         % Adjoint

            P_M = @(x) LORAKS_operators(x,Nx,Ny,Nseg,Rad,LORAKS_type,[]);
            Ph_M = @(x) LORAKS_operators(x,Nx,Ny,Nseg,Rad,-LORAKS_type,[]);
            P_M_low = @(x) LORAKS_operators(x,Nx_low,Ny_low,Nseg,Rad,LORAKS_type,[]);
            Ph_M_low = @(x) LORAKS_operators(x,Nx_low,Ny_low,Nseg,Rad,-LORAKS_type,[]);

            N1 = Nx; N2 = Ny; Nc2 = Nseg;
            N1_low = Nx_low; N2_low = Ny_low;

            ZD = @(x) padarray(reshape(x,[N1 N2 Nc2]),[2*Rad, 2*Rad], 'post');
            ZD_H = @(x) x(1:N1,1:N2,:,:);
            ZD_low = @(x) padarray(reshape(x,[N1_low N2_low Nc2]),[2*Rad, 2*Rad], 'post');
            ZD_H_low = @(x) x(1:N1_low,1:N2_low,:,:);

            if (1)
            fprintf('BUDA-LORAKS Reconstruction (Single-channel matrix) (1)... ');
            for iter = 1:max_iter
                z_cur = z_1;
                pz = B(z_1);
                MM = P_M(pz);

                Um = svd_left(MM);
                nmm = Um(:,rank+1:end)'; % null space
                Bhr = 0;

                if LORAKS_type == 1 % S
                    nf = size(nmm,1);
                    nmm = reshape(nmm,[nf, patchSize, 2*Nc2]);
                    nss_h = reshape(nmm(:,:,1:2:end)+1j*nmm(:,:,2:2:end),[nf, patchSize*Nc2]);
                    Nis = filtfilt(nss_h,'C',N1,N2,Nc2,Rad);
                    Nis2 = filtfilt(nss_h,'S',N1,N2,Nc2,Rad);
                    LhL = @(x) 2*Bh((ZD_H(ifft2(squeeze(sum(Nis.*repmat(fft2(ZD(B(x))),[1 1 1 Nc2]),3))))) ...
                        -(ZD_H(ifft2(squeeze(sum(Nis2.*repmat(conj(fft2(ZD(B(x)))),[1 1 1 Nc2]),3))))));
                end

                % data fitting
                M = @(x) AhA_1(x) + lambda*LhL(x);
                [z_1,~] = pcg(M, Ahd_1 + lambda*Bhr);

                t = (norm(z_cur-z_1)/norm(z_1));

                % display the status
                %             if ~rem(iter,1)
                %                 disp(['iter ' num2str(iter) ', relative change in solution: ' num2str(t)]);
                %             end

                % check for convergence
                if t < tol
                    disp('Convergence tolerance met: change in solution is small');
                    break;
                end

            end

            recon = reshape(z_1, [Nx Ny 2]);
            BUDA_1(:,:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_BUDA_1.mat'], 'BUDA_1','ndwi','gs_index','slc_index','-v7.3');
            end
            toc
            end

            if (1)
            tic
            fprintf('BUDA-LORAKS Reconstruction (Single-channel matrix) (2)... ');
            for iter = 1:max_iter
                z_cur = z_2;
                pz = B(z_2);
                MM = P_M(pz);

                Um = svd_left(MM);
                nmm = Um(:,rank+1:end)'; % null space
                Bhr = 0;

                if LORAKS_type == 1 % S
                    nf = size(nmm,1);
                    nmm = reshape(nmm,[nf, patchSize, 2*Nc2]);
                    nss_h = reshape(nmm(:,:,1:2:end)+1j*nmm(:,:,2:2:end),[nf, patchSize*Nc2]);
                    Nis = filtfilt(nss_h,'C',N1,N2,Nc2,Rad);
                    Nis2 = filtfilt(nss_h,'S',N1,N2,Nc2,Rad);
                    LhL = @(x) 2*Bh((ZD_H(ifft2(squeeze(sum(Nis.*repmat(fft2(ZD(B(x))),[1 1 1 Nc2]),3))))) ...
                        -(ZD_H(ifft2(squeeze(sum(Nis2.*repmat(conj(fft2(ZD(B(x)))),[1 1 1 Nc2]),3))))));
                end

                % data fitting
                M = @(x) AhA_2(x) + lambda*LhL(x);
                [z_2,~] = pcg(M, Ahd_2 + lambda*Bhr);

                t = (norm(z_cur-z_2)/norm(z_2));

                % display the status
                %             if ~rem(iter,1)
                %                 disp(['iter ' num2str(iter) ', relative change in solution: ' num2str(t)]);
                %             end

                % check for convergence
                if t < tol
                    disp('Convergence tolerance met: change in solution is small');
                    break;
                end

            end

            recon = reshape(z_2, [Nx Ny 2]);
            BUDA_2(:,:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_BUDA_2.mat'], 'BUDA_2','ndwi','gs_index','slc_index','-v7.3');
            end
            toc
            end

            if (1)
            tic
            fprintf('BUDA-LORAKS Reconstruction (Single-channel matrix) (low)... ');
            for iter = 1:max_iter
                z_cur = z_low;
                pz = B_low(z_low);
                MM = P_M_low(pz);

                Um = svd_left(MM);
                nmm = Um(:,rank+1:end)'; % null space
                Bhr = 0;

                if LORAKS_type == 1 % S
                    nf = size(nmm,1);
                    nmm = reshape(nmm,[nf, patchSize, 2*Nc2]);
                    nss_h = reshape(nmm(:,:,1:2:end)+1j*nmm(:,:,2:2:end),[nf, patchSize*Nc2]);
                    Nis = filtfilt(nss_h,'C',N1_low,N2_low,Nc2,Rad);
                    Nis2 = filtfilt(nss_h,'S',N1_low,N2_low,Nc2,Rad);
                    LhL = @(x) 2*Bh_low((ZD_H_low(ifft2(squeeze(sum(Nis.*repmat(fft2(ZD_low(B_low(x))),[1 1 1 Nc2]),3))))) ...
                        -(ZD_H_low(ifft2(squeeze(sum(Nis2.*repmat(conj(fft2(ZD_low(B_low(x)))),[1 1 1 Nc2]),3))))));
                end

                % data fitting
                M = @(x) AhA_low(x) + lambda*LhL(x);
                [z_low,~] = pcg(M, Ahd_low + lambda*Bhr);

                t = (norm(z_cur-z_low)/norm(z_low));

                % display the status
                %             if ~rem(iter,1)
                %                 disp(['iter ' num2str(iter) ', relative change in solution: ' num2str(t)]);
                %             end

                % check for convergence
                if t < tol
                    disp('Convergence tolerance met: change in solution is small');
                    break;
                end

            end

            recon = reshape(z_low, [Nx_low Ny_low 2]);
            BUDA_low(:,:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_BUDA_low.mat'], 'BUDA_low','ndwi','gs_index','slc_index','-v7.3');
            end
            toc
            end
            end


            %% BUDA-LORAKS (multi-channel matrix)

            if (0)
            tic
            fprintf('\n===== =============== =====\n');
            fprintf('Diff Dir: %d/%d, gSlider idx: %d/%d, Slice: %d/%d',ndwi,Ndir,gs_index,Ngslider,slc_index,Nslices*MB);
            fprintf('\n===== =============== =====\n');
            fprintf('BUDA-LORAKS Reconstruction (Multi-channel matrix)\n');

            LORAKS_type = 1; % Type 1: S-matrix (Support + Phase + parallel imaging)
            % Type 2: C matrix (Support + Parallel imaging)

            % load buda_recon
            z_1 = vect(BUDA_1(:,:,:,slc_index,ndwi,gs_index));         % initialization
            z_2 = vect(BUDA_2(:,:,:,slc_index,ndwi,gs_index));         % initialization
            z_low = vect(BUDA_low(:,:,:,slc_index,ndwi,gs_index));         % initialization
            % z = vect(sense_b0);   % initialization

            % parameters for tuning

            Rad = 3;                  % kernel radius (3~4 is recommended)
            rank = 100;             % rank for the matrix
            lambda = 1e-3;          % lambda

            %
            tol = 1e-4;
            max_iter = 5;  % change this based on necessary # iterations

            [in1,in2] = meshgrid(-Rad:Rad,-Rad:Rad);
            idx = find(in1.^2+in2.^2<=Rad^2);
            patchSize = numel(idx);

            B = @(x) vect(ft2(coil_sens.*reshape(x, [Nx Ny Nseg])));           % SENSE encoding
            Bh = @(x) vect(sum(conj(coil_sens).*ift2(reshape(x,[Nx, Ny, Nseg Nc])),4));         % Adjoint of SENSE encoding
            B_low = @(x) vect(ft2(coil_sens_low.*reshape(x, [Nx_low Ny_low Nseg])));           % SENSE encoding
            Bh_low = @(x) vect(sum(conj(coil_sens_low).*ift2(reshape(x,[Nx_low, Ny_low, Nseg Nc])),4));         % Adjoint of SENSE encoding

            P_M = @(x) LORAKS_operators(x,Nx,Ny,Nseg*Nc,Rad,LORAKS_type,[]);
            Ph_M = @(x) LORAKS_operators(x,Nx,Ny,Nseg*Nc,Rad,-LORAKS_type,[]);
            P_M_low = @(x) LORAKS_operators(x,Nx_low,Ny_low,Nseg*Nc,Rad,LORAKS_type,[]);
            Ph_M_low = @(x) LORAKS_operators(x,Nx_low,Ny_low,Nseg*Nc,Rad,-LORAKS_type,[]);

            N1 = Nx; N2 = Ny; Nc2 = Nc*Nseg;
            N1_low = Nx_low; N2_low = Ny_low;

            ZD = @(x) padarray(reshape(x,[N1 N2 Nc2]),[2*Rad, 2*Rad], 'post');
            ZD_H = @(x) x(1:N1,1:N2,:,:);
            ZD_low = @(x) padarray(reshape(x,[N1_low N2_low Nc2]),[2*Rad, 2*Rad], 'post');
            ZD_H_low = @(x) x(1:N1_low,1:N2_low,:,:);

            if (0)
            fprintf('BUDA-LORAKS Reconstruction (Multi-channel matrix) (1)... ');
            for iter = 1:max_iter
                z_cur = z_1;
                pz = B(z_1);
                MM = P_M(pz);

                Um = svd_left(MM);
                nmm = Um(:,rank+1:end)'; % null space
                Bhr = 0;

                if LORAKS_type == 1 % S
                    nf = size(nmm,1);
                    nmm = reshape(nmm,[nf, patchSize, 2*Nc2]);
                    nss_h = reshape(nmm(:,:,1:2:end)+1j*nmm(:,:,2:2:end),[nf, patchSize*Nc2]);
                    Nis = filtfilt(nss_h,'C',N1,N2,Nc2,Rad);
                    Nis2 = filtfilt(nss_h,'S',N1,N2,Nc2,Rad);
                    LhL = @(x) 2*Bh((ZD_H(ifft2(squeeze(sum(Nis.*repmat(fft2(ZD(B(x))),[1 1 1 Nc2]),3))))) ...
                        -(ZD_H(ifft2(squeeze(sum(Nis2.*repmat(conj(fft2(ZD(B(x)))),[1 1 1 Nc2]),3))))));
                end

                % data fitting
                M = @(x) AhA_1(x) + lambda*LhL(x);
                [z_1,~] = pcg(M, Ahd_1 + lambda*Bhr);

                t = (norm(z_cur-z_1)/norm(z_1));

                % display the status
                %             if ~rem(iter,1)
                %                 disp(['iter ' num2str(iter) ', relative change in solution: ' num2str(t)]);
                %             end

                % check for convergence
                if t < tol
                    disp('Convergence tolerance met: change in solution is small');
                    break;
                end
            end

            recon = reshape(z_1, [Nx Ny 2]);
            BUDA_LORAKS_1(:,:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_BUDA_LORAKS_1.mat'], 'BUDA_LORAKS_1','ndwi','gs_index','slc_index','-v7.3');
            end
            toc
            end

            if (1)
            tic
            fprintf('BUDA-LORAKS Reconstruction (Multi-channel matrix) (2)... ');
            for iter = 1:max_iter
                z_cur = z_2;
                pz = B(z_2);
                MM = P_M(pz);

                Um = svd_left(MM);
                nmm = Um(:,rank+1:end)'; % null space
                Bhr = 0;

                if LORAKS_type == 1 % S
                    nf = size(nmm,1);
                    nmm = reshape(nmm,[nf, patchSize, 2*Nc2]);
                    nss_h = reshape(nmm(:,:,1:2:end)+1j*nmm(:,:,2:2:end),[nf, patchSize*Nc2]);
                    Nis = filtfilt(nss_h,'C',N1,N2,Nc2,Rad);
                    Nis2 = filtfilt(nss_h,'S',N1,N2,Nc2,Rad);
                    LhL = @(x) 2*Bh((ZD_H(ifft2(squeeze(sum(Nis.*repmat(fft2(ZD(B(x))),[1 1 1 Nc2]),3))))) ...
                        -(ZD_H(ifft2(squeeze(sum(Nis2.*repmat(conj(fft2(ZD(B(x)))),[1 1 1 Nc2]),3))))));
                end

                % data fitting
                M = @(x) AhA_2(x) + lambda*LhL(x);
                [z_2,~] = pcg(M, Ahd_2 + lambda*Bhr);

                t = (norm(z_cur-z_2)/norm(z_2));

                % display the status
                %             if ~rem(iter,1)
                %                 disp(['iter ' num2str(iter) ', relative change in solution: ' num2str(t)]);
                %             end

                % check for convergence
                if t < tol
                    disp('Convergence tolerance met: change in solution is small');
                    break;
                end
            end

            recon = reshape(z_2, [Nx Ny 2]);
            BUDA_LORAKS_2(:,:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_BUDA_LORAKS_2.mat'], 'BUDA_LORAKS_2','ndwi','gs_index','slc_index','-v7.3');
            end
            toc
            end

            tic
            fprintf('BUDA-LORAKS Reconstruction (Multi-channel matrix) (low)... ');
            for iter = 1:max_iter
                z_cur = z_low;
                pz = B_low(z_low);
                MM = P_M_low(pz);

                Um = svd_left(MM);
                nmm = Um(:,rank+1:end)'; % null space
                Bhr = 0;

                if LORAKS_type == 1 % S
                    nf = size(nmm,1);
                    nmm = reshape(nmm,[nf, patchSize, 2*Nc2]);
                    nss_h = reshape(nmm(:,:,1:2:end)+1j*nmm(:,:,2:2:end),[nf, patchSize*Nc2]);
                    Nis = filtfilt(nss_h,'C',N1_low,N2_low,Nc2,Rad);
                    Nis2 = filtfilt(nss_h,'S',N1_low,N2_low,Nc2,Rad);
                    LhL = @(x) 2*Bh_low((ZD_H_low(ifft2(squeeze(sum(Nis.*repmat(fft2(ZD_low(B_low(x))),[1 1 1 Nc2]),3))))) ...
                        -(ZD_H_low(ifft2(squeeze(sum(Nis2.*repmat(conj(fft2(ZD_low(B_low(x)))),[1 1 1 Nc2]),3))))));
                end

                % data fitting
                M = @(x) AhA_low(x) + lambda*LhL(x);
                [z_low,~] = pcg(M, Ahd_low + lambda*Bhr);

                t = (norm(z_cur-z_low)/norm(z_low));

                % display the status
                %             if ~rem(iter,1)
                %                 disp(['iter ' num2str(iter) ', relative change in solution: ' num2str(t)]);
                %             end

                % check for convergence
                if t < tol
                    disp('Convergence tolerance met: change in solution is small');
                    break;
                end
            end

            recon = reshape(z_low, [Nx_low Ny_low 2]);
            BUDA_LORAKS_low(:,:,:,slc_index,ndwi,gs_index) = recon;
            if slc_index == Nslices*MB
                save(['diff_',num2str(ndwi),'_BUDA_LORAKS_low.mat'], 'BUDA_LORAKS_low','ndwi','gs_index','slc_index','-v7.3');
            end
            toc
            end

        end
    end
end