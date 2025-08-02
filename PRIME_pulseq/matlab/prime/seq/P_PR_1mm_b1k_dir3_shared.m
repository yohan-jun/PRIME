%% PRIME

clc;  close all;  clear;

diff_set        =   1;
seq_file        =   'P_PR_1mm_b1k_dir3_shared.seq';
save_seq_file   =   1;
seq_pns         =   0;
seq_plot        =   0;
seq_report      =   0;
check_freq      =   0;
test_check      =   0; % for debugging

bay             =   4; % Bay 4/8


%% Set system limits

B0field     =   2.89;
lims        =   mr.opts('MaxGrad',78,'GradUnit','mT/m',...
                        'MaxSlew',108,'SlewUnit','T/m/s',...
                        'rfRingdownTime',10e-6,'rfDeadtime',100e-6,'B0',B0field); % for readout
lims1       =   mr.opts('MaxGrad',68,'GradUnit','mT/m',...
                        'MaxSlew',100,'SlewUnit','T/m/s',...
                        'rfRingdownTime',10e-6,'rfDeadtime',100e-6,'B0',B0field); % for gz fat
lims2       =   mr.opts('MaxGrad',68,'GradUnit','mT/m',...
                        'MaxSlew',60,'SlewUnit','T/m/s',...
                        'rfRingdownTime',10e-6,'rfDeadtime',100e-6,'B0',B0field); % for diffusion gradients


%% Imaging parameter

table               =   [0 0 1; 0 1 0; 1 0 0;]; % test one axis
table               =   [0 0 0; 0 0 0; table]; % dummy + b0 + dwi v2

% load('diff_table.mat');
% diff_table          =   diff_table(15*(diff_set-1)+1:15*diff_set,:); % split diff directions
% table               =   [0 0 0; diff_table(:,2:4)]; % dummy + dwi (b0 included)
diffusion_count     =   size(table,1);

fov                 =   220e-3;
Nx                  =   220;
Nx_org              =   220;            % for gx pre-phasing calculation
Nx_low              =   110;
Nx_low_org          =   110;
Ny                  =   Nx;             % Define FOV and resolution
Ny_low              =   110;            % Define FOV and resolution
thickness           =   5*220/220*1e-3;           % slice thinckness
Nslices             =   26;
MB                  =   2;              % multi-band
CAIPI               =   1;              % CAIPI (0/1)
caipi_factor        =   2;
gamma               =   425800;

bval                =   [1000];
if length(bval) == 1
    bFactor         =   bval.*ones(1,size(table,1)); % s/mm^2
else
    bFactor         =   ones(1,size(table,1)); % s/mm^2
end
b0_idx              =   find(sum(table,2)==0);
bFactor(b0_idx)     =   0;

% for bb = 1:length(bval)
%     b_idx           =   find(diff_table(:,1)==bb); % v7
%     bFactor(b_idx+1)=   bval(bb);
% end

% In siemens's sequence, only b0 sequence uses crushers, so in our
% sequence, we introduce a factor to keep the echo time the same for b0/b1000
crusher_switch      =   0;
crusher_d           =   0.95e-3;
dur_rephase         =   7.6e-04;
recon               =   false;          % plot the traj for check
TE                  =   75e-3;
TE_low              =   95e-3;
TR                  =   3500e-3;

RSegment            =   2;              % multishot factor
R                   =   1.5;            % In-plane accerelation factor
R_low               =   1.0;            % In-plane accerelation factor

Echotimeshift       =   1;              % for multishot
buda                =   1;

ro_os               =   1.00;           % oversampling factor (in contrast to the product sequence we don't really need it)
ro_os_low           =   1.00;           % oversampling factor (in contrast to the product sequence we don't really need it)
readoutTime         =   9.0e-4;         % this controls the readout bandwidth
readoutTime_low     =   3.5e-4;         % this controls the readout bandwidth
cal_readoutTime_low =   1;              % calculate readoutTime_low based on readoutTime, R, R_low

partFourierFactor           =  6/8;     % partial Fourier factor: 1: full sampling 0: start with ky=0
partFourierFactor_fe        =  1;       % frequency-encoding partial Fourier factor
partFourierFactor_low_fe    =  1;       % frequency-encoding partial Fourier factor
Nx                  =   Nx * partFourierFactor_fe;
Nx_low              =   floor(Nx_low * partFourierFactor_low_fe);

tRFex               =   3e-3;           % sec
tRFref              =   3e-3;           % sec
sat_ppm             =   -3.45;
dur_rewinder        =   1.0e-03;

deltak              =   1/fov;
deltaky             =   RSegment*R*deltak;
deltaky_low         =   RSegment*R_low*deltak;
kWidth              =   Nx*deltak;
kWidth_org          =   Nx_org*deltak;
kWidth_low          =   Nx_low*deltak;
kWidth_low_org      =   Nx_low_org*deltak;
kz                  =   1/(Nslices/MB*thickness*caipi_factor);

rf_type             =   1; % 1 (gSlider), 2 (sinc)


if rf_type == 1
    %% Load RF
    load('rf_MB1_1mm_s26.mat') % gz_90, gz_180
    load('seq_blocks_MB2_1mm_s26_PR.mat') % rf, gz blocks

    rf_order = [1,2,3,4,5];
end


%% Interleaved slice factor

for islice_1 = 1:Nslices
    freqOffset_factor(islice_1)=islice_1-1-(Nslices-1)/2;
end

Nslices         = Nslices/MB;
% thick_factor    = freqOffset_factor(1:Nslices)*thickness;

% interleaved
slic_index      = [2:2:Nslices 1:2:Nslices];
if rf_type == 2 && MB == 2
    interleaved_freqOffset_factor = freqOffset_factor(Nslices/2+1:Nslices*3/2);
    interleaved_freqOffset_factor = interleaved_freqOffset_factor(slic_index);
else
    interleaved_freqOffset_factor = freqOffset_factor(slic_index);
end
thick_factor    = interleaved_freqOffset_factor*thickness;


if rf_type == 1 && MB == 2
    %% Excitation pulse p1, phase 1
    length_rf_90                            = length(rf(1).signal);
    if length_rf_90 < length(g_90_versed)
        gz_test                             = g_90_versed(1:length_rf_90);
    else
        gz_test                             = zeros(length_rf_90,1);
        gz_test(1:length(g_90_versed),:)    = g_90_versed;
    end
    gz_test(end)                            = 0;
    gz_test                                 = gamma*gz_test;

    for ii = 1:length_rf_90
        if (ii == 1)
            gz_ex_integral2(ii) = gz_test(end);
        else
            gz_ex_integral2(ii) = gz_ex_integral2(ii-1) + gz_test(end-ii+1);
        end
    end

    gz_ex_integral2 = gz_ex_integral2*2e-06;
    center_phase    = (gz_ex_integral2(end))/2.0;
    gz_ex_integral2 = -(gz_ex_integral2(end:-1:1) - center_phase);
    % now we need to pre-modulate each slice based on their position

    for islice_2 = 1:Nslices
        theta1(:,islice_2)  = gz_ex_integral2'.*thick_factor(islice_2);
        p1(:, islice_2)     = exp(1j*2*pi*theta1(:,islice_2));
    end


    %% Refocusing pulse p2, phase 2

    length_rf_180                               = length(rf180.signal);
    if length_rf_180 < length(g_180_versed)
        gz180_test                              = g_180_versed(1:length_rf_180);
    else
        gz180_test                              = zeros(length_rf_180,1);
        gz180_test(1:length(g_180_versed),:)    = g_180_versed;
    end
    gz180_test(end)                             = 0;
    gz180_test                                  = gamma*gz180_test;

    for ii=1:length_rf_180
        if (ii==1)
            gz_integral2(ii) = gz180_test(end);
        else
            gz_integral2(ii) = gz_integral2(ii-1) + gz180_test(end-ii+1);
        end
    end

    gz_integral2 = gz_integral2*2e-06;
    center_phase = (gz_integral2(end))/2.0;
    gz_integral2 = -(gz_integral2(end:-1:1) - center_phase);
    % now we need to calculate for each slice

    for islice_3 = 1:Nslices
        theta2(:,islice_3)  = gz_integral2'.*thick_factor(islice_3);
        p2(:, islice_3)     = exp(1j*2*pi*theta2(:,islice_3));
    end
end


%% Define Seq

seq             =   mr.Sequence(lims);      % Create a new sequence object
% rep 1 works as the dummy scan, and I also turn off the PE gradients to collect the reference
for rep = 1:diffusion_count
    if rep == 1
        pe_enable=0;           % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
    else
        pe_enable=1;
    end
    if test_check == 1
        Nmulti_end = 1;
    else
        Nmulti_end = RSegment;
    end
    for Nmulti = 1:Nmulti_end
        % fat sat
        sat_freq            =   sat_ppm*1e-6*lims.B0*lims.gamma;
        rf_fs               =   mr.makeGaussPulse_QL(110*pi/180,'system',lims1,'Duration',8e-3,...
                                         'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
        gz_fs               =   mr.makeTrapezoid('z',lims1,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

        if rf_type == 2
            % Create 90 degree slice selection pulse and gradient
            [rf, gz, gzReph]    =   mr.makeSincPulse_QL_1(pi/2,'system',lims,'Duration',tRFex,...
                'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);

            % Create 180 degree slice refocusing pulse and gradients
            [rf180, gz180]      =   mr.makeSincPulse(pi,'system',lims,'Duration',tRFref,...
                'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',pi/2,'use','refocusing');

            if MB == 2
                Bandsep = [0 Nslices]- Nslices/2;
                nBand = 0;
                rf_90 = mrz.makeMBPulse_less(rf.signal,'n_bands',nBand,'timeBwProduct',4,'band_sep',Bandsep);
                rf.signal=rf_90;

                rf_180 = mrz.makeMBPulse_less(rf180.signal,'n_bands',nBand,'timeBwProduct',4,'band_sep',Bandsep);
                rf180.signal=rf_180;
            end
        end

        % Spoiler in front of the sequence
        spoiler_amp = 3*8*42.58*10e2;
        est_rise = 500e-6;
        est_flat = 2500e-6;

        gp_r = mr.makeTrapezoid('x','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gp_p = mr.makeTrapezoid('y','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gp_s = mr.makeTrapezoid('z','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);

        gn_r = mr.makeTrapezoid('x','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gn_p = mr.makeTrapezoid('y','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gn_s = mr.makeTrapezoid('z','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);

        if test_check == 1
            igSlider_end = 1;
        else
            igSlider_end = length(rf);
        end
        for igSlider = 1:igSlider_end

            % first gz180 crusher
            gz180_crusher_1     =   mr.makeTrapezoid('z',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d);
            gz180_crusher_2     =   mr.makeTrapezoid('y',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d);
            gz180_crusher_3     =   mr.makeTrapezoid('x',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d);

            % second gz180 crusher
            gz180_c1            =   mr.makeTrapezoid('z',lims,'Amplitude',0*19.05*42.58*10e2,'Duration',crusher_d);
            gz180_c2            =   mr.makeTrapezoid('y',lims,'Amplitude',-1*19.05*42.58*10e2,'Duration',crusher_d);
            gz180_c3            =   mr.makeTrapezoid('x',lims,'Amplitude',1*19.05*42.58*10e2,'Duration',crusher_d);

            % define the output trigger to play out with every slice excitatuion
            trig                =   mr.makeDigitalOutputPulse('osc0','duration', 100e-6); % possible channels: 'osc0','osc1','ext1'

            % Phase blip in shortest possible time
            blip_dur            =   ceil(2*sqrt(deltaky/lims.maxSlew)/10e-6/2)*10e-6*2; % round-up the duration to 2x the gradient raster time

            % the split code below fails if this really makes a trpezoid instead of a triangle...
            gy                  =   mr.makeTrapezoid('y',lims,'Area',-deltaky,'Duration',blip_dur); % use negative blips to save one k-space line on our way towards the k-space center

            % readout gradient is a truncated trapezoid with dead times at the beginnig
            % and at the end each equal to a half of blip_dur
            % the area between the blips should be defined by kWidth
            % we do a two-step calculation: we first increase the area assuming maximum
            % slewrate and then scale down the amlitude to fix the area
            extra_area          =   blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
            gx                  =   mr.makeTrapezoid('x',lims,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
            actual_area         =   gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
            gx.amplitude        =   gx.amplitude/actual_area*kWidth;
            gx.area             =   gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
            gx.flatArea         =   gx.amplitude*gx.flatTime;

            extra_area_org      =   blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
            gx_org              =   mr.makeTrapezoid('x',lims,'Area',kWidth_org+extra_area_org);
            actual_area_org     =   gx_org.area-gx_org.amplitude/gx_org.riseTime*blip_dur/2*blip_dur/2/2-gx_org.amplitude/gx_org.fallTime*blip_dur/2*blip_dur/2/2;
            gx_org.amplitude    =   gx_org.amplitude/actual_area_org*kWidth_org;
            gx_org.area         =   gx_org.amplitude*(gx_org.flatTime + gx_org.riseTime/2 + gx_org.fallTime/2);
            gx_org.flatArea     =   gx_org.amplitude*gx_org.flatTime;

            % calculate ADC
            % we use ramp sampling, so we have to calculate the dwell time and the
            % number of samples, which are will be qite different from Nx and
            % readoutTime/Nx, respectively.
            adcDwellNyquist     =   deltak/gx.amplitude/ro_os;

            % round-down dwell time to 100 ns
            adcDwell            =   floor(adcDwellNyquist*1e7)*1e-7;
            adcSamples          =   floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4

            % MZ: no idea, whether ceil,round or floor is better for the adcSamples...
            adc                 =   mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);

            % realign the ADC with respect to the gradient
            time_to_center      =   adc.dwell*((adcSamples-1)/2+0.5); % Siemens samples in the center of the dwell period
            adc.delay           =   round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
            % this rounding actually makes the sampling points on odd and even readouts
            % to appear misalligned. However, on the real hardware this misalignment is
            % much stronger anyways due to the grdient delays
            adc.id              =   1;

            % FOV positioning requires alignment to grad. raster... -> TODO

            % split the blip into two halves and produce a combined synthetic gradient
            gy_parts                    =   mr.splitGradientAt(gy, blip_dur/2, lims);
            [gy_blipup, gy_blipdown,~]  =   mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
            gy_blipdownup               =   mr.addGradients({gy_blipdown, gy_blipup}, lims);

            % pe_enable support
            gy_blipup.waveform      =   gy_blipup.waveform*pe_enable;
            gy_blipdown.waveform    =   gy_blipdown.waveform*pe_enable;
            gy_blipdownup.waveform  =   gy_blipdownup.waveform*pe_enable;

            % phase encoding and partial Fourier
            Ny_pre      =   round((partFourierFactor-1/2)*Ny-1);  % PE steps prior to ky=0, excluding the central line
            Ny_pre      =   round(Ny_pre/RSegment/R);
            Ny_post     =   round(Ny/2+1); % PE lines after the k-space center including the central line
            Ny_post     =   round(Ny_post/RSegment/R);
            Ny_meas     =   Ny_pre+Ny_post;

            % pre-phasing gradients
            gxPre           = mr.makeTrapezoid('x',lims,'Area',-gx_org.area/2);
            gyPre           = mr.makeTrapezoid('y',lims,'Area',(Ny_pre*deltaky-(Nmulti-1)*R*deltak));
            [gxPre,gyPre]   = mr.align('right',gxPre,'left',gyPre);
            % relax the PE prepahser to reduce stimulation
            gyPre           = mr.makeTrapezoid('y',lims,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre));
            gyPre.amplitude = gyPre.amplitude*pe_enable;

            % post gradient to go back to k-space center
            gyPre_post_area = -1*(Ny_post-0.5)*deltaky;

            if(mod(Nmulti,2) ~= 0)
                gyPre_post_area = (Ny_post-1)*deltaky;
            end

            gxPre_post_area_nav = -1*(0.5*gx.flatArea + gx.amplitude*gx.fallTime/2); % partFourierFactor should be larger than 0.5
            gxPre_post_area     = -1*(abs(gxPre_post_area_nav*2) - (0.5*gx_org.flatArea + gx_org.amplitude*gx_org.fallTime/2)); % partFourierFactor should be larger than 0.5
            gxPre_post_area_nav = -1*(abs(gxPre_post_area_nav*2) - (0.5*gx_org.flatArea + gx_org.amplitude*gx_org.fallTime/2)); % partFourierFactor should be larger than 0.5

            if(mod(Ny_meas,2) == 0)
                gxPre_post_area     = 0.5*gx_org.flatArea + gx_org.amplitude*gx_org.riseTime/2;
            end

            gxPre_post      = mr.makeTrapezoid('x',lims,'Area', gxPre_post_area, 'Duration', 0.0012);
            gyPre_post      = mr.makeTrapezoid('y',lims,'Area', gyPre_post_area, 'Duration', 0.0012);
            gzPre_post      = mr.makeTrapezoid('z',lims,'Area', -kz*(caipi_factor-1)/2, 'Duration', 0.0012);

            if(mod(Ny_meas,2) ~= 0)
                gzPre_post.amplitude = -gzPre_post.amplitude;
            end

            [gxPre_post,gyPre_post,gzPre_post]=mr.align('right',gxPre_post,'left',gyPre_post,gzPre_post);
            gyPre_post = mr.makeTrapezoid('y',lims,'Area',gyPre_post.area,'Duration',mr.calcDuration(gxPre_post,gyPre_post,gzPre_post));
            gzPre_post = mr.makeTrapezoid('z',lims,'Area',gzPre_post.area,'Duration',mr.calcDuration(gxPre_post,gyPre_post,gzPre_post));

            gyPre_post.amplitude = gyPre_post.amplitude*pe_enable;
            gzPre_post.amplitude = gzPre_post.amplitude*pe_enable;

            if (buda)
                if (mod(Nmulti,2) == 0)
                    gyPre           =   mr.scaleGrad(gyPre,-1);
                    gy_blipdownup   =   mr.scaleGrad(gy_blipdownup,-1);
                    gy_blipdown     =   mr.scaleGrad(gy_blipdown,-1);
                    gy_blipup       =   mr.scaleGrad(gy_blipup,-1);
                end
            end


            %% 2nd echo

            blip_dur_low        =   ceil(2*sqrt(deltaky_low/lims.maxSlew)/10e-6/2)*10e-6*2;
            gy_low              =   mr.makeTrapezoid('y',lims,'Area',-deltaky_low,'Duration',blip_dur_low);

            if cal_readoutTime_low == 1
                readoutTime_low = R_low/R*(readoutTime+blip_dur)-blip_dur_low;
            end

            extra_area_low      =   blip_dur_low/2*blip_dur_low/2*lims.maxSlew; % check unit!;
            gx_low              =   mr.makeTrapezoid('x',lims,'Area',kWidth_low+extra_area_low,'duration',readoutTime_low+blip_dur_low);
            actual_area_low     =   gx_low.area-gx_low.amplitude/gx_low.riseTime*blip_dur_low/2*blip_dur_low/2/2-gx_low.amplitude/gx_low.fallTime*blip_dur_low/2*blip_dur_low/2/2;
            gx_low.amplitude    =   gx_low.amplitude/actual_area_low*kWidth_low;
            gx_low.area         =   gx_low.amplitude*(gx_low.flatTime + gx_low.riseTime/2 + gx_low.fallTime/2);
            gx_low.flatArea     =   gx_low.amplitude*gx_low.flatTime;

            extra_area_low_org      =   blip_dur_low/2*blip_dur_low/2*lims.maxSlew; % check unit!;
            gx_low_org              =   mr.makeTrapezoid('x',lims,'Area',kWidth_low_org+extra_area_low_org,'duration',readoutTime_low+blip_dur_low);
            actual_area_low_org     =   gx_low_org.area-gx_low_org.amplitude/gx_low_org.riseTime*blip_dur_low/2*blip_dur_low/2/2-gx_low_org.amplitude/gx_low_org.fallTime*blip_dur_low/2*blip_dur_low/2/2;
            gx_low_org.amplitude    =   gx_low_org.amplitude/actual_area_low_org*kWidth_low_org;
            gx_low_org.area         =   gx_low_org.amplitude*(gx_low_org.flatTime + gx_low_org.riseTime/2 + gx_low_org.fallTime/2);
            gx_low_org.flatArea     =   gx_low_org.amplitude*gx_low_org.flatTime;

            adcDwellNyquist_low =   deltak/gx_low.amplitude/ro_os_low;
            % round-down dwell time to 100 ns
            adcDwell_low        =   floor(adcDwellNyquist_low*1e7)*1e-7;
            adcSamples_low      =   floor(readoutTime_low/adcDwell_low/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
            % MZ: no idea, whether ceil,round or floor is better for the adcSamples...
            adc_low             =   mr.makeAdc(adcSamples_low,'Dwell',adcDwell_low,'Delay',blip_dur_low/2);
            % realign the ADC with respect to the gradient
            time_to_center_low  =   adc_low.dwell*((adcSamples_low-1)/2+0.5); % Siemens samples in the center of the dwell period
            adc_low.delay       =   round((gx_low.riseTime+gx_low.flatTime/2-time_to_center_low)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
            % adc_low.id          =   2; % for separate adc events

            gy_parts_low                        =   mr.splitGradientAt(gy_low, blip_dur_low/2, lims);
            [gy_blipup_low, gy_blipdown_low,~]  =   mr.align('right',gy_parts_low(1),'left',gy_parts_low(2),gx_low);
            gy_blipdownup_low                   =   mr.addGradients({gy_blipdown_low, gy_blipup_low}, lims);
            % pe_enable support
            gy_blipdown_low.waveform            =   gy_blipdown_low.waveform*pe_enable;
            gy_blipup_low.waveform              =   gy_blipup_low.waveform*pe_enable;
            gy_blipdownup_low.waveform          =   gy_blipdownup_low.waveform*pe_enable;

            Ny_low_pre                          =   round(Ny_low/2); % PE steps prior to ky=0, excluding the central line
            Ny_low_pre                          =   round(Ny_low_pre/RSegment/R_low);
            Ny_low_post                         =   round((partFourierFactor-1/2)*Ny_low+1); % PE lines after the k-space center including the central line
            Ny_low_post                         =   round(Ny_low_post/RSegment/R_low);
            Ny_low_meas                         =   Ny_low_pre+Ny_low_post;

            gxPre_low                           =   mr.makeTrapezoid('x',lims,'Area',-gx_low_org.area/2);
            gyPre_low                           =   mr.makeTrapezoid('y',lims,'Area',(Ny_low_pre*deltaky_low-(Nmulti-1)*R_low*deltak));
            [gxPre_low,gyPre_low]               =   mr.align('right',gxPre_low,'left',gyPre_low);
            gyPre_low                           =   mr.makeTrapezoid('y',lims,'Area',gyPre_low.area,'Duration',mr.calcDuration(gxPre_low,gyPre_low));
            gyPre_low.amplitude                 =   gyPre_low.amplitude*pe_enable;

            gyPre_low_post_area                 =   -1*(Ny_low_post-0.5)*deltaky_low;

            if(mod(Nmulti,2) ~= 0)
                gyPre_low_post_area             =   (Ny_low_post-1)*deltaky_low;
            end

            gxPre_low_post_area_nav             =   -1*(0.5*gx_low.flatArea + gx_low.amplitude*gx_low.fallTime/2); % partFourierFactor should be larger than 0.5
            gxPre_low_post_area                 =   -1*(abs(gxPre_low_post_area_nav*2) - (0.5*gx_low_org.flatArea + gx_low_org.amplitude*gx_low_org.fallTime/2)); % partFourierFactor should be larger than 0.5
            gxPre_low_post_area_nav             =   -1*(abs(gxPre_low_post_area_nav*2) - (0.5*gx_low_org.flatArea + gx_low_org.amplitude*gx_low_org.fallTime/2)); % partFourierFactor should be larger than 0.5

            if(mod(Ny_low_meas,2) == 0)
                gxPre_low_post_area             =   0.5*gx_low_org.flatArea + gx_low_org.amplitude*gx_low_org.riseTime/2;
            end

            gxPre_low_post                      =   mr.makeTrapezoid('x',lims,'Area', gxPre_low_post_area, 'Duration', 0.0012);
            gyPre_low_post                      =   mr.makeTrapezoid('y',lims,'Area', gyPre_low_post_area, 'Duration', 0.0012);
            gzPre_low_post                      =   mr.makeTrapezoid('z',lims,'Area', -kz*(caipi_factor-1)/2, 'Duration', 0.0012);

            if(mod(Ny_low_meas,2) ~= 0)
                gzPre_low_post.amplitude        =   -gzPre_low_post.amplitude;
            end

            [gxPre_low_post,gyPre_low_post,gzPre_low_post] =   mr.align('right',gxPre_low_post,'left',gyPre_low_post,gzPre_low_post);
            gyPre_low_post                      =   mr.makeTrapezoid('y',lims,'Area',gyPre_low_post.area,'Duration',mr.calcDuration(gxPre_low_post,gyPre_low_post,gzPre_low_post));

            gyPre_low_post.amplitude            =   gyPre_low_post.amplitude*pe_enable;
            gzPre_low_post.amplitude            =   gzPre_low_post.amplitude*pe_enable;
            
            if (buda)
                if (mod(Nmulti,2) == 0)
                    gyPre_low           =   mr.scaleGrad(gyPre_low,-1);
                    gy_blipdownup_low   =   mr.scaleGrad(gy_blipdownup_low,-1);
                    gy_blipdown_low     =   mr.scaleGrad(gy_blipdown_low,-1);
                    gy_blipup_low       =   mr.scaleGrad(gy_blipup_low,-1);
                end
            end
            
            
            %% gZBlip

            gzBlip                  = mr.makeTrapezoid('z',lims,'Area',kz,'Duration',blip_dur);
            gzBlip2                 = mr.makeTrapezoid('z',lims,'Area',kz*(caipi_factor-1),'Duration',blip_dur);
            
            gzBlip_parts1           = mr.splitGradientAt(gzBlip, blip_dur/2, lims);
            gzBlip_parts2           = mr.splitGradientAt(gzBlip2, blip_dur/2, lims);
            [gzBlip_l1,gzBlip_r1,~] = mr.align('right',gzBlip_parts1(1),'left',gzBlip_parts1(2),gx);
            [gzBlip_l2,gzBlip_r2,~] = mr.align('right',gzBlip_parts2(1),'left',gzBlip_parts2(2),gx);

            gzBlip_neg1             = mr.addGradients({mr.scaleGrad(gzBlip_r1,-1), mr.scaleGrad(gzBlip_l1,-1)}, lims);
            gzBlip_pos1             = mr.addGradients({mr.scaleGrad(gzBlip_r1,-1), mr.scaleGrad(gzBlip_l2,1)}, lims);
            gzBlip_neg2             = mr.addGradients({mr.scaleGrad(gzBlip_r2,1), mr.scaleGrad(gzBlip_l1,-1)}, lims);

            gzBlipPre               = mr.makeTrapezoid('z',lims,'Area',kz*(caipi_factor-1)/2,'Duration',mr.calcDuration(gyPre));
            gzBlipPre.amplitude     = gzBlipPre.amplitude*pe_enable;

            gzBlip_neg1.waveform    = gzBlip_neg1.waveform*pe_enable;
            gzBlip_neg1.first       = gzBlip_neg1.first*pe_enable;
            gzBlip_neg1.last        = gzBlip_neg1.last*pe_enable;

            gzBlip_pos1.waveform    = gzBlip_pos1.waveform*pe_enable;
            gzBlip_pos1.first       = gzBlip_pos1.first*pe_enable;
            gzBlip_pos1.last        = gzBlip_pos1.last*pe_enable;

            gzBlip_neg2.waveform    = gzBlip_neg2.waveform*pe_enable;
            gzBlip_neg2.first       = gzBlip_neg2.first*pe_enable;
            gzBlip_neg2.last        = gzBlip_neg2.last*pe_enable;

            gzBlip_l1.waveform      = gzBlip_l1.waveform*pe_enable;
            gzBlip_l1.first         = gzBlip_l1.first*pe_enable;
            gzBlip_l1.last          = gzBlip_l1.last*pe_enable;

            gzBlip_l2.waveform      = gzBlip_l2.waveform*pe_enable;
            gzBlip_l2.first         = gzBlip_l2.first*pe_enable;
            gzBlip_l2.last          = gzBlip_l2.last*pe_enable;

            gzBlip_r1.waveform      = gzBlip_r1.waveform*pe_enable;
            gzBlip_r1.first         = gzBlip_r1.first*pe_enable;
            gzBlip_r1.last          = gzBlip_r1.last*pe_enable;

            gzBlip_r2.waveform      = gzBlip_r2.waveform*pe_enable;
            gzBlip_r2.first         = gzBlip_r2.first*pe_enable;
            gzBlip_r2.last          = gzBlip_r2.last*pe_enable;


            %% gZBlip_low

            gzBlip_low                  = mr.makeTrapezoid('z',lims,'Area',kz,'Duration',blip_dur_low);
            gzBlip2_low                 = mr.makeTrapezoid('z',lims,'Area',kz*(caipi_factor-1),'Duration',blip_dur_low);

            gzBlip_low_parts1           = mr.splitGradientAt(gzBlip_low, blip_dur_low/2, lims);
            gzBlip_low_parts2           = mr.splitGradientAt(gzBlip2_low, blip_dur_low/2, lims);
            [gzBlip_low_l1,gzBlip_low_r1,~] = mr.align('right',gzBlip_low_parts1(1),'left',gzBlip_low_parts1(2),gx_low);
            [gzBlip_low_l2,gzBlip_low_r2,~] = mr.align('right',gzBlip_low_parts2(1),'left',gzBlip_low_parts2(2),gx_low);

            gzBlip_low_neg1             = mr.addGradients({mr.scaleGrad(gzBlip_low_r1,-1), mr.scaleGrad(gzBlip_low_l1,-1)}, lims);
            gzBlip_low_pos1             = mr.addGradients({mr.scaleGrad(gzBlip_low_r1,-1), mr.scaleGrad(gzBlip_low_l2,1)}, lims);
            gzBlip_low_neg2             = mr.addGradients({mr.scaleGrad(gzBlip_low_r2,1), mr.scaleGrad(gzBlip_low_l1,-1)}, lims);

            gzBlipPre_low               = mr.makeTrapezoid('z',lims,'Area',kz*(caipi_factor-1)/2,'Duration',mr.calcDuration(gyPre_low));
            gzBlipPre_low.amplitude     = gzBlipPre_low.amplitude*pe_enable;

            gzBlip_low_neg1.waveform    = gzBlip_low_neg1.waveform*pe_enable;
            gzBlip_low_neg1.first       = gzBlip_low_neg1.first*pe_enable;
            gzBlip_low_neg1.last        = gzBlip_low_neg1.last*pe_enable;

            gzBlip_low_pos1.waveform    = gzBlip_low_pos1.waveform*pe_enable;
            gzBlip_low_pos1.first       = gzBlip_low_pos1.first*pe_enable;
            gzBlip_low_pos1.last        = gzBlip_low_pos1.last*pe_enable;

            gzBlip_low_neg2.waveform    = gzBlip_low_neg2.waveform*pe_enable;
            gzBlip_low_neg2.first       = gzBlip_low_neg2.first*pe_enable;
            gzBlip_low_neg2.last        = gzBlip_low_neg2.last*pe_enable;

            gzBlip_low_l1.waveform      = gzBlip_low_l1.waveform*pe_enable;
            gzBlip_low_l1.first         = gzBlip_low_l1.first*pe_enable;
            gzBlip_low_l1.last          = gzBlip_low_l1.last*pe_enable;

            gzBlip_low_l2.waveform      = gzBlip_low_l2.waveform*pe_enable;
            gzBlip_low_l2.first         = gzBlip_low_l2.first*pe_enable;
            gzBlip_low_l2.last          = gzBlip_low_l2.last*pe_enable;

            gzBlip_low_r1.waveform      = gzBlip_low_r1.waveform*pe_enable;
            gzBlip_low_r1.first         = gzBlip_low_r1.first*pe_enable;
            gzBlip_low_r1.last          = gzBlip_low_r1.last*pe_enable;

            gzBlip_low_r2.waveform      = gzBlip_low_r2.waveform*pe_enable;
            gzBlip_low_r2.first         = gzBlip_low_r2.first*pe_enable;
            gzBlip_low_r2.last          = gzBlip_low_r2.last*pe_enable;


            %% Calculate delay times

            if (sum(table(rep,:))==0) % need to adjust delta TE1 and delta TE2 for b0, bcuz of the crusher
                crusher_switch=1;
            else
                crusher_switch=0;
            end

            % Calculate delay times
            durationToCenter        = (Ny_pre + 0.5) * mr.calcDuration(gx);
            durationToCenter2       = (Ny_post + 0.5) * mr.calcDuration(gx);
            durationToCenter_low    = (Ny_low_pre + 0.5) * mr.calcDuration(gx_low);

            if rf_type == 1
                rfCenterInclDelay       = rf(rf_order(igSlider)).delay + mr.calcRfCenter_QL(rf(rf_order(igSlider)));
                rf180centerInclDelay    = rf180.delay + mr.calcRfCenter_QL(rf180);
            elseif rf_type == 2
                rfCenterInclDelay       = rf.delay + mr.calcRfCenter(rf);
                rf180centerInclDelay    = rf180.delay + mr.calcRfCenter(rf180);
            end

            delayTE1                = ceil((TE/2 - mr.calcDuration(rf,gz) - mr.calcDuration(gzReph) + rfCenterInclDelay - rf180centerInclDelay)/lims.gradRasterTime)*lims.gradRasterTime;
            delayTE1                = delayTE1 - crusher_switch * mr.calcDuration(gz180_crusher_1);
            assert(delayTE1>=0);

            gxPre.delay=0;
            gyPre.delay=0;
            delayTE2                = ceil((TE/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - durationToCenter)/lims.gradRasterTime)*lims.gradRasterTime;
            delayTE2                = delayTE2 - mr.calcDuration(gxPre,gyPre);
            delayTE2                = delayTE2 - crusher_switch*mr.calcDuration(gz180_crusher_1);
            assert(delayTE2>=0);

            delayTE3                = ceil((TE_low/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - durationToCenter2)/lims.gradRasterTime)*lims.gradRasterTime;
            delayTE3                = delayTE3 - mr.calcDuration(gxPre_post,gyPre_post,gzPre_post) - mr.calcDuration(gz180_c1);
            delayTE3                = delayTE3 - 1e-4;
            assert(delayTE3>=0);

            delayTE4                = ceil((TE_low/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - durationToCenter_low)/lims.gradRasterTime)*lims.gradRasterTime;
            delayTE4                = delayTE4 - mr.calcDuration(gxPre_low,gyPre_low);
            delayTE4                = delayTE4 - mr.calcDuration(gz180_c1);
            assert(delayTE4>=0);

            [gxPre,gyPre]           = mr.align('right',gxPre,'left',gyPre);

            % diffusion weithting calculation
            % delayTE2 is our window for small_delta
            % delayTE1+delayTE2-delayTE2 is our big delta
            % we anticipate that we will use the maximum gradient amplitude, so we need
            % to shorten delayTE2 by gmax/max_sr to accommodate the ramp down
            small_delta = delayTE2 - ceil(lims2.maxGrad/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
            big_delta   = delayTE1 + mr.calcDuration(rf180,gz180); 
            % we define bFactCalc function below to eventually calculate time-optimal
            % gradients. for now we just abuse it with g=1 to give us the coefficient
            % g=sqrt(bFactor*1e6/bFactCalc(1,small_delta,big_delta)); % for now it looks too large!
            g           = sqrt(bFactor(1,rep)*1e6/bFactCalc(1,small_delta,big_delta));

            gr          = ceil(g/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
            gDiff       = mr.makeTrapezoid('z','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',lims2);


            %% split the diffusion gradient to 3 axis

            g_x=g.*table(rep,1);
            g_y=g.*table(rep,2);
            g_z=g.*table(rep,3);

            if ((sum(table(rep,:))==0)||(sum(table(rep,:))==1))
                g_xr=ceil(abs(g_x)/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
                g_yr=ceil(abs(g_y)/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
                g_zr=ceil(abs(g_z)/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
                gDiff_x=mr.makeTrapezoid('x','amplitude',g_x,'riseTime',g_xr,'flatTime',small_delta-g_xr,'system',lims2);
                gDiff_y=mr.makeTrapezoid('y','amplitude',g_y,'riseTime',g_yr,'flatTime',small_delta-g_yr,'system',lims2);
                gDiff_z=mr.makeTrapezoid('z','amplitude',g_z,'riseTime',g_zr,'flatTime',small_delta-g_zr,'system',lims2);
                duration_diff=max(max(mr.calcDuration(gDiff_x), mr.calcDuration(gDiff_y)), mr.calcDuration(gDiff_z));

            else % 3 axis, 'rotation'
                [azimuth,elevation,r] = cart2sph(g_x,g_y,g_z);
                polar= -(pi/2-elevation);

                Gr=mr.rotate('z',azimuth,mr.rotate('y',polar,gDiff));
                if size(Gr,2)==3
                    gDiff_x=Gr{1,2};
                    gDiff_y=Gr{1,3};
                    gDiff_z=Gr{1,1};
                else
                    if size(Gr,2)==2
                        diffusion_blank=find( table(rep,:)==0);
                        switch diffusion_blank
                            case 2
                                gDiff_x=Gr{1,2};
                                gDiff_z=Gr{1,1};
                                gDiff_y=gDiff; gDiff_y.channel='y'; gDiff_y.amplitude=0; gDiff_y.area=0; gDiff_y.flatArea=0;
                            case 1
                                gDiff_z=Gr{1,1};
                                gDiff_y=Gr{1,2};
                                gDiff_x=gDiff; gDiff_x.amplitude=0; gDiff_x.area=0; gDiff_x.flatArea=0;gDiff_x.channel='x';
                            case 3
                                gDiff_x=Gr{1,2};
                                gDiff_y=Gr{1,1};
                                gDiff_z=gDiff; gDiff_z.amplitude=0; gDiff_z.area=0; gDiff_z.flatArea=0;gDiff_z.channel='z';
                        end
                    end
                end
                duration_diff=mr.calcDuration(gDiff);
            end

            assert(duration_diff<=delayTE1);
            assert(duration_diff<=delayTE2);


            %% Calculate the echo time shift for multishot EPI

            actual_esp          = gx.riseTime + gx.flatTime + gx.fallTime;
            TEShift             = actual_esp/RSegment;
            TEShift             = round(TEShift,5);
            TEShift_before_echo = (Nmulti-1)*TEShift;
            if TEShift_before_echo == 0
                TEShift_before_echo = 0.00001;
            end
            TEShift_after_echo  = (RSegment-(Nmulti-1))*TEShift;
            dETS_before         = mr.makeDelay(TEShift_before_echo);
            dETS_after          = mr.makeDelay(TEShift_after_echo);

            actual_low_esp          = gx_low.riseTime + gx_low.flatTime + gx_low.fallTime;
            TEShift_low             = actual_low_esp/RSegment;
            TEShift_low             = round(TEShift_low,5);
            TEShift_low_before_echo = (Nmulti-1)*TEShift_low;
            if TEShift_low_before_echo == 0
                TEShift_low_before_echo = 0.00001; % apply the minimum duration for the no delay case
            end
            TEShift_low_after_echo  = (RSegment-(Nmulti-1))*TEShift_low;
            dETS_low_before         = mr.makeDelay(TEShift_low_before_echo);
            dETS_low_after          = mr.makeDelay(TEShift_low_after_echo);


            %% TR calculation

            TR_slc      =   TR/Nslices; % TR per slice
            delayTR     =   ceil((TR_slc - mr.calcDuration(gp_r) - mr.calcDuration(gn_r) - rfCenterInclDelay - TE - TE_low - Ny_low_post*mr.calcDuration(gx_low) - mr.calcDuration(gyPre_low_post))/seq.gradRasterTime)*seq.gradRasterTime;
            delayTR     =   round(delayTR,3); % pulseq will allow delay time at power -3
            dTR         =   mr.makeDelay(delayTR);


            %% Define Seq Block

            if test_check == 1
                slc_end = 1;
            else
                slc_end = Nslices;
            end
            for slc = 1:slc_end

                % 90 RF Module Define
                seq.addBlock(gp_r,gp_p,gp_s);
                seq.addBlock(rf_fs,gn_r,gn_p,gn_s);

                if rf_type == 1
                    if MB == 2
                        rf(rf_order(igSlider)).signal   = p1(:,slc).*rf(rf_order(igSlider)).signal;
                        rf180.signal                    = p2(:,slc).*rf180.signal;
                        rf180.phaseOffset               = pi/2;
                    elseif MB == 1
                        rf(igSlider).freqOffset         = max(gz.waveform)*thickness*interleaved_freqOffset_factor(slc);
                        rf(igSlider).phaseOffset        = 0;
                        rf180.freqOffset                = max(gz180.waveform)*thickness*interleaved_freqOffset_factor(slc);
                        rf180.phaseOffset               = pi/2;
                    end

                    rf(igSlider).signal(1)  = 0;
                    rf(igSlider).signal(end)= 0;

                    gz.waveform(1)          = 0;
                    gz.waveform(end)        = 0;
                    
                    rf180.signal(1)         = 0;
                    rf180.signal(end)       = 0;
                    
                    gz180.waveform(1)       = 0;
                    gz180.waveform(end)     = 0;
                end

                if rf_type == 2
                    rf.freqOffset           = gz.amplitude*thickness*interleaved_freqOffset_factor(slc);
                    rf180.freqOffset        =   gz180.amplitude*thickness*interleaved_freqOffset_factor(slc);
                    rf180.phaseOffset       =   pi/2;
                end
                
                if rf_type == 2
                    seq.addBlock(rf,gz,trig);
                elseif rf_type == 1
                    seq.addBlock(rf(rf_order(igSlider)),gz,trig);
                end
                seq.addBlock(gzReph);

                % diff gradient
                seq.addBlock(mr.makeDelay(delayTE1),gDiff_x,gDiff_y,gDiff_z);

                % 180 RF Module Define
                if (sum(table(rep,:))==0)
                    if(~recon)
                        seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                    end
                end

                seq.addBlock(rf180,gz180);

                if (sum(table(rep,:))==0)
                    if(~recon)
                        seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                    end
                end

                % diff gradient
                seq.addBlock(mr.makeDelay(delayTE2),gDiff_x, gDiff_y, gDiff_z);

                % for no blip
                if rep == 1
                    gy_blipup.last=0;
                    gy_blipdown.first=0;
                    gy_blipdownup.last=0;
                    gy_blipdownup.first=0;
                end
                
                if (Echotimeshift)
                    seq.addBlock(dETS_before); % echotimeshift for multishot
                end
                
                if CAIPI
                    gzPre = gzBlipPre; % gzReph
                end
                
                if CAIPI
                    seq.addBlock(mr.align('left',gyPre,gzPre,'right',gxPre));
                else
                    seq.addBlock(gxPre,gyPre);
                end

                for i = 1:Ny_meas % loop over phase encodes
                    gx_temp = mr.scaleGrad(gx,(-1)^(i-1));
                    if i == 1
                        if CAIPI
                            step = 0;
                            seq.addBlock(gx_temp,gy_blipup,mr.scaleGrad(gzBlip_l1,-1),adc); % Read the first line of k-space with a single half-blip at the end
                        else
                            seq.addBlock(gx_temp,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
                        end
                    elseif i==Ny_meas
                        if CAIPI
                            step = step + 1;
                            switch mod(step,caipi_factor)
                                case num2cell(1:caipi_factor-2)
                                    seq.addBlock(gx_temp,gy_blipdown,mr.scaleGrad(gzBlip_r1,-1),adc);
                                case caipi_factor-1
                                    seq.addBlock(gx_temp,gy_blipdown,mr.scaleGrad(gzBlip_r1,-1),adc);
                                case 0
                                    seq.addBlock(gx_temp,gy_blipdown,mr.scaleGrad(gzBlip_r2,1),adc);
                            end
                        else
                            seq.addBlock(gx_temp,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
                        end
                    else
                        if CAIPI
                            step = step+1;
                            switch mod(step,caipi_factor)
                                case num2cell(1:caipi_factor-2)
                                    seq.addBlock(gx_temp,gy_blipdownup,gzBlip_neg1,adc);
                                case caipi_factor-1
                                    seq.addBlock(gx_temp,gy_blipdownup,gzBlip_pos1,adc);
                                case 0
                                    seq.addBlock(gx_temp,gy_blipdownup,gzBlip_neg2,adc);
                            end
                        else
                            seq.addBlock(gx_temp,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                        end
                    end
                end

                seq.addBlock(gxPre_post,gyPre_post,gzPre_post); % go back to k-space center

                if (Echotimeshift)
                    seq.addBlock(dETS_after); % echotimeshift for multishot
                end


                %% 2nd echo

                gyPre_slt           =   gyPre_low;
                gy_blipup_slt       =   gy_blipup_low;
                gy_blipdown_slt     =   gy_blipdown_low;
                gy_blipdownup_slt   =   gy_blipdownup_low;
                gyPre_post_slt      =   gyPre_low_post;

                seq.addBlock(mr.makeDelay(delayTE3));
                seq.addBlock(gz180_c1,gz180_c2,gz180_c3);
                seq.addBlock(rf180,gz180);
                seq.addBlock(gz180_c1,gz180_c2,gz180_c3);

                seq.addBlock(mr.makeDelay(delayTE4));

                % for no blip
                if rep == 1
                    gy_blipup_slt.last      =   0;
                    gy_blipdown_slt.first   =   0;
                    gy_blipdownup_slt.last  =   0;
                    gy_blipdownup_slt.first =   0;
                end

                if (Echotimeshift)
                    seq.addBlock(dETS_low_before); % echotimeshift for multishot
                end

                if CAIPI
                    gzPre_low = gzBlipPre_low; % gzReph
                end

                if CAIPI
                    seq.addBlock(mr.align('left',gyPre_slt,gzPre_low,'right',gxPre_low));
                else
                    seq.addBlock(gxPre_low,gyPre_slt);
                end

                for i = 1:Ny_low_meas % loop over phase encodes
                    gx_low_temp = mr.scaleGrad(gx_low,(-1)^(i-1));
                    if i == 1
                        if CAIPI
                            step = 0;
                            seq.addBlock(gx_low_temp,gy_blipup_slt,mr.scaleGrad(gzBlip_low_l1,-1),adc_low); % Read the first line of k-space with a single half-blip at the end
                        else
                            seq.addBlock(gx_low_temp,gy_blipup_slt,adc_low); % Read the first line of k-space with a single half-blip at the end
                        end
                    elseif i==Ny_low_meas
                        if CAIPI
                            step = step + 1;
                            switch mod(step,caipi_factor)
                                case num2cell(1:caipi_factor-2)
                                    seq.addBlock(gx_low_temp,gy_blipdown_slt,mr.scaleGrad(gzBlip_low_r1,-1),adc_low);
                                case caipi_factor-1
                                    seq.addBlock(gx_low_temp,gy_blipdown_slt,mr.scaleGrad(gzBlip_low_r1,-1),adc_low);
                                case 0
                                    seq.addBlock(gx_low_temp,gy_blipdown_slt,mr.scaleGrad(gzBlip_low_r2,1),adc_low);
                            end
                        else
                            seq.addBlock(gx_low_temp,gy_blipdown_slt,adc_low); % Read the last line of k-space with a single half-blip at the beginning
                        end
                    else
                        if CAIPI
                            step = step+1;
                            switch mod(step,caipi_factor)
                                case num2cell(1:caipi_factor-2)
                                    seq.addBlock(gx_low_temp,gy_blipdownup_slt,gzBlip_low_neg1,adc_low);
                                case caipi_factor-1
                                    seq.addBlock(gx_low_temp,gy_blipdownup_slt,gzBlip_low_pos1,adc_low);
                                case 0
                                    seq.addBlock(gx_low_temp,gy_blipdownup_slt,gzBlip_low_neg2,adc_low);
                            end
                        else
                            seq.addBlock(gx_low_temp,gy_blipdownup_slt,adc_low); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                        end
                    end
                end

                seq.addBlock(gxPre_low_post,gyPre_post_slt,gzPre_low_post); % go back to k-space center

                if (Echotimeshift)
                    seq.addBlock(dETS_low_after); % echotimeshift for multishot
                end

                seq.addBlock(dTR);

                if rf_type == 1
                    clear rf rf180
                    load('seq_blocks_MB2_1mm_s26_PR.mat') % rf, gz blocks
                end

            end % slice loop
        end % gSlider loop
    end % multishot loop
    disp(['diffusion dir:', num2str(rep)])
    clear gDiff gDiff_x gDiff_y gDiff_z
end % diffusion loop


%% check whether the timing of the sequence is correct

[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end


%% do some visualizations

if seq_plot == 1
    tic
    fprintf('Plotting figures... ');

    seq.plot(); % plot sequence waveforms
    
    % trajectory calculation
    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();
    
    % plot k-spaces
    figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
    hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
    
    figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
    axis('equal'); % enforce aspect ratio for the correct trajectory display
    hold on; plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
    toc
end


%% prepare the sequence output for the scanner

seq.setDefinition('FOV', [fov fov thickness*Nslices*MB]);
seq.setDefinition('Name', 'epi-diff');

if seq_pns == 1
    tic
    fprintf('Checking PNS... ');
    % pns < 80\% (human)
    [pns_ok, pns_n, pns_c, tpns] = seq.calcPNS('file.asc');
    pns = max(pns_n(:));
    
    clear ktraj_adc t_adc ktraj t_ktraj t_excitation t_refocusing slicepos t_slicepos pns_n pns_c tpns
    toc
end

% if pns_ok
    if save_seq_file == 1
        tic
        fprintf('Saving .seq file... ');
        seq.write(seq_file);
        save(strcat(seq_file(1:end-4),'_param'));
        toc
    end
% end

% seq.install('siemens');
% seq.sound(); % simulate the seq's tone


%% check forbidden frequencies

if check_freq == 1
    tic
    fprintf('Checking frequencies... ');
    sys = lims;
    gradSpectrum;
    toc
end


%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  

if seq_report == 1
    tic
    fprintf('Checking test report... ');
    rep = seq.testReport;
    fprintf([rep{:}]);
    toc
end


%%

function b=bFactCalc(g, delta, DELTA)
% see DAVY SINNAEVE Concepts in Magnetic Resonance Part A, Vol. 40A(2) 39????5 (2012) DOI 10.1002/cmr.a
% b = gamma^2  g^2 delta^2 sigma^2 (DELTA + 2 (kappa - lambda) delta)
% in pulseq we don't need gamma as our gradinets are Hz/m
% however, we do need 2pi as diffusion equations are all based on phase
% for rect gradients: sigma=1 lambda=1/2 kappa=1/3
% for trapezoid gradients: TODO
sigma=1;
%lambda=1/2;
%kappa=1/3;
kappa_minus_lambda=1/3-1/2;
b= (2*pi * g * delta * sigma)^2 * (DELTA + 2*kappa_minus_lambda*delta);
end

function b = bFactCalc_trap(g, delta, DELTA, rup)
% handbook, page 287, example 9.2
b = (2*pi*g)^2 * (delta^2*(DELTA-delta/3)  +  rup^3/30-delta*rup^2/6);
end