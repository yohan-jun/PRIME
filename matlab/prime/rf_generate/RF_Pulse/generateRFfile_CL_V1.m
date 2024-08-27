%% generate RF pulse
% final version
% Fuyixue Wang, February, 2017
% Congyu Liao, Nov. 2017 
clear;
close all;
clc;
%%
% linear phase
profile_filename = 'Profiles/RF_5xstandard_1p00t_4p3mm_Verse1_CLv1_linPhase_gSliderRF_5x_exlinPhs.mat';% load the gSlider5x 90 and 180 pulses
 



 rf_filename = 'generated_file/test_QL_RF_5xstandard_1p03t_4p3mm_Verse1_CLv1_linPha_gSliderRF_5x_exlinPhs.float';

VerseDownFac = 1/2.2;  % verse factor
N_sets = 1; % sms scaling_fac=1.0;
N_shots=1;
%%
[rf_90_set_wRefocus,g_90_versed_10us,AsymmetryFac_90,rf_180_versed,g_180_versed_10us,AsymmetryFac_180,slThick]=PulsePrep_v3(profile_filename,VerseDownFac);
%%
% Length
[L_rf90,N_Enc]=size(rf_90_set_wRefocus);
%% add multi-shot version
N_Enc_new=N_Enc*N_shots;
rf_90_set_wRefocus_new=zeros(L_rf90,N_Enc_new);
for ii=1:N_Enc
    for jj=1:N_shots
            rf_90_set_wRefocus_new(:,(ii-1)*N_shots+jj)=rf_90_set_wRefocus(:,ii);
    end
end
if (0)
    for kk=1: N_Enc_new
        figure(110);
        subplot(N_Enc_new,1,kk);hold on;
        plot(real(rf_90_set_wRefocus_new(:,kk)));
        plot(imag(rf_90_set_wRefocus_new(:,kk)),'r');
        plot(abs(rf_90_set_wRefocus_new(:,kk)),'g');title('multi shot rf excitation')
        title('multi shot rf versed'); 
    end 
end

clear N_Enc rf_90_set_wRefocus ii jj kk
N_Enc=N_Enc_new;  rf_90_set_wRefocus=rf_90_set_wRefocus_new;
clear N_Enc_new rf_90_set_wRefocus_new
%%
L_rf180 = length(rf_180_versed);
% abs and phase of Encoding pulses
angle_rf90 = angle(rf_90_set_wRefocus(:));
index_rf90 = angle_rf90<0;
angle_rf90(index_rf90) = angle_rf90(index_rf90)+2*pi*0.999;
angle_rf90 = reshape(angle_rf90,[L_rf90,N_Enc]);
% other pulse
angle_rf180 = angle(rf_180_versed(:));
index_rf180 = angle_rf180<0;
angle_rf180(index_rf180) = angle_rf180(index_rf180)+2*pi*0.999;
% gradient
L_g90 = length(g_90_versed_10us);
L_g180 = length(g_180_versed_10us);
g_90 = zeros(L_g90,3);
g_90(:,3)=g_90_versed_10us(:);
g_180 = zeros(L_g180,3);
g_180(:,3)=g_180_versed_10us(:);
%%
scale_g_90 = max(abs(g_90(:)));
scale_g_180 = max(abs(g_180(:)));
g_90 = g_90./ scale_g_90;
g_180 = g_180./scale_g_180;
%%
fidx = fopen(rf_filename,'w');
% length
fwrite(fidx,slThick,'float32'); % Munimum allowed thickness
fwrite(fidx,N_sets,'float32');      % sms
fwrite(fidx,N_Enc,'float32');       % number of encoding
fwrite(fidx,L_rf90,'float32');        % number of points for 90 pulses
fwrite(fidx,L_rf180,'float32');       % number of points for 180 pulses
fwrite(fidx,L_g90,'float32');        % number of points for gradient of 90
fwrite(fidx,L_g180,'float32');       % number of points for gradient of 180 pulses
% gradient
fwrite(fidx,scale_g_90,'float32');
fwrite(fidx,scale_g_180,'float32');
fwrite(fidx,g_90(:),'float32');
fwrite(fidx,g_180(:),'float32');
%%

rf_normalized_90 = zeros(size(rf_90_set_wRefocus));
rf_normalized_180 = zeros(size(rf_180_versed));
scalFac_90 = zeros(N_Enc,1);
scalFac_180 = zeros(N_Enc,1);

% realpulse = rf_90_set_wRefocus(:,3)./max(abs(rf_90_set_wRefocus(:,3)));
% AmpIntegral_real = sqrt(sum(real(realpulse)).^2+sum(imag(realpulse)).^2);
% max_real = max(abs(rf_90_set_wRefocus(:,3)));
AmpIntegral_real = 761.4774;
max_real = 4.2795;
% [max(abs(rf_90_set_wRefocus)), max(abs(rf_180_versed))]./max(abs(rf_90_set_wRefocus(:,1)))
for dif_Enc=1:N_Enc
    rf_90_temp = rf_90_set_wRefocus(:,dif_Enc);
    normalized_temp = rf_90_temp./max(abs(rf_90_temp(:)));
    
    rf_normalized_90 (:,dif_Enc) = normalized_temp;
    scalFac_90(dif_Enc) = (sqrt(sum(real(normalized_temp)).^2+sum(imag(normalized_temp)).^2)/AmpIntegral_real*max(abs(rf_90_temp))/max_real*(100.96/90))*scaling_fac;
    
    rf_normalized_180 = rf_180_versed./max(abs(rf_180_versed(:)));
    scalFac_180(dif_Enc) = (sqrt(sum(real(rf_normalized_180)).^2+sum(imag(rf_normalized_180)).^2)/AmpIntegral_real*max(abs(rf_180_versed))/max_real*(100.96/180))*scaling_fac;
end
%% 
for dif_set = 1:N_sets
    fwrite(fidx,[AsymmetryFac_90*ones(N_Enc,1);AsymmetryFac_180*ones(N_Enc,1)],'float32'); % AsymmetryFactor
    fwrite(fidx,[scalFac_90;scalFac_180],'float32'); %scaling factor
    % pulses
    for dif_Enc=1:N_Enc
        rf_mag = abs(rf_normalized_90(:,dif_Enc));
        fwrite(fidx,rf_mag,'float32');
        fwrite(fidx,angle_rf90(:,dif_Enc),'float32');
    end

    for dif_Enc=1:N_Enc
        rfOther_mag = abs(rf_normalized_180);
        fwrite(fidx,rfOther_mag,'float32');
        fwrite(fidx,angle_rf180,'float32');
    end
end

fclose(fidx);




% 1.0000    0.8559    0.8220    0.8559    1.0000    2.7119
[max(abs(rf_90_set_wRefocus)), max(abs(rf_180_versed))]./max(abs(rf_90_set_wRefocus(:,1)))
