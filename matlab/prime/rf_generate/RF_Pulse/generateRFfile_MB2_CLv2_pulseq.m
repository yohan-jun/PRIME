%% generate RF pulse
% final version
% Kawin Setsompop: March, 2017
% Fuyixue Wang, February, 2017
% Congyu Liao, Nov, 2017
% Congyu Liao, May, 2021, for GE scanners
clear;
close all;
clc;
%%
% linear phase
% profile_filename = 'Profiles/RF_1xstandard_1p15t_1p0mm_Verse1_CLv1_linPhase_exlinPhs.mat';  % load the SLR 90 and 180 pulses
profile_filename = '/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/RF_5x_linPhase_1p03_5mm_rf_4us_gSliderRF_1129_5x_exlinPhs.mat';% load the gSlider5x 90 and 180 pulses
 

VerseDownFac = 1/2.2;  % verse factor
N_sets = 1; % SMS sets, set to 1
scaling_fac=1.0;  % B1 scaling
N_shots=1;  % add multi-shot 


%% do verse for MB1 pulse
[rf_90_set_wRefocus,g_90_versed_4us,AsymmetryFac_90,rf_180_versed,g_180_versed_4us,AsymmetryFac_180,slThick,nomFlips]=PulsePrep_v4_pulseq(profile_filename,VerseDownFac);
%
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

clear N_Enc rf_90_set_wRefocus ii jj kk
N_Enc=N_Enc_new;  rf_90_set_wRefocus=rf_90_set_wRefocus_new;
clear N_Enc_new rf_90_set_wRefocus_new

%% save it to GE scanner
% save_rf_g(g,rf,frq,ang,thk,isodelay, format, ts_rf, ts_g,root_fname)

% SS_SAVE - Save spectral-spatial pulse
% Uses Chuck Cunningham's format for GE systems, and creates associated
% .dat-file
% Pulse parameters saved in header for Varian fules
%
%  ss_save(g,rf,ang,thk, isodelay, format, fspec, a_angs)
%
%  g - in G/cm
%  rf - in G
%  ang - flip angle in radians
%  thk - thickness in cm
%  isodelay (optional) - delay from in-phase point to end of pulse (GE
%  definition)
%  format (optional) - 'GE' (default), 'Varian'
%  fspec (optional) - frequency bands (Hz) to write in file
%  a_angs (optional) - band amplidutes (radians) to write in file


%% generate MB=2 version
g_90_versed_4us = slewRate_smooth (g_90_versed_4us);  % smooth the gradients to meet the max slew rate
g_180_versed_4us = slewRate_smooth (g_180_versed_4us);

slice_gap =13*5; % slice gap mm 70 mm
ts_rf =  4*10^-6;  % rf dwelltime= 2us
amplscale = slThick/5;  % I want 5 mm slice thickness on scanner....
[rf_90_set_wRefocus_MB2, rf_180_versed_MB2,g_90_integral,g_180_integral] = generateMB2pulse_CLv1_pulseq(rf_90_set_wRefocus,g_90_versed_4us,rf_180_versed,g_180_versed_4us,slice_gap,ts_rf,amplscale);

frq = [];
thk = slThick/10; %mm-> cm
isodelay = []; % default
format ='GE'; % default

ts_g =  4*10^-6;  % gradient dwelltime = 4us


g_90_integral1 = sprintf('%.5f,',g_90_integral.*(10^3));  % the gradient integral of 90 (Gauss/cm*us) used for slice-selection using RF phases (omega/theata modulation)
g_180_integral1 = sprintf('%.5f,',g_180_integral.*(10^3));


