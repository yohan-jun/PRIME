function [rf_90_set_wRefocus,g_90_versed_4us,g_90_versed_rew_4us,AsymmetryFac_90,rf_180_versed,g_180_versed_4us,AsymmetryFac_180,slThick,nomFlips,g_90_versed,g_180_versed]= PulsePrep_v4_pulseq_YJ(filename,VerseDownFac)

% Kawin Setsompop: March 27th 2017
% Updated by Fuyixue Wang: Aug 31,2017
% Updated by Congyu Liao: Nov 10,2017, after ISMRM submission...
load(filename);

Ts = dt*10^-3;
rf_90_set = rfEncOut; clear rfEncOut
rf_180 = rfOtherOut; clear rfOtherOut
gAmp = gAmp/10; % in G/cm
nEncodes = size(rf_90_set,2);

gmax_90 = gAmp*1.0; %1.0;
gmax_180 = gAmp*1.0;
gmax_refoc = 2.5;
smax = 12000;  % original 12000 for siemens, use 10000 instead for GE scanners;
%% 90
g_90_noRefocus = gAmp*[ones(1,size(rf_90_set,1))];    % max gradient amplitude [G/cm]
g_Refocus = gAmp*ones(1,(round(size(rf_90_set,1)/2)));
g_90 = gAmp*[ones(1,size(rf_90_set,1)) -g_Refocus];    % max gradient amplitude [G/cm]

% rf_90_forVERSE =[floor(zeros((size(rf_90_set,1)-size(rf_180,2))/2,1)); rf_180.'; ceil(zeros((size(rf_90_set,1)-size(rf_180,2))/2,1))];
rf_90_forVERSE =[zeros(floor((size(rf_90_set,1)-size(rf_180,2))/2),1); rf_180.'; zeros(ceil((size(rf_90_set,1)-size(rf_180,2))/2),1)];
[rf_90_versed, g_90_noRefocus_versed] = VERSEwrapper(rf_90_forVERSE,g_90_noRefocus.',gmax_90,smax,Ts,VerseDownFac); % need to watch out for peak values bug
[rf_throwout, g_Refocus_versed] = VERSEwrapper(0.5*ones(round(size(rf_90_set,1)/2),1),g_Refocus.',gmax_refoc,smax,Ts,gmax_refoc/gAmp);
k = cumsum(g_90_noRefocus);
k_versed = cumsum(g_90_noRefocus_versed);
 g_90_versed_rew = [g_90_noRefocus_versed.' -g_Refocus_versed.'];  
g_90_versed = [g_90_noRefocus_versed.'];

for EncodingCount = 1:nEncodes
    rf_90 = rf_90_set(:,EncodingCount).';
    
    %integral = sum(real(rf_90));
    rf_90_versed = spline(k,rf_90,k_versed).*(g_90_noRefocus_versed/gAmp);
    
    
    rf_90_set_noRefocus(:,EncodingCount) = rf_90_versed.';
%     rf_90_set_wRefocus(:,EncodingCount) = cat(2,rf_90_versed.',zeros(1,round(length(rf_90)/2)));
    rf_90_set_wRefocus(:,EncodingCount) = rf_90_versed.';
    
    figure(1);
    subplot(nEncodes,1,EncodingCount);hold on;
    plot(real(rf_90));
    plot(imag(rf_90),'r');
    plot(abs(rf_90),'g');title('RF 90 excitation')
    title('rf90');
    
    figure(2);
    subplot(nEncodes,1,EncodingCount);hold on;
    plot(real(rf_90_versed));
    plot(imag(rf_90_versed),'r');
    plot(abs(rf_90_versed),'g');title('RF 90 excitation')
    title('rf90versed');
    
    figure(3);
    plot(g_90); hold on; plot(g_90_versed,'r');
    title('g90 and g90versed');
    legend('g90','g90versed');
end

AsymmetryFac_90 =  (length(g_90_noRefocus_versed)/2)/length(g_90_versed);

% oversample = 0.4*10^-5/Ts; % QL
oversample=10^-5/Ts;
g_90_versed_4us = g_90_versed(max(floor(oversample/2),1):oversample:end); % quick hack for now
g_90_versed_rew_4us=g_90_versed_rew(max(floor(oversample/2),1):oversample:end); % QL
correctRFlength = length(g_90_versed_4us)*oversample;
currentRFlength = length(rf_90_set_wRefocus(:,1));
L = correctRFlength -currentRFlength;
if L > 0
    rf_90_set_wRefocus = [rf_90_set_wRefocus;zeros(L,nEncodes)];
elseif L<0
    rf_90_set_wRefocus = rf_90_set_wRefocus(1:end+L,:);
end

%% 180

g_180 = gAmp*[ones(size(rf_180))];
[rf_180_versed, g_180_versed] = VERSEwrapper(rf_180.',g_180.',gmax_180,smax,Ts,VerseDownFac);


figure(4);
subplot(2,1,1);
plot(real(rf_180));hold on;
plot(imag(rf_180),'r');
plot(abs(rf_180),'g');title('RF 180 refoc')
subplot(2,1,2);
plot(real(rf_180_versed));hold on;
plot(imag(rf_180_versed),'r');
plot(abs(rf_180_versed),'g');title('RF 180 refoc versed')

figure(5);
plot(g_180); hold on; plot(g_180_versed,'r');
title('g180 and g180versed');
legend('g180','g180versed');

AsymmetryFac_180 =  0.5;
g_180_versed_4us = g_180_versed(max(floor(oversample/2),1):oversample:end); % quick hack for now
correctRFlength = length(g_180_versed_4us)*oversample;
currentRFlength = length(rf_180_versed);
L = correctRFlength -currentRFlength;
if L > 0
    rf_180_versed = [rf_180_versed;zeros(L,1)];
elseif L<0
    rf_180_versed = rf_180_versed(1:end+L,:);
end

end


