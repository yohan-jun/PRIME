
function [rfEncOutmb2, rfOtherOutmb2,g_90_integral,g_180_integral] = generateMB2pulse_CLv1_pulseq_YJ(rfEncOut,g_90,rfOtherOut,g_180,slice_gap,ts_rf,amplscale)
% Generate MB2 pulses for gSlider 
% Congyu Liao, Nov. 2017
% Congyu Liao, May. 2021



% ts_rf =  2*10^-6;  % rf dwelltime= 2us

GAM = 4257.59 ; % gammar  = 4257.59 Hz/Gauss
G = size(rfEncOut,2);

% tbG=12;
% tbOther=8;
% otherThickFactor = 1.00;


rfEncOut1=rfEncOut;
rfOtherOut1=rfOtherOut;

rfEncOut2=rfEncOut;
rfOtherOut2=rfOtherOut;

% sliceShift=70;%mm
% freqShift=tbG/(slThick)*slice_gap; %Hz*s, since tbG/T is Hz. This freqShift is in norminal unit (unit less)
% freqShiftOther=tbOther/(slThick)/otherThickFactor*slice_gap;

%gx = interp1(1:length(gx), gx, 0.25:0.5:length(gx), 'linear', 'extrap').';  % start at midpoint of 2us dwell = 1us

% g_90_intp = [interp1(1:length(g_90(:)),g_90,1:0.5:length(g_90(:))),0];
% g_180_intp = [interp1(1:length(g_180(:)),g_180,1:0.5:length(g_180(:))),0];
% for 4 us, QL
% g_90_intp = [interp1(1:length(g_90(:)),g_90,1:0.4:length(g_90(:))),0];
% g_180_intp = [interp1(1:length(g_180(:)),g_180,1:0.4:length(g_180(:))),0];
g_90_intp = [interp1(1:length(g_90(:)),g_90,1:0.19987:length(g_90(:))),0]; % YJ % use this
g_180_intp = [interp1(1:length(g_180(:)),g_180,1:0.1998:length(g_180(:))),0]; % YJ % use this
% g_90_intp = [interp1(1:length(g_90(:)),g_90,1:0.1996:length(g_90(:))),0]; % YJ % 1mm / 3ms
% g_180_intp = [interp1(1:length(g_180(:)),g_180,1:0.1994:length(g_180(:))),0]; % YJ 1mm / 3ms

% integral of G90
for ii = 1: length(g_90_intp)
    if(ii==1)
        g_90_integral(ii) = g_90_intp(end);
    else
        g_90_integral(ii) =  g_90_integral(ii-1)+g_90_intp(end-ii+1);
    end
end
g_90_integral = g_90_integral .* ts_rf ;
CenterPhase_90 = (g_90_integral(end))/2.0;
g_90_integral = -(g_90_integral(end:-1:1)-CenterPhase_90);

% integral of G180
for ii = 1: length(g_180_intp)
    if(ii==1)
        g_180_integral(ii) = g_180_intp(end);
    else
        g_180_integral(ii) =  g_180_integral(ii-1)+g_180_intp(end-ii+1);
    end
end
g_180_integral = g_180_integral .* ts_rf ;
CenterPhase_180 = (g_180_integral(end))/2.0;
g_180_integral = -(g_180_integral(end:-1:1)-CenterPhase_180);

% Nout=length(rfEncOut);
% phs_shift=repmat(exp(((0:1/Nout:1-1/Nout)-1/2+1/Nout/2)*2*pi*1i*freqShift).',1,G);
phs_shift = repmat(exp(1i*2*pi*(GAM/10) * g_90_integral* amplscale *slice_gap).',1,G);
% for gg = 1:G
%     phs_shift_(:,gg) = [interp1(1:length(phs_shift(:,gg)),phs_shift(:,gg),1:0.4999:length(phs_shift(:,gg))),0];
% end
% phs_shift = phs_shift_;
% clear phs_shift_
rfEncOutmb2=rfEncOut1+rfEncOut2.*phs_shift;

% NoutOther=length(rfOtherOut);
phs_shift = exp(1i*2*pi*(GAM/10) *g_180_integral* amplscale*slice_gap).';
% phs_shift_ = [interp1(1:length(phs_shift),phs_shift,1:0.4998:length(phs_shift)),0]';
% phs_shift = phs_shift_;
rfOtherOutmb2=rfOtherOut1+rfOtherOut2.*phs_shift;
% rfOtherOut = rfOtherOut.';




for EncodingCount = 1:G
    rf_90 = rfEncOutmb2(:,EncodingCount).';
    
    figure(6);
    subplot(G,1,EncodingCount);hold on;
    plot(real(rf_90));
    plot(imag(rf_90),'r');
    plot(abs(rf_90),'g');title('MB2 90 excitation')
    title('rf90 MB2');
    
    
    figure(7);
    plot(g_90);
    title('g90 MB2');
    legend('g90 MB2');
end

figure(8);
plot(real(rfOtherOutmb2));
plot(imag(rfOtherOutmb2),'r');
plot(abs(rfOtherOutmb2),'g');title('MB2 180 excitation')
title('rf180 MB2');

figure(9);
plot(g_180); 
title('g180 MB2');
legend('g180 MB2');



end