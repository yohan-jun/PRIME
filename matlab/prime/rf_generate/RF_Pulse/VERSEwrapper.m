function [rf_versed, g_versed] = VERSEwrapper(rf,g,gmax,smax,Ts,VERSEdownFac)

addpath(genpath('Library_RFpulse'))

%gmax = 3;                     % max gradient amplitude [G/cm]
%smax = 15000;                 % maximum slew rate [G/cm/s]
%Ts = 10*10^-6/oversample;     %sampling time [s]

    %NEED TO PAD RF AND GRADIENT!!!! Otherwise verse function do not work!!!
    rf = [ 0; 0 ; rf ; 0; 0];
    rf_max = VERSEdownFac*max(abs(rf));
    g = [ 0; 0 ; g ; 0; 0];

    
    [rf_versed,g_versed] = mintverse(rf,g,Ts,rf_max,gmax,smax,Ts);
    
    
    if(0)       
        sprintf('Hack for Brian VERSE code: make last 6 rf values the same as start (avoid occasional peaky error)')
        rf_versed(end-5:end) = rf_versed(6:-1:1);
        g_versed(end-5:end) = g_versed(6:-1:1);    
    end
    
    PulseDurationBeforeVERSE = (length(g))*Ts/10^-3; %in m
    PulseDurationAfterVERSE = (length(g_versed))*Ts/10^-3; %in m

end