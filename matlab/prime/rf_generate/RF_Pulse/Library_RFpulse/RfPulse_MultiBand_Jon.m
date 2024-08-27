% Kawin Setsompop Jan 09
% clean up 4th Aug 09
% clean up 20th Sep 09

function RfPulse_MultiBand

%close all

%% parameters
BlochSim = 1;

TBW  = 8;                           % time-bandwidth-product for the RF-sinc-pulse along z
type = 1;                           % 1 = 90� excitation, 2 = 180� spin-echo, 3 = 180� inversion
HF = 0;                             % HF design on/off
factor1 = i;                        % flip direction for first 180 (for diffusion)   
factor2 = -i;                       % flip direction for second 180 (for diffusion)

gmax = 0.5;%1;                      % max gradient amplitude [G/cm]
gmax_RefocusLobe = 3;
smax = 15000;                       % maximum slew rate [G/cm/s]
oversample = 5;                     % oversampling for G- and RF-wavefrom (on a 10us base) -- final sampling have to be divisible by 0.5us

zt = 0.7;                           % thickness of the slice, +-zt/2[cm]
SliceShift = 0;%4                   % slice-shift in z [cm]

Slice_Seperation = 6;               % seperation btw slices in z [cm]
position = 0;
%position = [-1 0 1];                    % positions of the slices (unit in Slice_Seperation)
%position = [-1/2 1/2];

VERSE = 1;                          % turn VERSE on (VERSE =1)
VERSEdownFac = 1/2;                 % How much to VERSE down the max RF amplitude
AmpCompensateForSideSlice = 0;      % compensate for side slices profile drop (due to rect sampling)

B0 = [0 100];                  % use to sim effect of B0 offres (in Hz)
%B0 = [75]; 
delay = 0; % in term of sampling points
SimulateFlip = 90;% degree %note: this is for simulation only, the actual flip stored will be e.g. 90 for excitation pulse

if rem(10/oversample,0.5) ~= 0
    disp('*******************************')
    disp('Sampling period not compatible with system!!!')
    disp('*******************************')
    keyboard
end
gamma_s = 4258;   % gamma^*[Hz/G]
Ts = 1*10^-5;     % base sampling time [s]

%% gradient calculation

% calculate the gradient
kzt = TBW/zt;       % the total kz extend that is traversed by each spoke in z-direction (total -kzt...+kzt)

C1 = linspace(0,-kzt/2, 256); C1 = C1+i*eps;
C2 = linspace(0,kzt, 256); C2 = C2+i*eps;
[k_temp1,time1,gz1,s_temp1] = minTimeGradient(C1, 0, 0, gmax_RefocusLobe, 0.95*smax/1000,Ts*10^3);
[k_temp2,time2,gz2,s_temp2] = minTimeGradient(C2, 0, 0, gmax, 0.95*smax/1000,Ts*10^3);

g_low  = [gz1;gz2]; % string them together
g_k = ([0; g_low] + [g_low; 0])/2;
k = cumsum(g_k)*gamma_s*Ts; k = k(1:end-1);k = k(end:-1:1,:);
g_low  = g_low(end:-1:1);

g_used = g_low(1:(length(gz2)));    % only this part does slice selection (RF on)
k_used = k(1:(length(gz2))); 

if type ~= 1; %i.e. no rewinder
    g_low = g_low(1:(length(gz2)));
    k     = k_used;
end

%% RF calculation
if type == 1;     a = 'ex';  flip = pi/2; 
elseif type == 2; a = 'se'; flip = pi;
elseif type == 3; a = 'inv';  flip = pi; 
end 

if HF == 1 % use Pauly's toolbox
    rf_ideal = dzrf(500,TBW,a,'ls');
    k_ideal = -(((1:500)-251)/250)*kzt/2;
    if type == 1;
        rf = interp1(k_ideal,rf_ideal,k_used,'spline');
        rf = rf.* abs(g_used);
        rf_sinc = [rf; zeros(length(gz1),1)];
    else
        rf = interp1(k_ideal,rf_ideal,k_used,'spline');
        rf_sinc = rf.*abs(g_low);
    end
else % low flip design
    rf_sinc = z_rf(kzt,k_used,g_used,TBW);
    if type == 1;
        rf_sinc = [rf_sinc; zeros(length(gz1),1)];
    end
end

% normalization of RF
rf_sum_ideal = flip/(2*pi*gamma_s*Ts); 
rf_sinc = rf_sinc* (rf_sum_ideal/sum(real(rf_sinc))); %Have to do "sum(real(rf_sinc))" as when rf goes neg flip is neg !!

%PulseDurationBeforeVERSE = (length(g_low))*Ts/10^-3 %in ms
PulseDurationBeforeVERSEwoRefoc = (length(gz2))*Ts/10^-3 %in ms

if VERSE == 1
    if (sum(g_used < 0) ~= 0)
        g_used = abs(g_used); % make sure G is not negative otherwise not work...
        sprintf('WARNING: G to VERSE has neg value- took "the abs" so not crash')
        %keyboard
    end
    
    %NEED TO PAD RF AND GRADIENT!!!! Otherwise verse function do not work!!!
    RfToVerse = [ 0; 0 ; rf_sinc(1:(length(gz2))) ; 0; 0];
    GzToVerse = [ 0; 0 ; g_used ; 0; 0];
    rf_max = VERSEdownFac*max(abs(rf_sinc));

    [rf_sincVersed,gz2_Versed] = mintverse(RfToVerse,GzToVerse,Ts,rf_max,gmax,smax,Ts);

    if(1)
        sprintf('Hack for Brian VERSE code: make last 3 rf values the same as start (avoid occasional peaky error)')
        rf_sincVersed(end-2:end) = rf_sincVersed(3:-1:1);
        gz2_Versed(end-2:end) = gz2_Versed(3:-1:1);
    end
    
    if type == 1;
        rf_sinc = [rf_sincVersed ;rf_sinc(length(gz2)+1:end)];
        g_low = [gz2_Versed ;g_low(length(gz2)+1:end)];
    else
        rf_sinc = rf_sincVersed;
        g_low = gz2_Versed;
    end
    gz2 = gz2_Versed;
    %PulseDurationAfterVERSE = (length(g_low))*Ts/10^-3 %in ms
    PulseDurationAfterVERSEwoRefoc = (length(gz2_Versed))*Ts/10^-3 %in ms
end

% implement oversampling & calculate k

g_high = reshape(repmat(g_low,[1,oversample]).',[],1);
g_k = ([0; g_high] + [g_high; 0])/2;
k = (cumsum(g_k(end:-1:1))*gamma_s*Ts/oversample);
k = k(1:end-1); k = k(end:-1:1,:);

if oversample ~= 1
    rf_sinc = Resample(rf_sinc,oversample,'spline'); 
    
    % method below also work but more noisy from noise in g (when
    % accounting rf for velocity)
    % % %     if type == 1
    % % %         rf_sinc = z_rf(kzt,k(1:length(gz2)*oversample),g_high(1:length(gz2)*oversample),TBW);
    % % %         rf_sinc = [rf_sinc; zeros(length(gz1)*oversample,1)];   % add the plateau
    % % %         rf_sum_ideal = flip/(2*pi*gamma_s*Ts/oversample); % normalization
    % % %         rf_sinc = rf_sinc* (rf_sum_ideal/sum(real(rf_sinc))); %Have to do "sum(real(rf_sinc))" as when rf goes neg flip is neg !!
    % % %     else
    % % %         rf_sinc = z_rf(kzt,k,g_high,TBW);
    % % %         rf_sum_ideal = flip/(2*pi*gamma_s*Ts/oversample); % normalization
    % % %         rf_sinc = rf_sinc* (rf_sum_ideal/sum(real(rf_sinc))); %Have to do "sum(real(rf_sinc))" as when rf goes neg flip is neg !!
    % % %     end
end

if type ~= 1; % not the same k-space definition for e.g. 180 
    k = k-kzt/2;
end   
  
% phase modulation for slice shift & outerslice compensation due to rect sampling
Phase_SliceShift = exp(i*2*pi*k*SliceShift);   % slice-shift in z
MultiSlice_Modulation = zeros(size(k));
for count = 1:length(position)
    rf_current = rf_sinc.*Phase_SliceShift.*exp(i*2*pi*k*Slice_Seperation*position(count));

    if AmpCompensateForSideSlice == 1
        disp('correct for RF envelope due to rect sampling')
        disp('Work for LF: Might not work perfectly for SLR design')
        LFscaling = 0.05/flip; %scale down rf to LF regime (0.05 radian)
        [angle1, angle2] = abr3D(real(rf_current)*LFscaling,imag(rf_current)*LFscaling,zeros(size(g_high)),zeros(size(g_high)),...
            g_high,0,0,SliceShift+Slice_Seperation*position(count),Ts/oversample*1e3*length(rf_current));
        CorrectionScaling = 0.05/abs(2*conj(angle1).*angle2)
        rf_indiv(:,count) = rf_current*CorrectionScaling;
    else
        rf_indiv(:,count) = rf_current;
    end
    
end
rf_sinc = sum(rf_indiv,2);

%% Plot results and Bloch Sim

% plot rf and gradient
figure(2);clf;
hold on; 
subplot(2,2,1); 
if HF == 1
    plot(abs(rf_sinc),'r');
else
    plot(abs(rf_sinc),'b');
    %plot((rf_sinc),'b');
end
grid on;
title('\bf{RF-excitation}','Color','k')

subplot(2,2,2);
plot(g_high,'r'); leg{1}=sprintf('g(t)'); hold on;
plot(k,'g');      leg{2}=sprintf('k(t)'); hold off;
legend(leg); grid on;
title('\bf{g- and k-trajectory}','Color','k')


rf_sinc = [zeros(delay,1); rf_sinc];
g_high = [g_high; zeros(delay,1)];

if BlochSim == 1
    %Bloch simulation
    Scaling = (SimulateFlip*pi/180)/flip;
    Bx = real(rf_sinc)*Scaling; By = imag(rf_sinc)*Scaling;
    gz = g_high; gy = zeros(size(gz)); gx = gy;
    x = 0; y = 0;
    z = -10:0.01:10;  % magnetization along z [cm]
    %z = -2:0.01:2;
    Bmain = 10; T2map = 100000;
   
    for B0_Count =1:length(B0)
        B0_Hz = B0(B0_Count);
        [angle1, angle2] = abr3D_complete(Bx,By,gx,gy,gz,x,y,z,Ts/oversample*1e3*length(Bx),B0_Hz,Bmain,T2map);
        
        if type ==1;      slice = 2*conj(angle1).*angle2;                               % paper Eq. 5
        elseif type == 2; slice = (angle2.*angle2)*i;                                   % paper Eq. 10
        elseif type == 3; slice = (angle1.*conj(angle1) - angle2.*conj(angle2))  -1 ;   % paper Eq. 6
        end
        subplot(2,1,2);hold on;
        if HF == 1
            plot(z,squeeze(real(slice)),'r--'); leg{1}=sprintf('real'); 
            plot(z,squeeze(imag(slice)),'r:');  leg{2}=sprintf('imag');
            plot(z,squeeze(abs(slice)),'r');    leg{3}=sprintf('abs');
            %legend(leg,'Location','NorthWest'); grid on;
            grid on;
            %title('\bf{m_xy (sinc slice selection HF)}','Color','r')
        else
            plot(z,squeeze(real(slice)),'k--'); leg{1}=sprintf('real');
            plot(z,squeeze(imag(slice)),'k:');  leg{2}=sprintf('imag');
            plot(z,squeeze(abs(slice)),'k');    leg{3}=sprintf('abs');
            grid on;
            %legend(leg,'Location','NorthWest'); grid on;
            %title('\bf{m_xy (sinc slice selection)}','Color','k')
        end
    end
   
    idealMxy_90 = zeros(size(z));
    for count = 1:length(position)
        currentpoint = position(count)*Slice_Seperation+SliceShift;
        idealMxy_current_90 = (z >= currentpoint - zt/2 -eps) .*(z <= currentpoint + zt/2 +eps);
        idealMxy_90 = idealMxy_90 + idealMxy_current_90;
    end
    plot(z,idealMxy_90*sin(flip*Scaling),'r--'); hold off;
    
end
ylim([-0.5 1.2])

rf_sinc = [ rf_sinc ; zeros(delay,1)];
g_high = [zeros(delay,1) ;g_high];
max_rf = max(abs(rf_sinc));

factor_for_overlaped = max(abs(rf_sinc))*4.258*2 % factor use to mult V_adj to get the voltage needed
factor_for_indiv1 = max(abs(rf_indiv(:,1)))*4.258*2 

maxRFadjVoltage = 1000/factor_for_overlaped

% AmpInt_overlaped = (pi/2)/( 2*pi*4258 * max(abs(rf_sinc)) * (10^-5)/oversample )
% AmpInt_indiv1 = (pi/2)/( 2*pi*4258 * max(abs(rf_indiv(:,1))) * (10^-5)/oversample )
% 
% SARampInt_Assume_1Vadj =  sum(( factor_for_overlaped*(abs(rf_sinc)/max_rf) ).^2)/oversample
% 
% format long
%     Normalized_AmpInt = sum(real(rf_sinc)./max(abs(rf_sinc)))* (10^-5)/oversample
% format short

%% Create Header 

fid = fopen('RFinfo.txt', 'wb');

fprintf(fid,'maxRFadjVoltage =  %g (assume 1000v max allow)\n', maxRFadjVoltage);
fprintf(fid,'factor_for_overlaped =  %g\n', factor_for_overlaped);
fprintf(fid,'factor_for_indiv1 =  %g\n', factor_for_indiv1);

fprintf(fid,'zt =  %g cm\n', zt);
fprintf(fid,'SliceShift =  %g cm\n', SliceShift);
fprintf(fid,'Slice_Seperation =  %g cm\n', Slice_Seperation);
fprintf(fid,'position =  %g cm\n', position);
fprintf(fid,'\n');

fprintf(fid,'VERSE =  %g\n', VERSE);
fprintf(fid,'VERSEdownFac =  %g\n', VERSEdownFac);

fprintf(fid,'TBW =  %g\n', TBW);

fprintf(fid,'gmax =  %g G/cm \n', gmax);
fprintf(fid,'smax =  %g G/cm/s \n', smax);
fprintf(fid,'oversample =  %g\n', oversample);
%fprintf(fid,'SARampInt_Assume_1Vadj =  %g\n', SARampInt_Assume_1Vadj);
if VERSE == 1
    fprintf(fid,'PulseDurationAfterVERSEwoRefoc = %g\n', PulseDurationAfterVERSEwoRefoc);
else
    fprintf(fid,'PulseDurationNotVERSEwoRefoc = %g\n', PulseDurationBeforeVERSEwoRefoc);
end
fprintf(fid,'AmpCompensateForSideSlice =  %g\n', AmpCompensateForSideSlice);

fclose(fid);

% total duration need to be multiple of 20 us (siemens constrain)
leftover = rem(length(rf_sinc),oversample*2);
if leftover ~= 0
    pad = oversample*2-leftover;
    rf_sinc = [rf_sinc; zeros(pad,1)];
    rf_indiv = [rf_indiv; zeros(pad,size(rf_indiv,2))];
    g_low = [g_low;0];
end

%if length(rf_sinc) > 2024*4
if length(rf_sinc) > 2024*2
    disp('RF pulse is longer than what the system allow')
    keyboard
end

if type == 1
    cd header
    create_header(rf_sinc,[zeros(length(g_low),2),g_low])
    for count = 1:length(position)
        eval([ 'cd ' num2str(count)])
        create_header(rf_indiv(:,count),[zeros(length(g_low),2),g_low])
        cd ..
    end
else 
    cd header_180   
    create_header_180_diffusion(rf_sinc,[zeros(length(g_low),2),g_low],factor1,factor2)
    for count = 1:length(position)
        eval([ 'cd ' num2str(count)])
        create_header_180_diffusion(rf_indiv(:,count),[zeros(length(g_low),2),g_low],factor1,factor2)
        cd ..
    end
end
cd ..



%% subfunction: z_rf()

function rf = z_rf(kzt,kz,gz,TBW)
% used to create the "sinc" RF-pulse along z, filtered by a cosine function

z = kz/(kzt/2);                                     % normalize
m = TBW/4;                                          % TBW defines number of zero-crossings
snc = sin(m*2*pi*z+0.00001)./(m*2*pi*z+0.00001);    % the sinc
rf = snc.*(0.54+0.46*cos(pi*z));                    % hanning window
rf = rf .* abs(gz);                                 % weighting due to the varying speed of k(t)

%% subfunction: resample()

function rf_high = Resample(rf,oversample,type)
% resample g with a factor [1,2,3,...]

x = 0:(length(rf))+1;
xi = 0.5+(1/oversample)/2:(1/oversample):length(rf)+0.5-(1/oversample)/2;
y = [0; rf; 0];
rf_high = interp1(x,y,xi,type).';



