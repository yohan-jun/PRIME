function [rf] = generate_rf_QL(use)
%GENERATE_RF_QL 

N = 256; % # time points in filter
G = 5; % gSlider factor
Gpulse = 'ex'; % 'ex' or 'se' gSlider encoding
tbG = 12%12; % overall tb product of encoding pulse; Should be > 2*G. If G/tbG
% is too high, you may get an error about non-increasing band edges from
% firls
tbOther =8%8; % tb product of non-encoding pulse
usecvx = false; % use Boyd's cvx toolbox for beta filter design
dt = 2e-3; % ms, final dwell time of pulses  default: 2e-3
T = 11; % ms, pulse duration of gSlider pulse; other pulse duration will be tbOther/tbG*T
slThick = 5.0; % mm, gSlider slice thickness
otherThickFactor = 1.03% factor to increase slice thickness of non-gSlider pulse
DFTphs = false; % do DFT phases
cancelAlphaPhs = true; % Design the excitation pulse's beta to cancel its associated alpha phase
doRootFlip = false; % root-flip the non-encoding pulse,

filename = 'RF_5x_linPhase'; % 
% filename = 'RF_5x_rootflip'; 

% and design the gSlider pulses to cancel the root-flipped phase profile.
% This requires that the non-encoding pulse have lower tb than the encoding
% pulse. Ideally, the non-encoding pulse should have 2x lower tb than the encoding pulse,
% otherwise there will be some increase ripple in the encoding pulse profile
if strcmp(Gpulse,'ex')
    bsf = sqrt(1/2); % excitation pulse
    d1 = 0.01;d2 = 0.01; % passband and stopband ripples of the overall profile
    d1 = sqrt(d1/2); % Mxy passband ripple
    d2 = d2/sqrt(2); % Mxy stopband ripple
    d1O = 0.01;d2O = 0.01; % passband and stopband ripples of the se profile
    phi = pi; % slice phase - e.g., for ex, will be pi, for se, will be pi/2
    cvx_osfact = 8;
elseif strcmp(Gpulse,'se')
    bsf = 1; % spin echo pulse
    d1 = 0.001;d2 = 0.01; % passband and stopband ripples of the overall profile
    d1 = d1/4;
    d2 = sqrt(d2);
    d1O = 0.01;d2O = 0.01; % passband and stopband ripples of the ex profile
    phi = pi/2; % slice phase
    cvx_osfact = 8;
end

% print out some info about what we are doing
fprintf('--------gSlider RF Pulse Design---------\n');
fprintf('Designing %s gSlider encoding pulses.\n',Gpulse);
fprintf('Number of sub-slices: %d\n',G);
fprintf('Slab time-bandwidth product: %g\n',tbG);
fprintf('Duration: %d ms\n',T);
fprintf('Slice Thickness: %g mm\n',slThick);
fprintf('Output dwell time: %g ms\n',dt);
fprintf('Time-bandwidth product of other pulse: %g\n',tbOther);
fprintf('Slice thickness ratio of other pulse: %g\n',otherThickFactor);

% design and simulate the other pulse
if strcmp(Gpulse,'ex')
    Gother = 'se';
else
    Gother = 'ex';
end
fprintf('Designing the %s (non-encoding) pulse\n',Gother);
[rfOther,bOther] = dzrf(N,tbOther,Gother,'ls',d1O,d2O);
% simulate pulse on a scaled grid that matches the encoding pulse
[apO,bpO] = abr(rfOther,(-N/2:1/8:N/2-1/8)*tbOther/tbG/otherThickFactor);
if strcmp(Gother,'se')
    MxyO = bpO.^2;
elseif strcmp(Gother,'ex')
    MxyO = 2*conj(apO).*bpO.*exp(1i*2*pi/N*N/2*(-N/2:1/8:N/2-1/8)'*tbOther/tbG/otherThickFactor);
end

rfEnc = zeros(N,G);
Mxy = zeros(N*8,G);
nomFlips = zeros(G,1);
for Gind = 1:G % sub-slice to design for

    fprintf('Designing pulse for sub-slice %d of %d\n',Gind,G);

    % design the beta filter
    if ~DFTphs
        if usecvx
            fprintf('Designing beta filter using cvx\n');
            if doRootFlip
                if strcmp(Gpulse,'ex')
                    phsFact = 2; % we have to square the 180 beta phase to
                    % cancel it with an EX pulse
                else
                    phsFact = 1/2; % we have to halve the 90 beta phase to
                    % cancel it with an SE pulse
                end
                b = bsf*gSliderBeta_cvx(N,G,Gind,tbG,d1,d2,phi,cvx_osfact,...
                    bOther,phsFact,tbOther/tbG/otherThickFactor);
            else
                b = bsf*gSliderBeta_cvx(N,G,Gind,tbG,d1,d2,phi,cvx_osfact);
            end
        else % use firls
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('Designing beta filter using firls\n');
            b = bsf*gSliderBeta(N,G,Gind,tbG,d1,d2,phi);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if doRootFlip
                if strcmp(Gpulse,'ex')
                    phsFact = 2; % we have to square the 180 beta phase to
                    % cancel it with an EX pulse
                else
                    phsFact = 1/2; % we have to halve the 90 beta phase to
                    % cancel it with an SE pulse
                end
                freqFactor = tbOther/tbG/otherThickFactor;
                BPhsMatch = exp(-1i*2*pi/N*(-N/2:N/2-1)'*freqFactor*...
                    (-N/2:N/2-1))*bOther(:);
                BPhsMatchOS = exp(-1i*2*pi/N*(-N/2:1/8:N/2-1/8)'*freqFactor*...
                    (-N/2:N/2-1))*bOther(:);
                %b = ifft(fft(b(:)).*ifftshift(exp(1i*angle(BPhsMatch(:).^phsFact))));
                target = fft(b(:),8*N).* ...
                    ifftshift(exp(1i*angle(BPhsMatchOS(:).^phsFact)));
                if strcmp(Gpulse,'ex') && cancelAlphaPhs
                    a = b2a(b);
                    target = target.*exp(1i*angle(fft(flipud(a(:)),8*N)));
                end
                b = lsqr(@(x,tflag)fftoversamp(x,N,8,tflag),target);
            end
        end
    else
        fprintf('Designing DFT beta filter using firls\n');
        phs = 2*pi/G*(ceil(-G/2):ceil(G/2)-1)*(Gind+ceil(-G/2)-1);
        if strcmp(Gpulse,'se'); phs = phs./2; end
        b = bsf*gSliderBetaDFT(N,phs,tbG,d1,d2);
    end

    % scale and solve for rf - note that b2rf alone doesn't work bc
    % b is not flipped correctly wrt a for non-symmetric profiles
    a = b2a(b);
    if ~doRootFlip && strcmp(Gpulse,'ex') && cancelAlphaPhs
        %b = ifft(fft(b(:).').*exp(1i*angle(fft(fliplr(a(:).')))));
        target = fft(b(:),8*N).*exp(1i*angle(fft(flipud(a(:)),8*N)));
        b = lsqr(@(x,tflag)fftoversamp(x,N,8,tflag),target);
    end
    rfEnc(:,Gind) = cabc2rf(a,fliplr(b(:).')); % iSLR for min-power pulse

    % simulate the pulse
    [ap,bp] = abr(rfEnc(:,Gind),-N/2:1/8:N/2-1/8);
    % calculate target flip angle of pulse in degrees
    nomFlips(Gind) = 2*asin(abs(bp(length(bp)/2+1)))*180/pi;
    if strcmp(Gpulse,'ex')
        Mxy(:,Gind) = 2*conj(ap).*bp.*exp(1i*2*pi/N*N/2*(-N/2:1/8:N/2-1/8)');
    elseif strcmp(Gpulse,'se')
        Mxy(:,Gind) = bp.^2;
    end

end


% plot all the profiles
zG = (-N/2:1/8:N/2-1/8)*slThick/tbG;

% interpolate pulses to target dwell time
Nout = T/dt;
rfEncOut = zeros(Nout,G);
for ii = 1:G
    rfEncOut(:,ii) = interp1((0:N)./N*T,[rfEnc(:,ii); rfEnc(end,ii)],(0:Nout-1)*dt,'spline',0);
    rfEncOut(:,ii) = rfEncOut(:,ii)./sum(abs(rfEncOut(:,ii)))*sum(abs(rfEnc(:,ii)));
end
Tother = T*tbOther/tbG/otherThickFactor;
   NoutOther = round(Tother/dt);
%   NoutOther = round(Tother/dt)-1; %QL especially for T 14ms, otherwises it will violate the time raster 
rfOtherOut = interp1((0:N)./N*Tother,[rfOther, rfOther(end)],(0:NoutOther-1)*dt,'spline',0);
rfOtherOut = rfOtherOut./sum(abs(rfOtherOut))*sum(abs(rfOther));

% convert to uT
% rfEncOut = rfEncOut./(2*pi*42.58*dt*10^-3);
% rfOtherOut = rfOtherOut./(2*pi*42.58*dt*10^-3);
%% QL
if use
    rf=rfEncOut;
else
    rf=rfOtherOut;
end
%% end of QL
end

