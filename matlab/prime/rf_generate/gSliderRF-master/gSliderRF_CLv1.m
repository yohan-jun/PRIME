close all;
% Script to design gSlider pulses.
addpath rf_tools/ % JP's tools: gets dinf, b2a, cabc2rf, abr...
addpath functions/ % 
N = 256; % # time points in filter
G = 5; % gSlider factor
Gpulse = 'ex'; % 'ex' or 'se' gSlider encoding
tbG = 12%12; % overall tb product of encoding pulse; Should be > 2*G. If G/tbG
% is too high, you may get an error about non-increasing band edges from
% firls
tbOther =8%8; % tb product of non-encoding pulse
usecvx = false; % use Boyd's cvx toolbox for beta filter design
dt = 2e-3; % ms, final dwell time of pulses
T = 11; % ms, pulse duration of gSlider pulse; other pulse duration will be tbOther/tbG*T 11
slThick = 5.0; % mm, gSlider slice thickness
otherThickFactor = 1.03% factor to increase slice thickness of non-gSlider pulse
DFTphs = false; % do DFT phases
cancelAlphaPhs = true; % Design the excitation pulse's beta to cancel its associated alpha phase
doRootFlip = false; % root-flip the non-encoding pulse,

filename = 'RF_5x_linPhase_1p03_5mm_rf_2us'; % 
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
if doRootFlip
    fprintf('Root-flipping the non-encoding pulse\n');
    if tbOther > tbG
        error 'Non-encoding tb must be < gSlider tb to root-flip'
    end
%     [rfOther,bOther] = rootFlip(bOther,d1O,pi/2+strcmp(Gother,'se')*pi/2,tbOther,0.6);
     [rfOther,bOther] = rootFlip(bOther,d1O,pi/2+strcmp(Gother,'se')*pi/2,tbOther);   % Jun Ma's modification

end
if strcmp(Gpulse,'se') && cancelAlphaPhs % this pulse is the ex pulse; cancel alpha phs
    rfOther = b2rf(ifft(fft(bOther(:).').*exp(-1i*angle(fft(fliplr(b2a(bOther)))))));
end
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
h1 = figure;
h2 = figure;
for ii = 1:G
    if strcmp(Gpulse,'ex')
        titleText = sprintf('EX profile; gSlider factor %d; sub-slice %d',G,ii);
        seSig = conj(Mxy(:,ii)).*MxyO;
    else
        titleText = sprintf('SE profile; gSlider factor %d; sub-slice %d',G,ii);
        seSig = Mxy(:,ii).*conj(MxyO);
    end

    % plot the first encoding pulse
    figure(h1);
    subplot(G*100 + 10 + ii),hold on
    plot(zG,abs(Mxy(:,ii)));
    plot(zG,real(Mxy(:,ii)));
    plot(zG,imag(Mxy(:,ii)));
    title(titleText);
    legend('|Mxy|','Mx','My');
    xlabel 'mm'
    axis([min(zG) max(zG) -1 1]);

    % plot the overall spin echo signal
    figure(h2);
    titleText = sprintf('Spin echo signal profile; gSlider factor %d; sub-slice %d',G,ii);
    subplot(G*100 + 10 + ii),hold on
    plot(zG,abs(seSig));
    plot(zG,real(seSig));
    plot(zG,imag(seSig));
    title(titleText);
    legend('|SE signal|','Re\{SE Signal\}','Im\{SE Signal\}');
    xlabel 'mm'
    axis([min(zG) max(zG) -1 1]);

end

% plot both the pulse profiles
if strcmp(Gpulse,'ex')
    plind = 1;
    plindO = 2;
    titleText = sprintf('EX profile; gSlider factor = %d; sub-slice 1',G);
    titleTextO = 'SE profile';
else
    plind = 2;
    plindO = 1;
    titleText = sprintf('SE pulse; gSlider factor = %d; sub-slice 1',G);
    titleTextO = 'EX profile';
end
% plot the first encoding pulse
figure;subplot(310+plind),hold on
z = (-N/2:1/8:N/2-1/8)*slThick/tbG;
plot(z,abs(Mxy(:,1)));
plot(z,real(Mxy(:,1)));
plot(z,imag(Mxy(:,1)));
title(titleText);
legend('|Mxy|','Mx','My');
xlabel 'mm'
axis([min(z) max(z) -1 1]);

% plot the other pulse
subplot(310+plindO),hold on
plot(z,abs(MxyO));
plot(z,real(MxyO));
plot(z,imag(MxyO));
title(titleTextO);
legend('|Mxy|','Mx','My');
xlabel 'mm'
axis([min(z) max(z) -1 1]);

% plot them together, to compare thicknesses
subplot(313),hold on
plot(z,abs(Mxy(:,1)));
plot(z,abs(MxyO));
legend('|Mxy|','|Mxy|, other');
xlabel 'mm'
axis([-2*slThick 2*slThick 0 1]);
title 'Both profile amplitudes'

% interpolate pulses to target dwell time
Nout = T/dt;
rfEncOut = zeros(Nout,G);
for ii = 1:G
    rfEncOut(:,ii) = interp1((0:N)./N*T,[rfEnc(:,ii); rfEnc(end,ii)],(0:Nout-1)*dt,'spline',0);
    rfEncOut(:,ii) = rfEncOut(:,ii)./sum(abs(rfEncOut(:,ii)))*sum(abs(rfEnc(:,ii)));
end
Tother = T*tbOther/tbG/otherThickFactor;
NoutOther = round(Tother/dt);
rfOtherOut = interp1((0:N)./N*Tother,[rfOther, rfOther(end)],(0:NoutOther-1)*dt,'spline',0);
rfOtherOut = rfOtherOut./sum(abs(rfOtherOut))*sum(abs(rfOther));

% convert to uT
rfEncOut = rfEncOut./(2*pi*42.58*dt*10^-3);
rfOtherOut = rfOtherOut./(2*pi*42.58*dt*10^-3);

gAmp = tbG/(T*10^-3)/(slThick*10^-3)/42580; % mT/m, gradient amplitude

% simulate the final slice profiles, after interpolation
z = linspace(-2,2,8*N)*slThick; % mm
MxyOut = zeros(8*N,G);
for ii = 1:G
    [ap,bp] = abr(2*pi*42.58*dt*10^-3*rfEncOut(:,ii),...
        2*pi*42580*dt*10^-3*gAmp*ones(Nout,1),z./1000);
    if strcmp(Gpulse,'ex')
        MxyOut(:,ii) = 2*conj(ap).*bp.*exp(1i*2*pi*42580*dt*10^-3*gAmp*Nout/2*z(:)./1000); 
    elseif strcmp(Gpulse,'se')
        MxyOut(:,ii) = bp.^2;
    end
end
[ap,bp] = abr(2*pi*42.58*dt*10^-3*rfOtherOut,...
    2*pi*42580*dt*10^-3*gAmp*ones(NoutOther,1),z./1000);
if strcmp(Gpulse,'ex')
    MxyOut = bsxfun(@times,conj(MxyOut),bp.^2);%2*conj(ap).*bp.*exp(1i*2*pi/N*N/2*(-N/2:1/8:N/2-1/8)');
elseif strcmp(Gpulse,'se')
    MxyOut = bsxfun(@times,conj(2*conj(ap).*bp.*exp(1i*2*pi*42580*dt*10^-3*gAmp*Nout/2*z(:)./1000)),MxyOut);
end
% calculate signal matrix
edges = zeros(G,2);
ftw = dinf(d1,d2)/tbG; % fractional transition width of the slab profile
for ii = 1:G
    Gcent = (ii-G/2-1/2)*slThick/G;        
    edges(ii,:) = [Gcent-(slThick/G/2-ftw*(slThick/2)) Gcent+(slThick/G/2-ftw*(slThick/2))];
end
encMtx = zeros(G);
for ii = 1:G % loop over pulses
    for jj = 1:G % loop over subslices
        encMtx(ii,jj) = sum(MxyOut(z > edges(jj,1) & z < edges(jj,2),ii));
    end
end

% plot the spin echo signal profiles of the output pulses, and 
h1 = figure;
for ii = 1:G

    % plot the overall spin echo signal
    figure(h1);
    titleText = sprintf('Spin echo signal profile of interpolated output pulse; sub-slice %d',ii);
    subplot(G*100 + 10 + ii),hold on
    plot(z,abs(MxyOut(:,ii)));
    plot(z,real(MxyOut(:,ii)));
    plot(z,imag(MxyOut(:,ii)));
    % also plot the band edges we used for encMtx calculation
    for jj = 1:G
        plot(edges(jj,1)*ones(10,1),-1:2/9:1,'--k');
        plot(edges(jj,2)*ones(10,1),-1:2/9:1,'--k');
    end
    title(titleText);
    legend('|SE signal|','Re\{SE Signal\}','Im\{SE Signal\}');
    xlabel 'mm'
    axis([min(z) max(z) -1 1]);

end
    
% plot the final interpolated output pulses
tRFEnc = (0:size(rfEncOut,1)-1)*dt;
tRFOther = (0:length(rfOtherOut)-1)*dt;
figure
negBnd = floor(min([real(rfEncOut(:));imag(rfEncOut(:))]));
posBnd = ceil(max(abs(rfEncOut(:))));
for ii = 1:G
  subplot(2,G,ii),hold on
  plot(tRFEnc,real(rfEncOut(:,ii)));
  plot(tRFEnc,imag(rfEncOut(:,ii)));
  plot(tRFEnc,abs(rfEncOut(:,ii)));
  title(sprintf('%s gSlider Pulse %d',Gpulse,ii));
  xlabel 'ms'
  ylabel 'uT'
  axis([0 max(tRFEnc) negBnd posBnd])
end

subplot(210+plindO),hold on
plot(tRFOther,real(rfOtherOut));
plot(tRFOther,imag(rfOtherOut));
plot(tRFOther,abs(rfOtherOut));
if strcmp(Gpulse,'ex')
    title(sprintf('SE pulse'));
else
    title(sprintf('EX pulse'));
end
xlabel 'ms'
ylabel 'uT'
legend('|RF|','Re\{RF\}','Im\{RF\}');

% write out the pulses, gradient amplitude, target slice thickness,
% expected flip at middle of slice
% TODO: Should we phase shift the 90 or 180 by pi/2 for CPMG?
scaleFactEnc = zeros(G,1);
foo = rfOtherOut./max(abs(rfOtherOut));
scaleFactOther = 100.96/(90+strcmp(Gother,'se')*90)*...
        (sqrt(sum(real(foo))^2 + sum(imag(foo))^2)/386.3952)...
        *(max(abs(rfOtherOut))/4.2616); 
for ii = 3:3
    foo = rfEncOut(:,ii)./max(abs(rfEncOut(:,ii)));
    scaleFactEnc(ii) = 100.96/(90+strcmp(Gother,'ex')*90)*...
        (sqrt(sum(real(foo))^2 + sum(imag(foo))^2)/386.3952)...
        *(max(abs(rfEncOut(:,ii)))/4.2616); 
end
%%
% scaling factors for Connectom scanner. This is referenced to the middle
% (purely real) pulse of a gSlider 5 linear phase pulse set. The first
% ratio is target flip of that pulse, divided by the flip the scanner wants
% to apply. The second ratio is the AmpInt of the current pulse, divided by
% the AmpInt of the reference pulse, the third ratio is the peak amplitude of
% the current pulse, divided by the peak amplitude of the reference pulse.
% if doRootFlip == true;suffix = '_rootFlipped';else suffix = 'linPhs';end;
%  if DFTphs == true;suffix2 = '_DFT';else suffix2 = '';end
%  save(sprintf('gSliderRF_%dx_%s%s%s',G,Gpulse,suffix,suffix2),'rfEncOut','rfOtherOut',...
%      'gAmp','slThick','dt','nomFlips','scaleFactEnc','scaleFactOther');
 
 if doRootFlip == true;suffix = '_rootFlipped';else suffix = 'linPhs';end;
if DFTphs == true
    suffix2 = '_DFT';
else
    suffix2 = '';
end
zOut=z./slThick;

% save(['./Profiles/',filename,'_MxyOut.mat'],'MxyOut','zOut');
save([filename,'_MxyOut.mat'],'MxyOut','zOut');
% save(sprintf(['./Profiles/',filename,'_gSliderRF_%dx_%s%s%s'],G,Gpulse,suffix,suffix2),'rfEncOut','rfOtherOut',...
%     'gAmp','slThick','dt','nomFlips','scaleFactEnc','scaleFactOther');
save(sprintf(['/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/',filename,'_gSliderRF_1129_%dx_%s%s%s'],G,Gpulse,suffix,suffix2),'rfEncOut','rfOtherOut',...
    'gAmp','slThick','dt','nomFlips','scaleFactEnc','scaleFactOther');

size(rfEncOut)
size(rfOtherOut)

% T*1000/2
% Tother*1000/2
