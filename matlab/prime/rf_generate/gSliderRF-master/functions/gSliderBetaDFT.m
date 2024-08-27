function b = gSliderBetaDFT(N,phi,tb,d1,d2)

% Script to design gSlider beta filters.
ftw = dinf(d1,d2)/tb; % fractional transition width of the slab profile
G = length(phi); % # sub-slices

%% Design a gSlider beta polynomial - can get negative
% band positions by flipping the waveform
if max(abs(phi)) == 0
    % No phase, simple uniform pulse
    f = [0 (1-ftw)*(tb/2) (1+ftw)*(tb/2) (N/2)]/(N/2);
    m = [1 1 0 0];
    w = [1 d1/d2];
    b = firls(N-1,f,m,w); % the filter
else
    % Each sub-slice has its own phase;
    % design two shifted filters that we can add to kill off the left band,
    % then demodulate the result back to DC
    shift = N/4;
    f = [0 shift-(1+ftw)*(tb/2)]; % left stopband
    for ii = 1:length(phi)
        Gcent = shift+(ii-G/2-1/2)*tb/G;
        f = [f Gcent-(tb/G/2-ftw*(tb/2)) Gcent+(tb/G/2-ftw*(tb/2))];
    end
    f = [f shift+(1+ftw)*(tb/2) (N/2)]./(N/2); % right stopband
    m = [0 0 kron(exp(1i*phi(:).'),[1 1]) 0 0];
    w = [d1/d2 ones(1,G) d1/d2];
    b1 = firls(N-1,f,m,w); % the filter
    b2 = firls(N-1,f,m,w,'h');
    b = (b1+1i*b2).*exp(-1i*2*pi/N*shift*(0:N-1))/2*exp(-1i*pi/N*shift);
end
