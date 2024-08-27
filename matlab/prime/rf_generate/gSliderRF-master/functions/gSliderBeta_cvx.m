function b = gSliderBeta_cvx(N,G,Gind,tb,d1,d2,phi,osfact,bPhsMatch,phsFact,freqFactor)

% Script to design a gSlider beta filter using cvx.
ftw = dinf(d1,d2)/tb; % fractional transition width of the slab profile

% set up target patterns and weights
nn = (-N/2*osfact:N/2*osfact-1)'/(N/2*osfact);  % -1 to 1 - profile indices

Gcent = (Gind-G/2-1/2)*tb/G/(N/2);
if Gind > 1 && Gind < G
    % separate transition bands for slab+slice
    % here we check that each band has at least 1 point in it;
    % if not, we add 1 point in the middle. 
    lbnd = -(1-ftw)*(tb/2)/(N/2);
    rbnd = Gcent-(tb/G/2+ftw*(tb/2))/(N/2);
    if max(double(nn >= lbnd & nn <= rbnd)) == 0
        nn = sort([nn; lbnd + (rbnd-lbnd)/2],'ascend');
    end
    lbnd = Gcent-(tb/G/2-ftw*(tb/2))/(N/2);
    rbnd = Gcent+(tb/G/2-ftw*(tb/2))/(N/2);
    if max(double(nn >= lbnd & nn <= rbnd)) == 0
        nn = sort([nn; lbnd + (rbnd-lbnd)/2],'ascend');
    end
    lbnd = Gcent+(tb/G/2+ftw*(tb/2))/(N/2);
    rbnd = (1-ftw)*(tb/2)/(N/2);
    if max(double(nn >= lbnd & nn <= rbnd)) == 0
        nn = sort([nn; lbnd + (rbnd-lbnd)/2],'ascend');
    end
    d = double(nn >= -(1-ftw)*(tb/2)/(N/2) & nn <= Gcent-(tb/G/2+ftw*(tb/2))/(N/2));
    d = d + exp(1i*phi)*double(nn >= Gcent-(tb/G/2-ftw*(tb/2))/(N/2) & nn <= Gcent+(tb/G/2-ftw*(tb/2))/(N/2));
    d = d + double(nn >= Gcent+(tb/G/2+ftw*(tb/2))/(N/2) & nn <= (1-ftw)*(tb/2)/(N/2));
elseif Gind == 1
    % slab and sub-slice share a left transition band
    lbnd = -(1-ftw)*(tb/2)/(N/2);
    rbnd = Gcent+(tb/G/2-ftw*(tb/2))/(N/2);
    if max(double(nn >= lbnd & nn <= rbnd)) == 0
        nn = sort([nn; lbnd + (rbnd-lbnd)/2],'ascend');
    end
    lbnd = Gcent+(tb/G/2+ftw*(tb/2))/(N/2);
    rbnd = (1-ftw)*(tb/2)/(N/2);
    if max(double(nn >= lbnd & nn <= rbnd)) == 0
        nn = sort([nn; lbnd + (rbnd-lbnd)/2],'ascend');
    end
    d = exp(1i*phi)*double(nn >= -(1-ftw)*(tb/2)/(N/2) & nn <= Gcent+(tb/G/2-ftw*(tb/2))/(N/2));
    d = d + double(nn >= Gcent+(tb/G/2+ftw*(tb/2))/(N/2) & nn <= (1-ftw)*(tb/2)/(N/2));
elseif Gind == G
    % slab and sub-slice share a right transition band
    lbnd = -(1-ftw)*(tb/2)/(N/2);
    rbnd = Gcent-(tb/G/2+ftw*(tb/2))/(N/2);
    if max(double(nn >= lbnd & nn <= rbnd)) == 0
        nn = sort([nn; lbnd + (rbnd-lbnd)/2],'ascend');
    end
    lbnd = Gcent-(tb/G/2-ftw*(tb/2))/(N/2);
    rbnd = (1-ftw)*(tb/2)/(N/2);
    if max(double(nn >= lbnd & nn <= rbnd)) == 0
        nn = sort([nn; lbnd + (rbnd-lbnd)/2],'ascend');
    end    
    d = double(nn >= -(1-ftw)*(tb/2)/(N/2) & nn <= Gcent-(tb/G/2+ftw*(tb/2))/(N/2));
    d = d + exp(1i*phi)*double(nn >= Gcent-(tb/G/2-ftw*(tb/2))/(N/2) & nn <= (1-ftw)*(tb/2)/(N/2));
end
s = nn <= -(1+ftw)*(tb/2)/(N/2) | nn >= (1+ftw)*(tb/2)/(N/2); % stopband mask
wts = 1./abs(nn).^0.5; % decaying ripple weights
wts = wts.*s;wts(isnan(wts)) = 0;wts = wts./max(wts);

if exist('bPhsMatch','var')
    % evaluate the passed-in beta filter at the frequencies in f,
    % to incorporate the same phase into the current filter, so that the
    % phase cancels if it is
    BPhsMatch = exp(-1i*2*pi/N*nn*(N/2)*freqFactor*(-N/2:N/2-1))*bPhsMatch(:);
    d = d.*exp(1i*angle(BPhsMatch(:).^phsFact));
end

% masked desired pattern
dd = d(s | logical(abs(d)));
% masked stopband error pattern
wtsd = wts(s | logical(abs(d)));

ddri = [real(dd);imag(dd)];
ddriabs = [abs(dd);abs(dd)];
wtsdri = [wtsd;wtsd];
A = exp(-1i*2*pi/N*nn*N/2*(-N/2:N/2-1));
Ad = A(s | logical(abs(d)),:);
At = A(~s & ~logical(abs(d)),:);
Adri = [real(Ad) -imag(Ad);imag(Ad) real(Ad)]; % split to real+imaginary for cvx
Atri = [real(At) -imag(At);imag(At) real(At)];

% use cvx to do the constrained optimization
cvx_begin
    variable delta(1)
    variable x(2*N)
    minimize( delta )
    subject to
        -delta*ddriabs <= Adri*x - ddri <= delta*ddriabs + delta*d2/d1*wtsdri
        -1+delta <= Atri*x <= 1+delta
cvx_end

% get complex filter back
b = x(1:N)+1i*x(N+1:end);

