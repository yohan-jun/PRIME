function [ConstPhaseDiff, SecondLineShiftBy, Line1_corrected] = mrir_artifact_ghost_compute_CorrelationMethod_subfunction(Line1,Line2,OS_factor)

% Kawin Setsompop
% 7/29/2012

if nargin == 2
    OS_factor = 5; % oversampling factor
end

% convention here is to fix Line1

Line1 = double(Line1); 
Line2 = double(Line2);


%% STEP 1: need to pad the data to get higher resolution to get aacurate calculation for the
% shift that needs to be corrected 
LineLength = size(Line1,1);
PadOnEachSide = floor(LineLength*(OS_factor-1)/2);

Line1_OS = mrir_fDFT(mrir_zeropad( mrir_iDFT(Line1,1),PadOnEachSide,'both'),1);
Line2_OS = mrir_fDFT(mrir_zeropad( mrir_iDFT(Line2,1),PadOnEachSide,'both'),1);

%Line1_OS = mrir_fDFT(padarray( mrir_iDFT(Line1,1),PadOnEachSide),1);
%Line2_OS = mrir_fDFT(padarray( mrir_iDFT(Line2,1),PadOnEachSide),1);

% might also want to chop just so look at only center portion??? (maybe for speed)?

%figure(1) 
%subplot(2,1,1); plot(abs(Line1_OS)); hold on; plot(abs(Line2_OS),'r')
%subplot(2,1,2); plot(abs(xcorr(Line1_OS,Line2_OS)))


%% STEP 2: calculate the shift based on cross-correlation (this is better
% than using conv as it does not wrap and only use points that exist on
% both shifted dataset for a given shift)

[value,shift] = max(abs(  conv(Line1_OS,conj(Line2_OS(end:-1:1))) ));
%[value,shift] = max(abs(xcorr(Line1_OS,Line2_OS)));
SecondLineShiftBy = shift-size(Line1_OS,1);


%% STEP 3: fix the shift (might also want to do conjugation to the shifted point)
Line1_OSshifted = circshift(Line1_OS,-SecondLineShiftBy);

if(0)
    % use conjugate symetry approximation here.....
    if SecondLineShiftBy<0
        Line1_OSshifted(1:SecondLineShiftBy) = conj(Line1_OSshifted(1:SecondLineShiftBy));
    else
        Line1_OSshifted(end-SecondLineShiftBy:end) = conj(Line1_OSshifted(end-SecondLineShiftBy:end));
    end
end

%figure; plot(abs(Line1_OSshifted)); hold on; plot(abs(Line2_OS),'r');plot(abs(Line1_OS),'g'); 

%% STEP 4: calculate the constant phase diff between the lines and fix it 
PhaseDiff = angle(Line1_OSshifted./Line2_OS); 
IntensityMask = abs(Line1_OSshifted) > max(abs(Line1_OSshifted))* 0.1 ;
PhaseDiffSubset = PhaseDiff(IntensityMask); 

ConstPhaseDiff = lscov(ones(size(PhaseDiffSubset)),PhaseDiffSubset,abs(Line1_OSshifted(IntensityMask))); % weighted fit where weighting is based on abs(Line1_shifted)

Line1_OScorrected = Line1_OSshifted/exp(sqrt(-1)*ConstPhaseDiff);

%% STEP 5: down sample the final data by 3

%Line1_corrected = 3*Line1_OScorrected(1:OS_factor:end); % work only for even number of points along the line

Image_Line1_OScorrected = mrir_iDFT(Line1_OScorrected,1);
Line1_corrected = mrir_fDFT(Image_Line1_OScorrected(1+PadOnEachSide:LineLength+PadOnEachSide),1);

if (0)
%     figure
%     plot(abs(Line1_OScorrected)); hold on; plot(abs(Line2_OS),'r')
    
    figure;
    subplot(3,1,1); plot(abs(Line1_corrected)); hold on; plot(abs(Line2),'r'); title('Magnitude aligning'); 

    subplot(3,1,2); plot(PhaseDiff); hold on; plot(IntensityMask,'r'); title('PhaseDiff and IntensityMask'); 

    subplot(3,1,3);  plot(PhaseDiffSubset);
    hold on; plot(ConstPhaseDiff*ones(size(PhaseDiffSubset)),'r'); title('PhaseDiffSubset and fitted constant'); 
    
end





