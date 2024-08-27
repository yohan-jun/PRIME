function raw_corr = mrir_artifact_ghost_correct_CorrelationMethod_v3(raw_roft, linear_fit_coeff,OS_factor,ShiftFracForLine1)

% Kawin Setsompop
% 8/2/2012

% linear_fit_coeff here will be constant phase (row1) and secondlineshift
% (row2)

% input data needs to be in fourier space

% 8/16/2011
% much faster than version 2 and can do non-interger shift as the shift is
% done as a linear phase in the image domain

% 8/17/2011
% add the replication of linear_fit_coeff for cases with many reps


if nargin == 2
    OS_factor = 5; % oversampling factor
    ShiftFracForLine1 = 0.5; % the fraction of the shift that will be applied to the odd line group
end

if ( mrir_ice_dimensions(raw_roft, 'seg') < 2 ),
    error('uncorrected data contains only one segment');
end;

%% Take from Jon's old ghost correction code 

  % some sequences (e.g., "ge_functionals") collects reference lines for
  % phase correction every slice.


  % each data dimension will require its own unique phase correction lines,
  % except for the segments dimension (unless truly multi-shot or segmented
  % acquisition) and partitions dimension.

  %   dimension 08: Segments
  %   dimension 09: Partitions

  dims_raw = size(raw_roft);
  dims_raw(end+1:16) = 1;
  
  if (  ( mrir_ice_dimensions(raw_roft, 'rep') > 1 )  &&  ( prod(dims_raw([3:7,10:16])) ~= size(linear_fit_coeff,2) )  ),
      
      %%% SPECIAL CASE: repetitions
      linear_fit_coeff = reshape(linear_fit_coeff, [2 1 dims_raw(3:6) 1 1 1 dims_raw(10:16)]);
      %                                            1 2 3 4 5 6                                    7 8 9 0 1 2 3 4 5 6
      linear_fit_coeff = repmat(linear_fit_coeff, [1 1 1 1 1 1 mrir_ice_dimensions(raw_roft, 'rep') 1 1 1 1 1 1 1 1 1]);
      linear_fit_coeff = reshape(linear_fit_coeff, 2, []);
      
  end;
  
  
  % if user supplies a single slope and offset for entire data set,
  % replicate the coefficients over all dimensions except Seg and Par
  if ( size(linear_fit_coeff, 2) == 1 );
      linear_fit_coeff = repmat(linear_fit_coeff, [1, prod(dims_raw([3:7,10:16]))]);
  end;

%% 

raw_roft = double(raw_roft);

% Kawin Setsompop 9/07/2011

if sum(raw_roft(:,1,1,1,1,1,1,1,1,1,1,1),1) ~= 0
    fwd_lines = raw_roft(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:);
    rev_lines = raw_roft(:,2:2:end,:,:,:,:,:,2,:,:,:,:,:,:,:,:);
else
    fwd_lines = raw_roft(:,2:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:);
    rev_lines = raw_roft(:,1:2:end,:,:,:,:,:,2,:,:,:,:,:,:,:,:);
end

fwd_linesReshape = reshape(fwd_lines, size(fwd_lines,1), size(fwd_lines,2),[]);
fwd_linesReshapeCorrected = fwd_linesReshape;
rev_linesReshape = reshape(rev_lines, size(rev_lines,1), size(rev_lines,2),[]);
rev_linesReshapeCorrected = rev_linesReshape;


LineLength = size(raw_roft,1);
for DataSetCount = 1:size(fwd_linesReshape,3)
    Line1 = fwd_linesReshape(:,:,DataSetCount);
    Line2 = rev_linesReshape(:,:,DataSetCount);
    
    ConstPhaseDiff = linear_fit_coeff(1,DataSetCount);
    Total_SecondLineShiftBy = linear_fit_coeff(2,DataSetCount); 
    Line1Shift = -Total_SecondLineShiftBy*ShiftFracForLine1;
    Line2Shift = Total_SecondLineShiftBy + Line1Shift;
   
    Phase1 = exp(sqrt(-1)*2*pi* (-ceil((LineLength-1)/2):floor((LineLength-1) /2))/LineLength *(Line1Shift/OS_factor) ).';
    fwd_linesReshapeCorrected(:,:,DataSetCount) = mrir_fDFT(mrir_iDFT(Line1,1).*repmat(Phase1,1,size(Line1,2)),1)  /exp(sqrt(-1)*ConstPhaseDiff);
    Phase2 = exp(sqrt(-1)*2*pi* (-ceil((LineLength-1)/2):floor((LineLength-1) /2))/LineLength *(Line2Shift/OS_factor) ).';
    rev_linesReshapeCorrected(:,:,DataSetCount) = mrir_fDFT(mrir_iDFT(Line2,1).*repmat(Phase2,1,size(Line2,2)),1);
end

raw_corr = raw_roft;
if sum(raw_roft(:,1,1,1,1,1,1,1,1,1,1,1),1) ~= 0
    dims = size(raw_corr(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:));
    raw_corr(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:) = reshape(fwd_linesReshapeCorrected,dims);
    dims = size(raw_corr(:,2:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:));
    raw_corr(:,2:2:end,:,:,:,:,:,2,:,:,:,:,:,:,:,:) = reshape(rev_linesReshapeCorrected,dims);
else
    dims = size(raw_corr(:,2:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:));
    raw_corr(:,2:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:) = reshape(fwd_linesReshapeCorrected,dims);
    dims = size(raw_corr(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:));
    raw_corr(:,1:2:end,:,:,:,:,:,2,:,:,:,:,:,:,:,:) = reshape(rev_linesReshapeCorrected,dims);
end


if (0)
    raw_corrCollapsed = sum(raw_corr,8);
    raw_corrCollapsedUnCorrect = sum(raw_roft,8);
    Slc = 1;
    FigShift = 20; 
    
    % look at resulting images    
    figure(100+FigShift); 
    subplot(2,3,1); imagesc(mrir_array_combine(mrir_image_crop(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsedUnCorrect(:,:,:,1,1,1,1,1,1,Slc))))),0)); title('Uncorrected')
    subplot(2,3,2); imagesc(mrir_array_combine(mrir_image_crop(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsed(:,:,:,1,1,1,1,1,1,Slc))))),0)); title('Corrected')
    subplot(2,3,3); imagesc(mrir_array_combine(mrir_image_crop(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsed(:,1:2:end,:,1,1,1,1,1,1,Slc))))),0)); title('Odd line')
  
    for CoilCount =1:32
        figure(101+FigShift); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsedUnCorrect(:,:,CoilCount,1,1,1,1,1,1,Slc))))),[0 10]); title('Uncorrected')
        figure(102+FigShift); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsed(:,:,CoilCount,1,1,1,1,1,1,Slc))))),[0 10]); title('Corrected')
        figure(103+FigShift); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsedUnCorrect(:,1:2:end,CoilCount,1,1,1,1,1,1,Slc))))),[0 10]); title('Oddlines');
        %figure(104); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsed(:,2:2:end,CoilCount,1,1,1,1,1,1,Slc))))),[0 50]);        
        %figure(104); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsedUnCorrect(:,:,CoilCount,1,1,1,1,1,1,Slc))))),[0 50]);

    end
    
    % look at k-space data
    
    CoilToLookAt = [7, 9, 23,27]; 
    centerKxline = size(raw_corrCollapsed,1) *(1/2)   +1; 
    centerKyline = size(raw_corrCollapsed,2) *(1/3)   +1; % 6/8 PF
    windowKx = -ceil(size(raw_corrCollapsed,1) *(1/8)):ceil(size(raw_corrCollapsed,1) *(1/8)) ;
    windowKy = -ceil(size(raw_corrCollapsed,2) *(1/8)):ceil(size(raw_corrCollapsed,2) *(1/8)) ;
    
    for CoilCount = 1:length(CoilToLookAt)
        figure(201+FigShift); subplot(2,3,CoilCount); imagesc((abs(raw_corrCollapsed(windowKx+centerKxline,windowKy+centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc))));
        %figure(202); subplot(2,3,CoilCount); imagesc((abs(raw_corrCollapsedUnCorrect(windowKx+centerKxline,windowKy+centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc))));
        
        figure(203+FigShift); subplot(2,3,CoilCount); imagesc(angle(raw_corrCollapsed(windowKx+centerKxline,windowKy+centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));
        figure(204+FigShift); subplot(2,3,CoilCount); imagesc(angle(raw_corrCollapsedUnCorrect(windowKx+centerKxline,windowKy+centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));
                
%         figure(205); subplot(2,3,CoilCount); plot(abs(raw_corrCollapsed(windowKx+centerKxline,centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));
%         hold on; subplot(2,3,CoilCount); plot(abs(raw_corrCollapsed(windowKx+centerKxline,centerKyline+1,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)),'r'); 
%         
%         figure(206); subplot(2,3,CoilCount); plot(angle(raw_corrCollapsed(windowKx+centerKxline,centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));
%         hold on; subplot(2,3,CoilCount); plot(angle(raw_corrCollapsed(windowKx+centerKxline,centerKyline+1,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)),'r'); 
%                 
%          figure(207); subplot(2,3,CoilCount); plot(abs(raw_corrCollapsed(centerKxline,:,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));       
%          figure(208); subplot(2,3,CoilCount); plot(angle(raw_corrCollapsed(centerKxline,:,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));

    end
     
end


  