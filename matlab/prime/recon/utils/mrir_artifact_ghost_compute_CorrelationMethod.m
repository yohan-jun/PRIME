function linear_fit_coeff = mrir_artifact_ghost_compute_CorrelationMethod(phascor,TakeMedian,OS_factor)

% Kawin Setsompop
% 8/2/2012

% linear_fit_coeff here will be constant phase (row1) and secondlineshift
% (row2)

% input data needs to be in fourier space

if nargin == 1
    TakeMedian = 0; %2
    OS_factor = 5;
end

dims = size(phascor);
if length(dims) < 10
    dims(end+1:10) = 1;
end

% ceil is for backward compatibility (if NSeg is 1)
dims(mrir_DIM_SEG) = ceil( dims(mrir_DIM_SEG)/2 );


% separate out forward and reverse reference lines (assuming polarity
% reverses every line)

% NOTE: convention is that segment "0" is a normal line and segment "1"
% is a reflected line (but it doesn't matter which is which as long as
% normal and reflected are in different segments)
%                          1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6

if sum(phascor(:,1,1,1,1,1,1,1,1,1,1,1),1) ~= 0
    ref_FWD = mean(phascor(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:), 2);  % line "0" is forward
    ref_REV = mean(phascor(:,2:2:end,:,:,:,:,:,2,:,:,:,:,:,:,:,:), 2);  % line "0" is forward
else
    ref_FWD = mean(phascor(:,2:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:), 2);  % line "0" is forward
    ref_REV = mean(phascor(:,1:2:end,:,:,:,:,:,2,:,:,:,:,:,:,:,:), 2);  % line "0" is forward
end


% collapse two sets of reference lines each into an Nreadout x Nlines array
ref_REV = reshape(ref_REV, [dims(1), prod(dims(3:end))]);
ref_FWD = reshape(ref_FWD, [dims(1), prod(dims(3:end))]);
linear_fit_coeff = zeros(2,size(ref_REV,2));

for DataSetCount = 1:size(ref_REV,2)
    DataSetCount;
    [linear_fit_coeff(1,DataSetCount), linear_fit_coeff(2,DataSetCount), Line1_corrected{DataSetCount}] = ...
        mrir_artifact_ghost_compute_CorrelationMethod_subfunction(ref_FWD(:,DataSetCount),ref_REV(:,DataSetCount),OS_factor);
end
         

if TakeMedian == 1 % median
    linear_fit_coeff = reshape( repmat(median( reshape(linear_fit_coeff, [2 1 dims(3:7) 1 1 dims(10)]), 7) , [1 1 1 1 1 1 dims(7) 1 1 1] ) , 2,[]) ;
elseif TakeMedian == 2 % moving median
    
    TR_window = 1:10;
    
    coeff_matrix = reshape(linear_fit_coeff, [2 1 dims(3:7) 1 1 dims(10)]);
    coeff_matrix2 = coeff_matrix;
    
    CurrentGroupIndex = TR_window;
    if dims(7) > length(TR_window)
        for count = 1:dims(7)
            shift = count - ceil(length(TR_window)/2);
            if (shift > 0) && (CurrentGroupIndex(end) < dims(7))
                CurrentGroupIndex = TR_window + shift;
            end
            coeff_matrix_currentGroup = reshape(coeff_matrix(:,:,:,:,:,:, CurrentGroupIndex ,:,:,:),[2 1 dims(3:6) length(TR_window) 1 1 dims(10)]);
            coeff_matrix2(:,:,:,:,:,:,count,:,:,:) = median(coeff_matrix_currentGroup, 7);
        end
        linear_fit_coeff = reshape( coeff_matrix2, 2,[]) ;
        clear coeff_matrix2 coeff_matrix coeff_matrix_currentGroup
    else
        linear_fit_coeff = reshape( repmat(median( reshape(linear_fit_coeff, [2 1 dims(3:7) 1 1 dims(10)]), 7) , [1 1 1 1 1 1 dims(7) 1 1 1] ) , 2,[]) ;
    end
elseif TakeMedian == 3 % take median over all the slices and time
    b = repmat(median( repmat(median( reshape(linear_fit_coeff, [2 1 dims(3:7) 1 1 dims(10)]), 7) , [1 1 1 1 1 1 dims(7) 1 1 1]),10), [1 1 1 1 1 1 1 1 1 dims(10)]);
    linear_fit_coeff = reshape(b, 2,[]) ; clear b
end    

if(0)
    linear_fit_coeff_reshape = reshape(linear_fit_coeff,2,dims(3),dims(7),dims(10));
    figure(111); clf; figure(112); clf; 
    %Slice = 8;
    Rep = 1;
    Ncoils = dims(3);
    for count = 1:Ncoils
        figure(111); hold on; plot(squeeze(linear_fit_coeff_reshape(1,count,Rep,:))); title('constant phase term')
        figure(112); hold on; plot(squeeze(linear_fit_coeff_reshape(2,count,Rep,:))); title('pixel Shift (with OS9)')
    end
end
    
%%%%%%%%%
  
