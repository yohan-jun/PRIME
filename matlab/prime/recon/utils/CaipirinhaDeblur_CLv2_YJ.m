function [Kcorrected] = CaipirinhaDeblur_CLv2_YJ(K, prot, PhaseShiftBtwSimulSlices)

% Kawin Setsompop 
% Oct 19 2009

% now work for all acquisition orientations
% assume that the slices are already sorted (so can use para in original prot to calc stuff) 

% SliceSep in mm


% correction of SlicePos value due to rotation of principle gradient axis
% if size(K,10) > 1;
% 
%     % When the value for these coordinate is very small, e.g.
%     % sPosition_dSag = 1.343452e-9 then the read_meas_dat function that readin the data would not
%     % recognize it and will leave the array empty so fix it here
%     if isempty(prot.sSliceArray(1).sPosition_dSag) && isempty(prot.sSliceArray(2).sPosition_dSag)
%         for count = 1:length(prot.sSliceArray(1:end))
%             prot.sSliceArray(count).sPosition_dSag = 0;
%             prot.sSliceArray(count).sNormal_dSag = 0;
%         end
%     end
%     if isempty(prot.sSliceArray(1).sPosition_dCor) && isempty(prot.sSliceArray(2).sPosition_dCor)
%         for count = 1:length(prot.sSliceArray(1:end))
%             prot.sSliceArray(count).sPosition_dCor = 0;
%             prot.sSliceArray(count).sNormal_dCor = 0;
%         end
%     end
%     if isempty(prot.sSliceArray(1).sPosition_dTra) && isempty(prot.sSliceArray(2).sPosition_dTra)
%         for count = 1:length(prot.sSliceArray(1:end))
%             prot.sSliceArray(count).sPosition_dTra = 0;
%             prot.sSliceArray(count).sNormal_dTra = 0;
%         end
%     end
% 
%     NormalVec = [prot.sSliceArray(1).sNormal_dSag, prot.sSliceArray(1).sNormal_dCor, prot.sSliceArray(1).sNormal_dTra].';   
%     Pos(:,1) = [prot.sSliceArray(1:end).sPosition_dSag].';
%     Pos(:,2) = [prot.sSliceArray(1:end).sPosition_dCor].';
%     Pos(:,3) = [prot.sSliceArray(1:end).sPosition_dTra].';
%     SlicePos = Pos*NormalVec;% +4.3/2;
% else
%     keyboard('only single slice data so cant determine correctionFac if Gz is rotated')
% end


SlicePos = [size(K,10)-0.5:-1:0];
%SlicePos  = SlicePos(end:-1:1);

% PhaseShiftPerMM = (PhaseShiftBtwSimulSlices/SliceSep);
PhaseShiftPerMM = (PhaseShiftBtwSimulSlices/size(K,10));

Kcorrected = (zeros(size(K)));

if PhaseShiftBtwSimulSlices ~= 0
    for SlcCount = 1: size(K,10)
        if abs(PhaseShiftBtwSimulSlices) == pi % FOV/2 shift
%             Kcorrected(:,1:2:end,:,:,:,:,:,:,:,SlcCount) =  K(:,1:2:end,:,:,:,:,:,:,:,SlcCount);
%             Kcorrected(:,2:2:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:2:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*(SlicePos(SlcCount)));            
            Kcorrected(:,2:2:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:2:end,:,:,:,:,:,:,:,SlcCount);
            Kcorrected(:,1:2:end,:,:,:,:,:,:,:,SlcCount) =  K(:,1:2:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*(SlicePos(SlcCount)));
        elseif abs(PhaseShiftBtwSimulSlices) == 2*pi/3 % FOV/3 shift
            Kcorrected(:,1:3:end,:,:,:,:,:,:,:,SlcCount) =  K(:,1:3:end,:,:,:,:,:,:,:,SlcCount);
            Kcorrected(:,2:3:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:3:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,3:3:end,:,:,:,:,:,:,:,SlcCount) =  K(:,3:3:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*2*PhaseShiftPerMM*SlicePos(SlcCount));
        elseif abs(PhaseShiftBtwSimulSlices) == pi/2 % FOV/4 shift
            Kcorrected(:,1:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,1:4:end,:,:,:,:,:,:,:,SlcCount);
            Kcorrected(:,2:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:4:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,3:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,3:4:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*2*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,4:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,4:4:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*3*PhaseShiftPerMM*SlicePos(SlcCount));
        end
    end
else
    Kcorrected = K;
end



