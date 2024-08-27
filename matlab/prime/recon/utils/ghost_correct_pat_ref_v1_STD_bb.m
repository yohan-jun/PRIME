function [PAT_ref_PATnavCor_kspace] = ghost_correct_pat_ref_v1_STD_bb(prot, PAT_ref,PAT_nav, AccY)

%% ________________________________________________________________________
%% Ghost correction for FLEET Ref Data
%% ________________________________________________________________________
%%
%% ________________________________________________________________________
%%

  
disp('Step:  Ghost Correction start')

sz_patrefscan = size(PAT_ref);
nPhases = sz_patrefscan(2);
sy = 1; % start line: 3

lin_fit_RefnavCor = mrir_artifact_ghost_compute_CorrelationMethod(PAT_nav);


Ref=zeros(size(PAT_ref));
Ref_RefnavCor=zeros(size(PAT_ref));
for R=0: AccY-1
    R
    Ref(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:)=PAT_ref(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:);

    Ref_RefnavCor(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:)= mrir_artifact_ghost_correct_CorrelationMethod_v3(Ref(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:), lin_fit_RefnavCor);
end

data_hybrid = mrir_iDFT_freqencode(Ref_RefnavCor); 

clear Ref_RefnavCor

if prot.alRegridMode==1    % if Regrid Mode =1, that means there is no trapezoid gridding.
    data_hybrid_tp=data_hybrid;
else
    data_hybrid_tp = mrir_regrid_trapezoid(data_hybrid, prot);
end

Ref_RefnavCor_kspace = mrir_fDFT_freqencode(data_hybrid_tp);

PAT_ref_PATnavCor_kspace=sum(Ref_RefnavCor_kspace,8);

disp('Step:  Ghost Correction complete')

end