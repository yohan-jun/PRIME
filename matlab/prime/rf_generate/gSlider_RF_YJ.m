% From Congyu Liao, PhD
clc;close all;clear all;

lims = mr.opts('MaxGrad',78,'GradUnit','mT/m',...
    'MaxSlew',200,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', 2.89); % prisma
thickness=5*220/200*1e-3;            % slice thinckness
% for igSlider=1:5
% nomFlips=[90 90 90 90 90]'; % tell the scanner we want 90 degree pulse Qiang Liu test 1
% 
% [rf(igSlider), gz,rf180,gz180] = mr.makeArbitraryRf_QL_versed((pi*nomFlips(igSlider)./180), igSlider, 'system',lims,'timeBwProduct',12, ...
% 'SliceThickness',thickness);
% 
% end

% multiband2
for igSlider=1:5
nomFlips=[90 90 90 90 90]'; 

[rf(igSlider), gz,gzReph, rf180,gz180] = mr.makeArbitraryRf_QL_versed_mb2_YJ((pi*nomFlips(igSlider)./180), igSlider, 'system',lims,'timeBwProduct',12, ...
'SliceThickness',thickness);

end

save('1p1mm/seq_blocks_MB2_1p1mm_032624_s26.mat','rf','rf180','gz','gz180','gzReph');

% singleband 4 us
% for igSlider=1:5
% nomFlips=[90 90 90 90 90]'; % tell the scanner we want 90 degree pulse Qiang Liu test 1
% 
% [rf(igSlider), gz, gzReph, rf180,gz180] = mr.makeArbitraryRf_QL_versed_4us((pi*nomFlips(igSlider)./180), igSlider, 'system',lims,'timeBwProduct',12, ...
% 'SliceThickness',thickness);
% 
% end