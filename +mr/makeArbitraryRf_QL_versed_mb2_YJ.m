function [rf, gz, gzReph, rf180, gz180] = makeArbitraryRf_QL_versed_mb2_YJ(flip, igSlider, varargin)
%makeArbitraryRf Create an RF pulse with the given pulse shape.
%   rf=makeArbitraryRf(singal, flip) Create RF pulse with complex signal 
%   and given flip angle (in radians)
%
%   rf=makeArbitraryRf(..., 'FreqOffset', f,'PhaseOffset',p)
%   Create block pulse with frequency offset and phase offset.
%
%   [rf, gz]=makeArbitraryRf(..., 'Bandwidth', bw, 'SliceThickness', st) 
%   Create RF pulse and corresponding slice select gradient. The bandwidth
%   of the pulse must be given for the specified shape.
%
%   See also  Sequence.makeSincPulse, Sequence.addBlock

validPulseUses = mr.getSupportedRfUse();

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeArbitraryRf';
    
    addRequired(parser, 'flipAngle', @isnumeric);
    addOptional(parser, 'system', mr.opts(), @isstruct);
    addParamValue(parser, 'freqOffset', 0, @isnumeric);
    addParamValue(parser, 'phaseOffset', 0, @isnumeric);
    addParamValue(parser, 'timeBwProduct', 0, @isnumeric);
    addParamValue(parser, 'bandwidth', 0, @isnumeric);
    % Slice params
    addParamValue(parser, 'maxGrad', 0, @isnumeric);
    addParamValue(parser, 'maxSlew', 0, @isnumeric);
    addParamValue(parser, 'sliceThickness', 0, @isnumeric);
    % Delay
    addParamValue(parser, 'delay', 0, @isnumeric);
    addParamValue(parser, 'dwell', 0, @isnumeric); % dummy default value
    % whether it is a refocusing pulse (for k-space calculation)
    addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));

end
parse(parser,flip,varargin{:});
opt = parser.Results;


if opt.dwell==0
    opt.dwell=opt.system.rfRasterTime;
end
%% from generateRFfile_CL_V1

load('1p1mm/rf_MB2_1p1mm_032624_s26.mat')


%% 90 pulse
signal=rf_90_set_wRefocus_MB2(2:end-1,:);
signal = signal./abs(sum(signal.*2*opt.dwell))*flip/(2*pi); % original!
% array_scale=[1.068 1.261 1.122 1.261 1.068];
% array_scale=[1 9/8 9/8 9/8 1];
array_scale=[1 1 1 1 1];
signal=signal*array_scale(igSlider);
signal=signal(:,igSlider);
rf.signal=signal(1:end);
dwell_1=2*opt.system.rfRasterTime;
rf.t=([1:length(rf.signal)]-0.5)*dwell_1;
rf.shape_dur=length(rf.signal)*dwell_1;
rf.delay = opt.delay;
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;
rf.deadTime = opt.system.rfDeadTime;
rf.ringdownTime = opt.system.rfRingdownTime;
rf.type = 'rf';
rf.delay = opt.delay;
if rf.deadTime > rf.delay
    rf.delay = rf.deadTime;
end
% g_90_versed_4us=g_90_versed_4us(1:end-1);% 1148
g_90_versed_4us=g_90_versed_4us(1:1148);% 1148 YJ use this
% g_90_versed_4us=g_90_versed_4us(1:329);% 1mm 3ms case
gz=mr.makeArbitraryGrad('z', g_90_versed_4us*425800,'delay',rf.delay,'system',opt.system); 
gz.first=0;
gz.last=0;
g_90_rewinder=g_90_versed_rew_4us(1149:end); % use this
% g_90_rewinder=g_90_versed_rew_4us(330:end); % 1mm 3ms case
g_90_rewinder(end) = 0; % YJ
gzReph=mr.makeArbitraryGrad('z', g_90_rewinder*425800,'system',opt.system);
gzReph.first=0;
gzReph.last=0;
%% 180 pulse
signal=rf_180_versed_MB2(2:end-1);
flip_180=2*flip;
signal = signal./abs(sum(signal.*2*opt.dwell))*flip_180/(2*pi); % 180 flip
% array_scale=[1.068 1.261 1.122 1.261 1.068];
array_scale_180=1/1.00;
signal=signal*array_scale_180;
rf180.signal=signal(1:end);
dwell_1=2*opt.system.rfRasterTime;
rf180.t=([1:length(rf180.signal)]-0.5)*dwell_1;
rf180.shape_dur=length(rf180.signal)*dwell_1;
rf180.delay = opt.delay;
rf180.freqOffset = opt.freqOffset;
rf180.phaseOffset = opt.phaseOffset;
rf180.deadTime = opt.system.rfDeadTime;
rf180.ringdownTime = opt.system.rfRingdownTime;
rf180.delay = opt.delay;
rf180.type = 'rf';
if rf180.deadTime > rf180.delay
    rf180.delay = rf180.deadTime;
end
g_180_versed_4us=g_180_versed_4us(1:end-1);
gz180=mr.makeArbitraryGrad('z', g_180_versed_4us*425800,'delay',rf180.delay,'system',opt.system); 
gz180.first=0;
gz180.last=0;

end
