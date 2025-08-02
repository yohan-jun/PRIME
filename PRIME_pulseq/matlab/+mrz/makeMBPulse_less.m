function [pulseout,p] = makeMBPulse_less(pulsein,varargin)
%makeMBPulse make an SMS pulse inspired from sigpy.mri.rf.multiband.mb_rf
%added from mb_rf tool:
%nBands = 2 different slice position for debugging
%nBands = 0 for arbitrary array of distances and phases
%
%       mb_rf(pulse_in, n_bands=3, band_sep=5, phs_0_pt='None')
%          Args:
%              pulse_in (array): samples of single-band RF pulse.
%              n_bands (int): number of bands.
%              band_sep (float): normalized slice separation.
%              phs_0_pt (string): set of phases to use. Can be 'phs_mod' (Wong),
%                 'amp_mod' (Malik), 'quad_mod' (Grissom), or 'None'
% 
%          band_sep = slice_sep/slice_thick*tb, where tb is time-bandwidth product
%          of the single-band pulse
% 
%          Returns:
%              array: multibanded pulse out
%              p:     modulation functions to check the maximum
% 
%          References:
%              Wong, E. (2012). 'Optimized Phase Schedules for Minimizing Peak RF
%              Power in Simultaneous Multi-Slice RF Excitation Pulses'. Proc. Intl.
%              Soc. Mag. Reson. Med., 20 p. 2209.
%              Malik, S. J., Price, A. N., and Hajnal, J. V. (2015). 'Optimized
%              Amplitude Modulated Multi-Band RF pulses'. Proc. Intl. Soc. Mag.
%              Reson. Med., 23 p. 2398.
%
%       2022.03.18

validPulseUses = mr.getSupportedRfUse();

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeMBPulse';
    
    % RF params
    addRequired(parser, 'pulsein', @isnumeric);
    addOptional(parser, 'system', mr.opts(), @isstruct);
    addParamValue(parser, 'timeBwProduct', 4, @isnumeric);   
    addParameter(parser,'phs_0_pt','None',@(x) any(validatestring(x,{'None','phs_mod','amp_mod','quad_mod'})));
    addParameter(parser,'n_bands',3,@isnumeric);
    addParameter(parser,'band_sep',5,@isnumeric);
    addParameter(parser,'phs',[0],@isnumeric);
    % Slice params
    addParamValue(parser, 'maxGrad', 0, @isnumeric);
    addParamValue(parser, 'maxSlew', 0, @isnumeric);
    addParamValue(parser, 'sliceThickness', 0, @isnumeric);
    addParamValue(parser, 'delay', 0, @isnumeric);
    addParamValue(parser, 'dwell', 0, @isnumeric); % dummy default value    
    % whether it is a refocusing pulse (for k-space calculation)
    addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));
end

parse(parser, pulsein, varargin{:});
opt = parser.Results;

if opt.dwell==0
    opt.dwell=opt.system.rfRasterTime;
end
nBands  = opt.n_bands;
tb      = opt.timeBwProduct;
BandSep = opt.band_sep;
method  = opt.phs_0_pt;
if opt.phs == 0
    if nBands>0
        phs = zeros(nBands,1);
    else
        phs = zeros(length(BandSep),1);
    end
else
    phs = opt.phs;
    if nBands>0
        if length(phs)~=nBands
            error('phase array should be equal to nBands')
        end
    else
        if length(phs)~=length(BandSep)
            error('phase array should be equal to nBands')
        end
    end
end


n = length(pulsein);
p = 0;
if nBands==2
    for ii = 1:nBands
        p = p + exp(1j*2*pi*BandSep*tb*(-n/2:n/2-1)/n*(ii-1)+1j*phs(ii));
    end
elseif nBands == 0
%     BandSep = sort(BandSep);
%     phs = sort(phs);
    for ii = 1:length(BandSep)
        p = p + exp(1j*2*pi*BandSep(ii)*tb*(-n/2:n/2-1)/n+1j*phs(ii));
    end    
else
    for ii = 1:nBands
        p = p + exp(1j*2*pi*BandSep*tb*(-n/2:n/2-1)/n*(ii-(nBands+1)/2)+1j*phs(ii));%not in the center, exp(1j*2*pi*BandSep*tb*(-n/2+0.5:n/2-0.5)/n*(ii-(nBands+1)/2)+1j*phs(ii));
    end
end

pulseout = p.*pulsein;
