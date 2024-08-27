% Kawin Setsompop
% 15 Aug 2009

% Hack function to use the new minTimeGradient (Lustig) instead of mintimegrad code (Brian)
% the old code is "mintimegradOriginal.m"

% IMPORTNAT NOTE!!!!!

% This new gradient calculation function cannot do a couple of
% things that the original code can do:
% 1) it cannot be use to rampdown the gradient (in which case we dont care
% about the moment (i.e. kspace) term and we just want to ramp down as
% quickly as possible. This is the case when we want to finish up the
% excitation spiral trajectory
% 2) it cannot be use in the case of 2D or 3D trajectory where g0 and gf are non-zero and different for e.g. gx and gy, 
% as the code will only take in one g0 and gf values
%
% For these cases replace "mintimegrad" with "mintimegradOriginal"


function [g] = mintimegrad(Nest,g0,gf,moment,T,gmax,smax,t0,type)

sprintf('Using new minTimeGradient Code')

if size(moment) == 1
    C = linspace(0,moment, 256); C = C+i*eps;
else
    if sum(abs(g0) + abs(gf)) ~= 0
        sprintf('Lustig minTimeGrad: can only take in single g0 and gf values')
        keyboard
    end
    if size(moment) == 2
        C = linspace(0,moment(1), 256) + i*linspace(0,moment(2), 256);
    else
        C = [linspace(0,moment(1),256).', linspace(0,moment(2),256).', linspace(0,moment(3),256).'] ;
    end
end
[k,time,g,s] = minTimeGradient(C,g0(1), gf(1), gmax, 0.98*smax/1000,T*10^3);
