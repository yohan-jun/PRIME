%
% 	function [bv,gv] = mintverse(b,g,dt,bmax,gmax,smax,dtout)
%
%	MEX function version of minimum-time VERSE algorithm.
%
%	INPUT:
%		b 	= RF pulse (real or complex, arb units)
%		g 	= gradient pulse (arb units)
%		dt	= time steps of each "block" (arb units)
%		bmax 	= maximum RF amplitude, in units of b.
%		gmax 	= maximum Gradient amplitude in units of g.
%		smax 	= maximum Gradient slew rate in units of g per unit dt
%		dtout   = output sampling step. (units of dt).
%		emax	= maximum RF energy (units if b*b*dt).  
%				(optional, no constraint if omitted)
%
%	OUTPUT:
%		bv	= minimum-time VERSE RF pulse.
%		gv	= minimum-time gradient waveform.
%


function [bv,gv] = mintverse(varargin)

% 	Arguments past dt are not required - mex function
%	will set defaults.  Thus we use the "varargin" clause.
%
%	Note that calling mintverse and mintverse_mex should be
%	equivalent.
%

[bv,gv] = mintverse_mex(varargin{:});



