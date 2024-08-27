%
% 	function [bv] = verse(b,g,gv)
%
%	MEX function version of VERSE algorithm.
%
%	INPUT:
%		b 	= RF pulse (real or complex, arb units)
%		g 	= gradient pulse (arb units)
%		gv	= VERSE gradient pulse (same units as g)
%
%	OUTPUT:
%		bv	= VERSE RF pulse.
%
%	NOTES:
%		1) All waveforms are assumed to have the same time scale.
%		2) This is just a resampling of the time function.  


function [bv] = verse(b,g,gv)

bv = verse_mex(b,g,gv);		% Just call the MEX function if this
				% function is invoked.




