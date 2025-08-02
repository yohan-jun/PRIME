function data_deapodized = mrir_regrid_trapezoid_rolloff(data_roft, trapezoid)
%MRIR_REGRID_TRAPEZOID_ROLLOFF
%
% data_deapodized = mrir_regrid_trapezoid_rolloff(data_roft, trapezoid)
%
% deapodization


%  "ice_trapezoid_rolloff.m"

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/03
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%



  %// Find the length of the vectors.
  lLenX = size(data_roft, 1);
  if ( lLenX ~= trapezoid.NxRegrid ),
    error('length of line for regridding does not match regrid function!');
  end;

  denom = abs(trapezoid.regrid_rolloff).^2;

  
  data_deapodized = data_roft(:,:).*conj(repmat(trapezoid.regrid_rolloff,[1,size(data_roft(:,:),2)]))./(repmat(denom,[1,size(data_roft(:,:),2)]));


  data_deapodized = reshape(data_deapodized, size(data_roft));

  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
