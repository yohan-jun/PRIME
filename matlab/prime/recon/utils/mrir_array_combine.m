function img_combine = mrir_array_combine(img_multichan, METHOD, varargin)
%MRIR_ARRAY_COMBINE
%
% img_combine = mrir_array_combine(img_multichan, METHOD);

% TODO: return measure of phase compatible with root-sum-of-squares
% combination, e.g., homodyne detection (Noll et al., 1991)

% TODO: implement all special cases of optimal SNR combination, including
% normalizations.

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/27
% $Id: mrir_array_combine.m,v 1.2 2007/05/03 20:43:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % method 0: root-sum-of-squares
  % method 1: noise-weighted root-sum-of-squares
  % method 3: B1-weighted sum, normalization: uniform noise
  % method 4: B1-weighted sum, normalization: uniform signal
  % method 5: optimal, normalization: uniform noise
  % method 6: optimal, normalization: uniform signal
  % method 9: adaptive combine
  

  if ( size(img_multichan, 3) == 1 ),
    img_combine = img_multichan;
  end;
  
  
  switch METHOD,
   case 0,
    img_combine = mrir_array_combine_rss(img_multichan);
   otherwise,
    error('unknown combination method: "%d"', METHOD);
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_combine.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: