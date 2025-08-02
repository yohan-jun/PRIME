function S_hat = mrir_array_whitening_apply(S, W)
%MRIR_ARRAY_WHITENING_APPLY
%
% S_hat = mrir_array_whitening_apply(S, W)
%
% S_hat = S * W.'
%
% see also MRIR_ARRAY_WHITENING_OPERATOR, MRIR_ARRAY_TRANSFORM.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/aug/03
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  dims = size(S);
  dims(end+1:16) = 1;

  Ncha = dims(3);

  % final size is given by number of rows in W
  dims(3) = size(W, 1);

  if ( size(W,1) ~= Ncha || size(W,2) ~= Ncha ),
    error('size of whitening operator does not match data dimensions');
  end;
  
  Nrem = prod(dims([1:2,4:end]));

  % permute and reshape n-dimension array to Ncha x Nrem matrix
  s = reshape(permute(S, [3, 1, 2, 4:16]), Ncha, []);

  if ( (size(s,1) ~= Ncha) || (size(s,2) ~= Nrem) ),
    error('incorrect dimensions found in resized matrix');
  end;

  s_hat = (s.' * W.').';

  S_hat = ipermute(reshape(s_hat, dims([3, 1, 2, 4:16])), [3, 1, 2, 4:16]);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_transform.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
