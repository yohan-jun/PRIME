function b = mrir_zeropad(a, padsize, direction)
%MRIR_ZEROPAD pad an array with zeros
%
% B = MRIR_ZEROPAD(A,PADSIZE) pads array A with PADSIZE(k) number of zeros
% along the k-th dimension of A. PADSIZE should be a vector of positive
% integers.
%
% B = MRIR_ZEROPAD(A,PADSIZE,DIRECTION) pads A in the direction specified by
% the string DIRECTION. DIRECTION can be one of the following strings.
%
%       string values for DIRECTION
%       'pre'         pads before the first array element along each
%                     dimension .
%       'post'        pads after the last array element along each
%                     dimension.
%       'both'        pads before the first array element and after the
%                     last array element along each dimension.

% (this function is a subset of matlab's PADARRAY, which is a part of the
% Image Processing toolbox---and therefore eats a precious license when its
% called.)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2009/dec/04
% $Id: mrir_zeropad.m,v 1.2 2011/03/28 04:14:47 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % preprocess the padding size
  if ( numel(padsize) < ndims(a) ),
    padsize           = padsize(:);
    padsize(ndims(a)) = 0;
  end;

  numDims = numel(padsize);

  % form index vectors to subsasgn input array into output array.
  % also compute the size of the output array.
  idx   = cell(1,numDims);
  sizeB = zeros(1,numDims);
  for k = 1:numDims,
    M = size(a,k);
    switch direction,
     case 'pre',
      idx{k}   = (1:M) + padsize(k);
      sizeB(k) = M + padsize(k);

     case 'post',
      idx{k}   = 1:M;
      sizeB(k) = M + padsize(k);

     case 'both',
      idx{k}   = (1:M) + padsize(k);
      sizeB(k) = M + 2*padsize(k);
    end;
  end;

  % initialize output array with the padding value and make sure the
  % output array is the same type as the input
  b         = cast( repmat(0, sizeB), class(a) );
  b(idx{:}) = a;

  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_zeropad.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
