function [W, varargout] = mrir_array_whitening_operator(noise, varargin);
%MRIR_ARRAY_WHITENING_OPERATOR
%
% W = mrir_array_whitening_operator(noise, METHOD);
%
% W = mrir_array_whitening_operator(covmtx, METHOD);
%
% operator defined such that W'*W = inv(covmtx)
%
%
%  to apply to [Nsamp X Nchan] sensitivity matrix "B":
%
%   Bhat = B * W.';
%
%
%  to apply to [Nchan x 1] single pixel array matrix "s":
%
%   shat = W * s;
%
%
%  see also:  mrir_array_coord_transform

% example:
%
%  noise = randn(100, 4);
%  x = sin(3*(1:100)/2/pi);
%  A = dftmtx(100);
%  W = mrir_array_whitening_operator(noise, 'svd');
%  xhat = x * W.';
%  Ahat = W * A * inv(W);

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/mar/14
% $Id: mrir_array_whitening_operator.m,v 1.3 2011/03/28 04:14:46 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  METHOD = 'svd';        % default
  if ( nargin >= 2 ),
    if ( ~isempty(varargin{1}) ),
      METHOD = varargin{1};
    end;
  end;


  %==--------------------------------------------------------------------==%

  % whitening is based on the covariance matrix

  if ( ndims(noise) == 2 && size(noise,1) == size(noise,2) && ...
       all(imag(diag(noise)) < eps(class(noise))) ),
    % user passed covariance matrix as input
    covmtx = noise;
  elseif ( ndims(noise) == 2 && size(noise,1) > size(noise,2) ),
    % user passed matrix of data where columns are examples (as in "cov.m");
    covmtx = cov(noise).';
  else,
    % user passed ICE-style data matrix, where channels are stored in
    % dimension 3
    covmtx = mrir_array_stats_matrix(double(noise), 'cov');
  end;

  NCha = size(covmtx, 1);

  % SVD-based method is preferable in cases where the covariance matrix
  % is not-invertible or near-singular.


  switch lower(METHOD),

   case 'svd',

    [U, S, V] = svd(covmtx);
    s = diag(S);

    W = diag( sqrt(1./s) ) * U';

   case 'chol',

    invcovmtx = inv(covmtx);
    % due to roundoff error, imag part of diag might not == 0, so fix
    invcovmtx(logical(eye(size(invcovmtx)))) = abs(diag(invcovmtx));

    W = chol(invcovmtx);

   case 'sqrt',

    invcovmtx = inv(covmtx);
    % due to roundoff error, imag part of diag might not == 0, so fix
    invcovmtx(logical(eye(size(invcovmtx)))) = abs(diag(invcovmtx));

    % here, W = W'
    W = sqrtm(invcovmtx);

   otherwise,
    error('unknown whitening method: "%s"', METHOD);
  end;


  if ( nargout >= 2 ),

    dims = size(noise);

    % by convention, ICE data is stored as cols x rows x chan, so reshape to
    % samples x chan
    whitened = reshape(noise, [], dims(end)) * W.';

    % the covariance matrix of the decorrelated data should be the identity!
    decorrmtx = mrir_array_stats_matrix(double(reshape(whitened, [], 1, NCha)), 'cov');

    varargout{1} = decorrmtx;
  end;

  if ( nargout >= 3 )
    % as a sanity check..."invmtx" should be the same as "covmtx"
    % ( W'*W == invcovmtx )
    invmtx = inv(W) * inv(W');

    varargout{2} = invmtx;
  end;

  % RMS error
  whitening_error = sum(sum( (abs(W * covmtx * W') - eye(NCha)).^2 )) / NCha / NCha;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_whitening_operator.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
