function [img_combine_rss, varargout] = mrir_array_combine_rss(img_multichan, varargin)
%MRIR_ARRAY_COMBINE_RSS
%
% img_combine_rss = mrir_array_combine_rss(img_multichan);

% TODO: return measure of phase compatible with root-sum-of-squares
% combination, e.g., homodyne detection (Noll et al., 1991)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/27
% $Id: mrir_array_combine_rss.m,v 1.1 2008/01/27 00:27:37 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  alpha = 1;
  
  if ( nargin >= 2 ),
    covmtx = varargin{1};
  
    % arbitrary constant (normalizes image to approx same level as RSS combo)
    alpha = sqrt(mean(diag(covmtx)));

    W = mrir_array_whitening_operator(covmtx, 'svd');
    img = mrir_array_whitening_apply(img_multichan, W);
    
  else,
    img = img_multichan;
  end;
  
  img_combine_rss = alpha * squeeze(sqrt(sum(abs(img).^2,3)));
%  img_combine_phz = angle(sum(img_multichan, 3));
  
%  img_combine_cplx = img_combine_rss .* exp(i * img_combine_phz);
  
%  if ( nargout > 1 ),
    
%    varargout{1} = img_combine_phz;
%    varargout{2} = img_combine_cplx;
    
%  end;
  
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_combine_rss.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: