function roft = mrir_iDFT_freqencode(raw, varargin)
%MRIR_IDFT_FREQENCODE  inverse Discrete Fourier Transform along freq encode
%
% roft = mrir_iDFT_freqencode(raw)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/30
% $Id: mrir_iDFT_freqencode.m,v 1.1 2007/09/17 20:18:36 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Npoint = size(raw, 1);
  if ( nargin >= 2 && ~isempty(varargin{2-1}) ),
    Npoint = varargin{2-1};
  end;

  roft = mrir_iDFT(raw, 1, Npoint);

  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_iDFT_freqencode.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
