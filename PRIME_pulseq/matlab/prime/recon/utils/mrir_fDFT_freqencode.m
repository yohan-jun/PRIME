function ksp = mrir_fDFT_freqencode(img, varargin)
%MRIR_FDFT_FREQENCODE  forward Discrete Fourier Transform along freq encode
%
% ksp = mrir_fDFT_freqencode(img)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/30
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Npoint = size(img, 1);
  if ( nargin >= 2 ),
    Npoint = varargin{2-1};
  end;

  %ksp = ifftshift(fft(fftshift(img, 1), Npoint, 1), 1) / Npoint;

  ksp = fftshift(fft(ifftshift(img, 1), Npoint, 1), 1) / Npoint; % Kawin Sep 2013
  return;


  %************************************************************************
  %%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
