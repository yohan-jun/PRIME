function ksp = mrir_fDFT_phasencode(img, coordinate_str, varargin)
%MRIR_FDFT_PHASENCODE  forward Discrete Fourier Transform along phase encode
%
% ksp = mrir_fDFT_phasencode(img, coordinate_str)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/24
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % assume primary phase encode direction by default
  if ( ~exist('coordinate_str', 'var') ),
    coordinate_str = 'lin';
  end;

  switch lower(coordinate_str(1:3)),
   case 'lin',
    dim = 2;
   case 'par',
    dim = 9;
   otherwise,
    dim = 0;
  end;

  if ( size(img, dim) == 1 ),
    disp(sprintf('input "%s" contains no data long dimension "%s"', ...
                 inputname(1), coordinate_str));
  end;

  Npoint = size(img, dim);
  if ( nargin >= 3 ),
    Npoint = varargin{3-2};
  end; 
  

  %ksp = ifftshift(fft(fftshift(img, dim), Npoint, dim), dim) / Npoint;
  ksp = fftshift(fft(ifftshift(img, dim), Npoint, dim), dim) / Npoint; % Kawin Sep 2013
  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
  