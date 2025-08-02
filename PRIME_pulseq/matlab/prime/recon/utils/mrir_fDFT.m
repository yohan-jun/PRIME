function ksp = mrir_fDFT(img, dim, varargin)
%MRIR_FDFT  forward Discrete Fourier Transform
%
% ksp = mrir_fDFT(img)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/30
% $Id: mrir_fDFT.m,v 1.1 2009/01/30 06:09:47 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  if ( nargin >= 3 ),
    Npoint = varargin{1};
  else,
    Npoint = size(img, dim);
  end;
  
  
  %==--------------------------------------------------------------------==%

  if ( Npoint > size(img, dim) ),

    % because k-space ordering requires FFT shifting, we must zero pad
    % before the shift, thus simply computing an N-point FFT will not
    % suffice.
    pad_dims = zeros(1, 16);

    pad_dims(dim) = Npoint - size(img, dim);

    % by convention, always pad AFTER the data
    img = padarray(img, pad_dims, 'post');

  end;

  %ksp = ifftshift(fft(fftshift(img, dim), Npoint, dim), dim) / Npoint;

  ksp = fftshift(fft(ifftshift(img, dim), Npoint, dim), dim) / Npoint; % Kawin Sep 2013
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_fDFT.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

  
  
  
  % note on scaling:
  %
  %  in here and in "mrir_iDFT.m", i have reversed the MATLAB convention for
  %  scaling the fft. MATLAB does not like the 'unitary' normalization that
  %  places the \sqrt{N} before both the fDFT and iDFT. see the help for
  %  "fft.m":
  % 
  %%       For length N input vector x, the DFT is a length N vector X,
  %%       with elements
  %%                        N
  %%          X(k) =       sum  x(n)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.
  %%                       n=1
  %%       The inverse DFT (computed by IFFT) is given by
  %%                        N
  %%          x(n) = (1/N) sum  X(k)*exp( j*2*pi*(k-1)*(n-1)/N), 1 <= n <= N
  %%                       k=1
  %
  %  it would have been better if i originally stuck with the unitary
  %  version, but other code (e.g., the SNR calculation code) assumes the
  %  scaling convention established in this function, so i am afraid to
  %  change it now. i have given the 1/N to the forward, and 1 to the
  %  inverse, so that still the composition of the two yields the identity,
  %  but this makes comparison with operators using the unitary
  %  scaling---like the DFT matrix.
