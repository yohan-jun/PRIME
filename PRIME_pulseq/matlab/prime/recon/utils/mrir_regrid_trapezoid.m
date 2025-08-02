function varargout = mrir_regrid_trapezoid(img, prot, varargin)
%MRIR_REGRID_TRAPEZOID
%
% img_deapod = mrir_regrid_trapezoid(img, prot)
%
% function contains an fDFT, so the regridding step can be applied anywhere
% in the reconstruction process after the first iDFT.
%
% parameters required: aflRegridADCDuration, alRegridRampupTime, alRegridRampdownTime

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/30
% $Id: mrir_regrid_trapezoid.m,v 1.1 2007/01/30 22:21:05 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  if ( nargin >= 3 ),
    prot_trapezoid = varargin{1};
  else,
    prot_trapezoid = [];
  end;

  DO__DEAPOD = 1;
  if ( nargin >= 4 ),
    DO__DEAPOD  = varargin{2};
  end;
  
  
  
  %==--------------------------------------------------------------------==%

  % before becoming vulnerable to obscure error messages, make sure that
  % prot contains all the relevant parameter values
  if ( isempty(prot.aflRegridADCDuration) || ...
       isempty(prot.alRegridRampupTime)   || ...
       isempty(prot.alRegridRampdownTime) ),
    error('protocol is missing parameter values needed for regridding');
  end;


  % ASSUME: waveform is symmetric, i.e., rise time == fall time

  if ( isempty(prot_trapezoid) ),
    % calculate regriding parameters from sequence protocol parameter values
    prot_trapezoid = mrir_regrid_trapezoid_prep(prot, size(img, 1));
  end;
  
  % project image into k-space for proper resampling
  raw = mrir_fDFT_freqencode(img);

  raw_regrid = mrir_regrid_trapezoid_apply(raw, prot_trapezoid);

  % project data back into x-space for deapodization
  img_regrid = mrir_iDFT_freqencode(raw_regrid);

  if ( DO__DEAPOD ),
    img_deapod = mrir_regrid_trapezoid_rolloff(img_regrid, prot_trapezoid);
    varargout{1} = img_deapod;
  else,
%    disp('skipping deapodization');
    varargout{1} = img_regrid;
  end;
  
  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/EPI_3D/MATLAB/mrir_regrid_trapezoid.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
