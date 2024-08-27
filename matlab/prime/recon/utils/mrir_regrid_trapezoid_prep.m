function trapezoid = mrir_regrid_trapezoid_prep(MrProt,lNx)
%MRIR_REGRID_TRAPEZOID_PREP
%
% trapezoid = mrir_regrid_trapezoid_prep(MrProt, lNx)

% this code was poached directly from "ice_trapezoid_init.m", by Fa-Hsuan Lin.

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/03
% $Id: mrir_regrid_trapezoid_prep.m,v 1.2 2008/03/30 23:03:09 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  REGRID_NONE        = hex2dec('01');
  REGRID_TRAPEZOIDAL = hex2dec('02');
  REGRID_SINUSOIDAL  = hex2dec('04');

  trapezoid.alRegridMode = MrProt.alRegridMode;

  if ( ~ismember(MrProt.alRegridMode, [REGRID_NONE, REGRID_TRAPEZOIDAL, REGRID_SINUSOIDAL]) ),
    error('unknown regridding mode: [%d]', MrProt.alRegridMode);
  end;

  if ( MrProt.alRegridMode == REGRID_NONE ),
    warning('MRIR:recon:regrid', 'regridding disabled from pulse sequence; skipping...');
    return;
  end;


  %==--------------------------------------------------------------------==%

  warning('off', 'MRIR:recon:regrid');

  % magic numbers from Joe Mandeville's original code
  warning('MRIR:recon:regrid', '*magic numbers* in use');
  regrid_width = 3.5;
  regrid_sigma = 4.91;

  if ( isempty(MrProt.alRegridFlattopTime) ), MrProt.alRegridFlattopTime = 0; end;

  trapezoid.NxRegrid = lNx;

  %// Determine the k-space coordinates as defined by the a sampling scheme
  %// that starts at time lDelaySampling and continues for time MrProt.aflRegridADCDurationDuration
  %// as a trapezoidal waveform plays out using the ramp and flat time as passed.
  if ( MrProt.alRegridMode == REGRID_SINUSOIDAL ),
    renormalize = (MrProt.alRegridFlattopTime) + 4/pi*(MrProt.alRegridRampupTime)*sin(pi/2*((MrProt.aflRegridADCDuration)/2-(MrProt.alRegridFlattopTime)/2)/(MrProt.alRegridRampupTime));
  else
    renormalize = MrProt.aflRegridADCDuration - ( (MrProt.aflRegridADCDuration-MrProt.alRegridFlattopTime)*(MrProt.aflRegridADCDuration-MrProt.alRegridFlattopTime) ) / 4/(MrProt.alRegridRampupTime);
  end;

  renormalize = renormalize / (trapezoid.NxRegrid-1);
  delta_time = (MrProt.aflRegridADCDuration) / (trapezoid.NxRegrid-1);

  kRaw = zeros(trapezoid.NxRegrid,1);
  %// Assign positive k-space values.  Note that the max value is (trapezoid.NxRegrid-1)/2.
  for j = trapezoid.NxRegrid/2:trapezoid.NxRegrid-1,
    time = delta_time * ( j-trapezoid.NxRegrid/2 + 0.5 );
    if ( time < MrProt.alRegridFlattopTime/2 ),
      kRaw(j+1) = time;
    else,
      if ( MrProt.alRegridMode == REGRID_SINUSOIDAL ),
        kRaw(j+1) = (MrProt.alRegridFlattopTime)/2 + 2/pi*(MrProt.alRegridRampupTime)*sin(pi/2*(time-MrProt.alRegridFlattopTime/2.)/MrProt.alRegridRampupTime);
      else,
        kRaw(j+1) = time - ( (time-MrProt.alRegridFlattopTime/2)*(time-MrProt.alRegridFlattopTime/2) ) /2/MrProt.alRegridRampupTime;
      end;
    end;
    kRaw(j+1) = kRaw(j+1)./renormalize;
  end;
  %// Reflect the line to assign negative k-space values.
  for j = 1:trapezoid.NxRegrid/2,
    kRaw(j) = - kRaw(trapezoid.NxRegrid-j+1);
  end;

  %//
  %// Locate the neighbors for each k-space point that fall within a given width.
  %//
  trapezoid.regrid_numberNeighbors = zeros(trapezoid.NxRegrid,1);
  trapezoid.regrid_density         = zeros(trapezoid.NxRegrid,1);
  trapezoid.regrid_rolloff         = zeros(trapezoid.NxRegrid,1);
  trapezoid.regrid_workingSpace    = zeros(trapezoid.NxRegrid,1);
  for  j = 0:trapezoid.NxRegrid-1,
    %       //
    %       // Define the roll-off function.
    %       //
    x = -0.5 + j/(trapezoid.NxRegrid-1);      %// Because delta_k = 1, DELTA_x = 1/delta_k = 1
    trapezoid.regrid_rolloff(j+1) = ice_inverse_kaiser_bessel(x, regrid_width, regrid_sigma);
    %       //
    %       // Locate the neighbors for each k-space point that fall within a given width.
    %       //
    k = kRaw(j+1);
    kNew = j - trapezoid.NxRegrid/2 + 0.5;   %// e.g., -63.5 to 63.5 for 128 steps, in uniform 1-unit increments
    trapezoid.regrid_numberNeighbors(j+1) = 0;
    trapezoid.regrid_density(j+1) = 0;
    for jn = 0:trapezoid.NxRegrid-1,

      kn = kRaw(jn+1);
      if ( abs(k-kn) <= regrid_width/2.0 ),
        trapezoid.regrid_density(j+1) = trapezoid.regrid_density(j+1)+ice_kaiser_bessel( (k-kn), regrid_width, regrid_sigma );
      end;
      if ( abs(kNew-kn) <= regrid_width/2.0 ),
        trapezoid.regrid_numberNeighbors(j+1) = trapezoid.regrid_numberNeighbors(j+1)+1;
      end;
    end;
  end;

  %// Find the maximum number of neighbors, in order to allocate memory.
  maxNeigh = -1;
  for j =0:trapezoid.NxRegrid-1,
    if ( trapezoid.regrid_numberNeighbors(j+1) > maxNeigh ), maxNeigh = trapezoid.regrid_numberNeighbors(j+1); end;
  end;

  %// Allocate 2D matrices.
  trapezoid.regrid_neighbor = zeros(trapezoid.NxRegrid,maxNeigh);
  trapezoid.regrid_convolve = zeros(trapezoid.NxRegrid,maxNeigh);

  %// Fill the 2D matrices.
  for j = 0:trapezoid.NxRegrid-1,
    %// Repeat the above block of code, and store matrix data.
    kNew = j - trapezoid.NxRegrid/2 + 0.5;   %// e.g., -63.5 to 63.5 for 128 steps, in uniform 1-unit increments
    trapezoid.regrid_numberNeighbors(j+1) = 0;

    for jn = 0:trapezoid.NxRegrid-1,

      kn = kRaw(jn+1);
      if ( abs(kNew-kn) <= regrid_width/2 ),
        trapezoid.regrid_neighbor(j+1,trapezoid.regrid_numberNeighbors(j+1)+1) = jn;
        trapezoid.regrid_convolve(j+1,trapezoid.regrid_numberNeighbors(j+1)+1) = ice_kaiser_bessel( (kNew-kn), regrid_width, regrid_sigma );
        trapezoid.regrid_numberNeighbors(j+1) = trapezoid.regrid_numberNeighbors(j+1)+1;
      end;
    end;
  end;

  return;

  

  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_regrid_trapezoid_prep.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
