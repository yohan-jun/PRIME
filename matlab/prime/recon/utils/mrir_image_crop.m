function img_cropped = mrir_image_crop(img_oversampled, varargin)
%MRIR_IMAGE_CROP  crop image volume reconstructed from oversampled k-space
%
% cropped = mrir_image_crop(img_oversampled)
% cropped = mrir_image_crop(img_oversampled, prot)
% cropped = mrir_image_crop(img_oversampled, os_factor_freqencode, os_factor_phasencode, os_factor_partencode)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/16
% $Id: mrir_image_crop.m,v 1.3 2009/02/15 00:48:06 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global VERBOSE; if ( isempty(VERBOSE) ), VERBOSE = 0; end;
  

  %==--------------------------------------------------------------------==%

  % NOTE: since reconstructed volumes are "squeezed" and tend to have
  % partitions along the THIRD dimension (instead of the ninth), assume that
  % this is the case. this assumption will be violated if (a) there are more
  % than one repetition, which would push the partitions to the fourth
  % dimension after squeezing, or (b) if no squeeze has been applied, which
  % is the case before coil combination occurs.  
  
  os_factor_freqencode = [];
  os_factor_phasencode = [];
  os_factor_partencode = [];

  if ( nargin >= 2 && isstruct(varargin{1}) ),
    prot = varargin{1};

    os_factor_freqencode = prot.flReadoutOSFactor;
    os_factor_phasencode = prot.dPhaseOversamplingForDialog + 1;
    os_factor_partencode = prot.dSliceOversamplingForDialog + 1;
  else,

    if ( nargin >= 2 )
      os_factor_freqencode = varargin{1};
    end;
    if ( nargin >= 3 )
      os_factor_phasencode = varargin{2};
    end;
    if ( nargin >= 4 )
      os_factor_partencode = varargin{3};
    end;

  end;

  % by default, only oversampling is only along first dimension (i.e.,
  % freqencode).

  if ( isempty(os_factor_freqencode) ),
    os_factor_freqencode = 2;  % default
  end;
  if ( isempty(os_factor_phasencode) ),
    os_factor_phasencode = 1;  % default
  end;
  if ( isempty(os_factor_partencode) ),
    os_factor_partencode = 1;  % default
  end;


  dims = size(img_oversampled);

  if ( os_factor_freqencode > 1 && mod( round(dims(1) / os_factor_freqencode), 1) ~= 0 ),
    error('samples are not integer multiple of oversampling factor along frequency-encoded direction?!');
  end;

  if ( os_factor_phasencode > 1 && mod( round(dims(2) / os_factor_phasencode), 1) ~= 0 ),
    error('samples are not integer multiple of oversampling factor along phase-encoded direction?!');
  end;

  if ( os_factor_partencode > 1 && mod( round(dims(3) / os_factor_partencode), 1) ~= 0 ),
    error('samples are not integer multiple of oversampling factor along partition-encoded direction?!');
  end;
  

  %==--------------------------------------------------------------------==%

  % center of k-space is at N/2+1, i.e., there are more negative
  % k-space indices than positive.

  img_cropped = img_oversampled;

  if ( os_factor_freqencode > 1 ),

    newdim = round(dims(1)/os_factor_freqencode);

    if ( VERBOSE ),
      dstr = sprintf('cropping FREQENCODE dimension:  %4d -> %4d', dims(1), newdim);
      disp(sprintf('==> [%s]: %s', mfilename, dstr));
    end;

    index_extract = round([1:newdim] + dims(1) * (os_factor_freqencode-1)/(2*os_factor_freqencode));

    %                                     1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
    img_cropped = img_cropped(index_extract,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);

  end;

  if ( os_factor_phasencode > 1 ),

    newdim = round(dims(2)/os_factor_phasencode);
    
    if ( VERBOSE ),
      dstr = sprintf('cropping PHASENCODE dimension:  %4d -> %4d', dims(2), newdim);
      disp(sprintf('==> [%s]: %s', mfilename, dstr));
    end;
    
    index_extract = round([1:newdim] + dims(2) * (os_factor_phasencode-1)/(2*os_factor_phasencode));

    %                         1             2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
    img_cropped = img_cropped(:,index_extract,:,:,:,:,:,:,:,:,:,:,:,:,:,:);

  end;

  if ( os_factor_partencode > 1 ),

    newdim = round(dims(3)/os_factor_partencode);

    if ( VERBOSE ),
      dstr = sprintf('cropping PARTENCODE dimension:  %4d -> %4d', dims(3), newdim);
      disp(sprintf('==> [%s]: %s', mfilename, dstr));
    end;
      
    index_extract = round([1:newdim] + dims(3) * (os_factor_partencode-1)/(2*os_factor_partencode));

    %                         1 2             3 4 5 6 7 8 9 0 1 2 3 4 5 6
    img_cropped = img_cropped(:,:,index_extract,:,:,:,:,:,:,:,:,:,:,:,:,:);

  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_crop.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
