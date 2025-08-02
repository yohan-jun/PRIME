function img_ordered = mrir_image_slice_deinterleave(img_interleaved)
%MRIR_IMAGE_SLICE_DEINTERLEAVE  sort interleaved slices into ascending order
%
% img_ordered = mrir_image_slice_deinterleave(img_interleaved)
%
% assumes that SLC dimension is 10, following ICE convention.
% 
% (wrapper around "tdr_sliceorder")
  
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/feb/05
% $Id: mrir_image_slice_deinterleave.m,v 1.1 2007/10/03 05:22:44 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Nslc = mrir_ice_dimensions(img_interleaved, 'slc');
  %                             1 2 3 4 5 6 7 8 9        1 2 3 4 5 6
  img_ordered = img_interleaved(:,:,:,:,:,:,:,:,:, ...
                                tdr_sliceorder(Nslc, 2), :,:,:,:,:,:);

  
  return;
  
  

  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_image_slice_deinterleave.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
