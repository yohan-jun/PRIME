function varargout = mrir_ice_dimensions(data, varargin);
%MRIR_ICE_DIMENSIONS  list extent of array along the 16 canonical ICE dimensions
%
% dimensions = mrir_ice_dimensions(data, varargin)
%
%
% example:
%
%   mrir_ice_dimensions(meas.data);
%   mrir_ice_dimensions([]);          % dump dimension ordering
%
%   mrir_ice_dimensions(meas.data, 'Cha');  % size along "Cha" dimension
    
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/18
% $Id: mrir_ice_dimensions.m,v 1.2 2011/03/28 04:14:46 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  dimensions = struct;

  dim_list = {'NColMeas', 'NLinMeas', 'NChaMeas', 'NSetMeas', 'NEcoMeas', ...
              'NPhsMeas', 'NRepMeas', 'NSegMeas', 'NParMeas', 'NSlcMeas', ...
              'NIdaMeas', 'NIdbMeas', 'NIdcMeas', 'NIddMeas', 'NIdeMeas', ...
              'NAveMeas'};

  % initialize
  for ind = 1:length(dim_list),
    dimensions = setfield(dimensions, dim_list{ind}, 0);
  end;


  % data_dims <= length(dim_list)
  data_dims = length(size(data));

  if ( nargin >= 2 ),

    dim_str = varargin{1};

    index = find(ismember(lower(dim_list), lower(dim_str)));
    
    if ( isempty(index) ),
      dim_list_abbrev = regexprep(dim_list, {'^N', 'Meas$'}, '');
      index = find(ismember(lower(dim_list_abbrev), lower(dim_str)));
    end;
    
    if ( isempty(index) ),
      error('dimension string "%s" not recognized', dim_str);
    end;

    varargout{1} = size(data, index);
    return;

  end;


  for ind = 1:data_dims,
    dimensions = setfield(dimensions, dim_list{ind}, size(data, ind));
  end;


  if ( nargout > 0 ),
    varargout{1} = dimensions;
    return;
  end;

  fprintf(1, '\n');
  for ind = 1:length(dim_list),
    disp(sprintf('  [#%02d]  %s: %d', ind, dim_list{ind}, getfield(dimensions, dim_list{ind})));
  end;
  fprintf(1, '\n');

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_ice_dimensions.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
