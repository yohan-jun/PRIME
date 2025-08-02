function varargout = mrir_filter_raw_apodize_1d(raw, dim_apod, varargin)
%MRIR_FILTER_RAW_APODIZE_1D
%
% raw_apod = mrir_filter_raw_apodize_1d(raw, dim_apod, varargin)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/mar/20
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  % siemens uses two methods: hanning window and the SinusFilter, which
  % is ostensibly the tukey taper
  WINDOW_TYPE = 'tukey';


  %==--------------------------------------------------------------------==%


  dims = size(raw);
  N = dims(dim_apod);

  k = -N/2+1:1:+N/2;


  switch WINDOW_TYPE,

   case {'tukey', 'sinus', 'taper'},

    % choose a tukey tapering window [a.k.a. half-cosine taper] with
    % maximum unattenuated frequency at N/2*(1-k_taper_percent)
    k_taper_percent =  0.15;

    if ( nargin >= 3 ),
      k_taper_percent = varargin{1};
    end;

    if ( k_taper_percent == 0 ),
      % caller requests to skip apodization

      if ( nargout > 0 ),
        varargout{1} = raw;
      end;

      return;
    end;

    k_taper_cutoff = (N/2) * (1 - k_taper_percent);
    w_apod = tukeywin(N, k_taper_percent);

   case {'hanning'},

    w_apod = 0.5 * (1 + cos(pi*k/N))

    % TODO: siemens splits the hanning window according to the desired slope
   
   
   
   case {'none'},
    return;
   otherwise,
    error('!');
  end;


  if ( nargout == 0 ),
    figure; plot(k, w_apod); axis tight;
    return;
  end;


  %==--------------------------------------------------------------------==%

  perm = ones(1,16);
  perm(dim_apod) = N;

  w_apod = reshape(w_apod, perm);

  rep = dims;
  rep(dim_apod) = 1;

  W_apod = repmat(w_apod, rep);

  raw_apod = raw .* W_apod;


  %==--------------------------------------------------------------------==%

  if ( nargout > 0 ),
    varargout{1} = raw_apod;
  end;


  return;


  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
