function [noise_bandwidth, noise_bandwidth_chan] = mrir_noise_bandwidth(noise, varargin)
%MRIR_NOISE_BANDWIDTH
%
% bandwidth = mrir_noise_bandwidth(noise)
%
% 
% references:
%
%  Kellman P, McVeigh ER.
%  Image reconstruction in SNR units: a general method for SNR measurement.
%  Magn Reson Med. 2005 Dec;54(6):1439-47.
%  PMID: 1626157  
  
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jun/01
% $Id: mrir_noise_bandwidth.m,v 1.2 2008/02/28 05:38:20 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  DISPLAY = 0;
  if ( nargin >= 2 ),
    DISPLAY = varargin{1};
  end;


  %==--------------------------------------------------------------------==%

  dims = size(noise);
  dims(end+1:16) = 1;
   
  % reshape to samples X chan
  order = 1:length(dims);
  order(3) = length(dims);
  order(end) = 3;

  % move channel index to end
  noise = permute(noise, order);
  noise = reshape( noise, [], dims(3) );

  % compute noise bandwidth (Kellman & McVeigh, 2005)
  N = fft(reshape(noise, dims(1), [], dims(3)), [], 1);
  N_power_spectrum = abs(N).^2;
  N_power_spectrum_chan = squeeze(mean(N_power_spectrum, 2));
  N_power_spectrum_chan_normalized = N_power_spectrum_chan./repmat(N_power_spectrum_chan(1, :), [dims(1), 1]);
  noise_bandwidth_chan = mean(N_power_spectrum_chan_normalized, 1);

  
  %  N = fft(reshape(noise, dims(1), []), [], 1);
  %  N_power_spectrum_avg = mean( abs(N).^2, 2);
  N_power_spectrum_avg = mean(N_power_spectrum_chan, 2);
  N_power_spectrum_avg_normalized = N_power_spectrum_avg/N_power_spectrum_avg(1);
  noise_bandwidth = mean(N_power_spectrum_avg_normalized);


  if ( DISPLAY ),
    figure('name', mfilename);
    rh = rectangle('position', [0.25*dims(1), 0.001, 0.5*dims(1), 1.098]);
    set(rh, 'FaceColor', 0.9*[1 1 1], 'LineStyle', 'none');
    hold on;
    plot(fftshift(N_power_spectrum_avg_normalized));
    set(gca, 'XTick', [1, dims(1)*[0.25 0.50 0.75 1.0]], 'XTickLabel', ...
             [-0.5, -0.25 0, +0.25 +0.5]);
    set(gca, 'XLim', [1, dims(1)], 'YLim', [0, 1.1]);
    set(gca, 'YGrid', 'on');
    xlabel('frequency (normalized)');
    ylabel('power (DC normalized)');
    title(sprintf('average of normalized noise power spectrum, BW=%2.3f', noise_bandwidth));
    box on;
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_noise_bandwidth.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:



%%%  N = 10000;
%%%
%%%  noise_true = complex(randn(N,1),randn(N,1)) * 3.0;
%%%  noise_freq = fft(noise_true);
%%%
%%%  T = (tukeywin(N, 0.5) + 0.4) / 1.4;
%%%
%%%  noise_filt = abs((noise_freq) .* T) .* exp(i*angle(noise_freq));
%%%
%%%  BW = sum(T.^2)/N
%%%
%%%  noise_meas = (ifft(noise_filt));
%%%
%%%
%%%  var(real(noise_meas)) / var(real(noise_true))

