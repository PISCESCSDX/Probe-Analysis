load eva_eva.mat
% variables: EP, EX, I, tt

% PLOT BUSINESS
%==========================================================================

% BUILD THE FOLNAME (for prints)
  filestring='mex_sydw_';
  filetree = fnametree(pwd);
  folname = [filetree{2}(1:8) filetree{1}];
% SET FONTSIZE
  fonts = 14;
% SET TIME INTERVALL
  tint = 2;
  % CREATE TIME INTERVALL END-INDEX
  t = tt*1e3;
  tend = max(find(t<t(1)+tint))+1;
% SET UPPER FREQUENCY
  f_up = 30e3;
% PARAMETER for Axes-Position
  y0=0.25; yd=0.98;  


% EMISSIVE PROBE
%===============
% EP.AC  timerows of AC
% EP.DC  timerows of DC
% EP.f  frequency axis
% EP.amp(:,i)  Power Spectral Amplitude
  % FREQUENCY AXES
    f_ind = find( EP.f <= f_up ); % indices
    frq = EP.f(f_ind);  % frq-axes
    ampS = PwrLog(EP.amp(f_ind, :));
  %----------------------------------
  % PLOT of 3 timerows in one diagram
  %----------------------------------
  fig(12,16,1); p_y=3; i=0;
  ttlim = [min(min(EP.DC(1:tend,:))) max(max(EP.DC(1:tend,:)))];
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  subplottt(tt, EP.DC(:,i), 1, tint, ttlim, 'V_{pl} [V]');
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  subplottt(tt, EP.DC(:,i), 1, tint, ttlim, 'V_{pl} [V]');
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  subplottt(tt, EP.DC(:,i), 1, tint, ttlim, 'V_{pl} [V]');
%   %----------------------------------
%   % PLOT of 3 timerows in one diagram
%   %----------------------------------
%   fig(12,16,1); p_y=3; i=0;
%   ttlim = [min(min(EP.AC)) max(max(EP.AC))];
%   i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
%   subplottt(tt, EP.AC(:,i), 1, tint, ttlim, 'V_{pl} [V]');
%   i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
%   subplottt(tt, EP.AC(:,i), 1, tint, ttlim, 'V_{pl} [V]');
%   i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
%   subplottt(tt, EP.AC(:,i), 1, tint, ttlim, 'V_{pl} [V]');  
% PRINT
print('-deps2c', '-r300', [filestring 'EPtt_' folname '.eps']);
  %----------------------------------
  % PLOT of the 3 frequency spectra
  %----------------------------------
  fig(12,16,2); p_y=3;
% find the fft-limits
    for i=1:size(ampS,2)
%       ampSsmooth(:,i) = smoothspec(ampS(:,i), 20, 0.05, 20, 100, 1);
      ampSsmooth(:,i) = ampS(:,i);
    end;
  fftlim = [min(min(ampSsmooth)) max(max(ampSsmooth))];
% MAKE THE PLOT
  i=0;
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  % subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth(:,i), 1, fftlim);  
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  % subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth(:,i), 1, fftlim);
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  %  subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth(:,i), 1, fftlim);
% PRINT
print('-deps2c', '-r300', [filestring 'EPfft_' folname '.eps']);


% INTERFEROMETER
%===============
% INTF.ndl{i} line integrated density
% radial plot density: INTF.r{i}, INTF.n{i}


% CURRENT MONITOR
%================
% EXCM timerows of signal
% EX.CMamp amplitudes (single value)
% EX.CMf frequencies (single value)