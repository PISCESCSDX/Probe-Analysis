load eva_eva.mat
% variables: BP, EP, EX, I, TP, tt

% PLOT BUSINESS
%==========================================================================

% BUILD THE FOLNAME (for prints)
  filestring='exc_sytu_';
  filetree = fnametree(pwd);
  folname = [filetree{2}(1:8) filetree{1}];
% SET FONTSIZE
  fonts = 14;
% SET TIME INTERVALL
  tint = 2;
% SET UPPER FREQUENCY
  f_up = 50e3;
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
  ttlim = [min(min(EP.AC)) max(max(EP.AC))];
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  subplottt(tt, EP.AC(:,i), 1, tint, ttlim, 'V_{pl} [V]');
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  subplottt(tt, EP.AC(:,i), 1, tint, ttlim, 'V_{pl} [V]');
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  subplottt(tt, EP.AC(:,i), 1, tint, ttlim, 'V_{pl} [V]');
% PRINT
print('-deps2c', '-r300', ['EPtt_' folname '.eps']);
  %----------------------------------
  % PLOT of the 3 frequency spectra
  %----------------------------------
  fig(12,16,2); p_y=3;
% find the fft-limits
    for i=1:size(ampS,2)
      ampSsmooth(:,i) = smoothspec(ampS(:,i), 20, 0.05, 20, 100, 1);
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


% BDOT PROBE
%===========
  f_ind = find( BP.f <= f_up ); % indices
  frq = BP.f(f_ind);  % frq-axes
  ampS1 = PwrLog(BP.amp1(f_ind, :));
  ampS2 = PwrLog(BP.amp2(f_ind, :));
  ampS3 = PwrLog(BP.amp3(f_ind, :));
% BP.amp1, BP.amp3, BP.amp3
  %----------------------------------
  % PLOT of the 3 frequency spectra
  %----------------------------------
% find the fft-limits
    for i=1:size(ampS1,2)
      sig = smoothspec(ampS1(:,i), 20, 0.03, 20, 100, 0);
      ampSsmooth1(:,i) = sig;
      sig = smoothspec(ampS2(:,i), 20, 0.03, 20, 100, 0);
      ampSsmooth2(:,i) = sig;
      sig = smoothspec(ampS3(:,i), 20, 0.03, 20, 100, 0);
      ampSsmooth3(:,i) = sig;      
    end;
% PLOT BDOT-PROBE 1
  fig(12,16,3); p_y=3;
  fftlim = [min(min(ampSsmooth1)) max(max(ampSsmooth1))];
% MAKE THE PLOT
  i=0;
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  % subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth1(:,i), 1, fftlim);  
  %
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  % subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth1(:,i), 1, fftlim);
  %
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  %  subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth1(:,i), 1, fftlim);
% PRINT
print('-deps2c', '-r300', [filestring 'BP1fft_' folname '.eps']);  
%
% PLOT BDOT-PROBE 2
  fig(12,16,4); p_y=3;
  fftlim = [min(min(ampSsmooth2)) max(max(ampSsmooth2))];
% MAKE THE PLOT
  i=0;
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  % subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth2(:,i), 1, fftlim);  
  %
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  % subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth2(:,i), 1, fftlim);
  %
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  %  subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth2(:,i), 1, fftlim);
% PRINT
print('-deps2c', '-r300', [filestring 'BP2fft_' folname '.eps']);
%
% PLOT BDOT-PROBE 2
  fig(12,16,5); p_y=3;
  fftlim = [min(min(ampSsmooth3)) max(max(ampSsmooth3))];
% MAKE THE PLOT
  i=0;
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  % subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth3(:,i), 1, fftlim);  
  %
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  % subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth3(:,i), 1, fftlim);
  %
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  %  subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth3(:,i), 1, fftlim);
% PRINT
print('-deps2c', '-r300', [filestring 'BP2fft_' folname '.eps']);


% TRANSPORT PROBE
%================
  % FREQUENCY AXES
    f_ind = find( TP.f <= f_up ); % indices
    frq = TP.f(f_ind);  % frq-axes
    ampS = PwrLog(TP.amp(f_ind, :));
  % FREQUENCY AXES
  %----------------------------------
  % PLOT of 3 timerows in one diagram
  %----------------------------------
  fig(12,16,6); p_y=3; i=0;
  ttlim = [min(min(TP.GammaTilde)) max(max(TP.GammaTilde))];
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  subplottt(tt, TP.GammaTilde(:,i), 1, tint, ttlim, '\Gamma [a.u.]');
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  subplottt(tt, TP.GammaTilde(:,i), 1, tint, ttlim, '\Gamma [a.u.]');
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  subplottt(tt, TP.GammaTilde(:,i), 1, tint, ttlim, '\Gamma [a.u.]');
% PRINT
print('-deps2c', '-r300', [filestring 'TPtt_' folname '.eps']);
  %----------------------------------
  % PLOT of the 3 frequency spectra
  %----------------------------------
  fig(12,16,7); p_y=3;
% find the fft-limits
    for i=1:size(ampS,2)
      ampSsmooth(:,i) = smoothspec(ampS(:,i), 20, 0.05, 20, 100, 1);
    end;
  fftlim = [min(min(ampSsmooth)) max(max(ampSsmooth))];
% MAKE THE PLOT
  i=0;
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  % subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth(:,i), 1, fftlim);  
  %
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  % subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth(:,i), 1, fftlim);
  %
  i=i+1; axes('Position', [0.20 ((p_y-i)/p_y+(i+2)/50) 0.7 0.22]);
  %  subplotfft(frq, ampS(:,i), 1, fftlim);
  subplotfft(frq, ampSsmooth(:,i), 1, fftlim);
% PRINT
print('-deps2c', '-r300', [filestring 'TPfft_' folname '.eps']);


% INTERFEROMETER
%===============
% INTF.ndl{i} line integrated density
% radial plot density: INTF.r{i}, INTF.n{i}


% EXCITER-CURRENTS
%=================
% EX.I(shot-no, time, coil/electrode): currents /A
% EX.EXCfrq: frq of maxamp of all 8
% EX.EXCamp: max ampat frq of all 8
% EX.EXCpha: phase of frq of all 8
  %----------------------------------
  % PLOT of 3 timerows in one diagram
  %----------------------------------
  p_y=1; fig(7,16/3*p_y,8); i=0; tint=2;
  ttlim = [min(min(min(EX.I(:,:,1)))) max(max(max(EX.I(:,:,1))))];
  i=i+1; axes('Position', [0.25 ((y0+p_y-i)*yd)/p_y 0.7 0.7/p_y]);
  subplottt(tt, EX.I(1,:,1), 1, tint, ttlim, 'I_{ex} [A]');
print('-deps2c', '-r300', ['EXtt_' folname '1.eps']);  
%   p_y=1; fig(5,16/3*p_y,9); i=0; tint=2;
%   ttlim = [min(min(min(EX.I(:,:,1)))) max(max(max(EX.I(:,:,1))))];
%   i=i+1; axes('Position', [0.20 ((y0+p_y-i)*yd)/p_y 0.7 0.7/p_y]);
%   subplottt(tt, EX.I(2,:,1), 1, tint, ttlim, 'I_{ex} [A]');
% print('-deps2c', '-r300', ['EXtt_' folname '2.eps']);
  p_y=1; fig(7,16/3*p_y,10); i=0; tint=2;
  ttlim = [min(min(min(EX.I(:,:,1)))) max(max(max(EX.I(:,:,1))))];
  i=i+1; axes('Position', [0.25 ((y0+p_y-i)*yd)/p_y 0.7 0.7/p_y]);
  subplottt(tt, EX.I(3,:,1), 1, tint, ttlim, 'I_{ex} [A]');
print('-deps2c', '-r300', ['EXtt_' folname '3.eps']);


% CURRENT MONITOR
%================
% EXCM timerows of signal
% EX.CMamp amplitudes (single value)
% EX.CMf frequencies (single value)