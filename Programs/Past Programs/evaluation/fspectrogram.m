function [tvec fvec spec] = fspectrogram(tt, sig, fwinpts)
%==========================================================================
%function [tvec fvec fmat] = fspectrogram(tt, sig, fwinpts)
%--------------------------------------------------------------------------
% FSPECTROGRAM calculates the fft frequency spectogram of a signal "sig".
% (20091125-Mi-12:23 C. Brandt)
%--------------------------------------------------------------------------
% IN: tt: time vector
%     sig: time series
%     fwinpts: samples per f-window (dividable by 2!)
%OUT: tvec: time step vector
%     fvec: frequency axis
%     spec: time-fspec-matrix [dBV]
% EX: [tvec, fvec, spec] = fspectrogram(tt, sig, fwinpts)
%==========================================================================

% Window overlap
  olap = 0.5;
% Time interval
  dt = tt(2)-tt(1);
% Calculate steps
  dwin = round(olap*fwinpts);
  steps = floor((length(tt)-fwinpts)/dwin);
% Time vector
  tvec = (( 1:steps )-1).*(dwin-1)*dt;


disp('calculate fspectra ...');
spec = zeros(floor(fwinpts/2),steps);
for i=1:steps
  ind = (1:fwinpts)+(i-1)*dwin;
  [fvec amp] = fftspec(tt(ind), sig(ind), 1, 0.1);
  spec(:,i) = amp;
end

end