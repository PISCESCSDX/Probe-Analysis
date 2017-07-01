function [tvec mvec spec] = mspectrogram(tt, sig, pts)
%==========================================================================
%function [tvec fvec fmat] = fspectrogram(tt, sig, pts)
%--------------------------------------------------------------------------
% MSPECTROGRAM calculates the mode spectogram of a 2d-signal "sig".
% (20100906-Mi-19:23 C. Brandt)
%--------------------------------------------------------------------------
% IN: tt: time vector
%     sig: 2d space-time series
%     pts: distance between ffts
%OUT: tvec: time step vector
%     mvec: frequency axis
%     spec: time-mspec-matrix [dBV]
% EX: [tvec, mvec, spec] = fspectrogram(tt, sig, fwinpts)
%==========================================================================

% Time vector
  tvec = tt(1:pts:end);
  steps = floor( (length(tt)+1)/pts )+1;

spec = zeros(length(sig(1,:))/2, steps);
disp('calculate mspectra ...');
for i=1:steps
  sigi = sig((i-1)*pts+1, :);
  ls = length(sigi); lwave = ( ( 1:ls )' )/ls;
  [mvec amp] = fftspec(lwave, sigi', 1, 0.1);
  spec(:,i) = amp;
end

end