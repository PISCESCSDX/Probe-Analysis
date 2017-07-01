function [freq mvec kfsp] = kfspec(A, tt, mlim, fend, wl, ol)
%==========================================================================
%function [freq mvec kfsp] = kfspec(A, tt, mlim, fend, wl, ol)
%--------------------------------------------------------------------------
% IN: A    matrix [n x 64] n sample, 64 channels (2pi)
%     tt   vector [n x 1] time in s
%     mlim [mode_min mode_max]
%     fend end of frequency axis
%     wl   length of the window (in % of original data)
%     ol   overlap of windows (0<ol<1)
% OUT: freq frequency axis
%      mvec vector with mode numbers
%      kfsp real part of the kfspectrum
%--------------------------------------------------------------------------
% EX:  [freq mvec kfsp] = kfspec(A, tt, [-8 8], 25e3, 40, 0.5);
%==========================================================================

if nargin < 6; ol = 0.5; end
if nargin < 5; wl = 40;  end

% calculate kf-spectrum, kf is already the real part
    [freq kfsp] = kfmean(A, tt, wl, ol);
% calculate the frequency axis
    fend_ind = min(find(freq>fend))-1; %#ok<MXFND>
    freq = freq(1:fend_ind);
% if you want the other part of the freq behind fNy use:
%   freq = freq(end-fend_ind+1:end);
% create mode number vector
    mvec = (mlim(1):1:mlim(2));
% kfsp looks now: mode x freq; mode = 0 ... 32 -32 .. -1
% bring the central mode (m=0) to the center of kfsp
    szkf = floor(size(kfsp, 1)/2);
    kfsp = circshift(kfsp, szkf);
% calculate the amplitudes of kfsp (take abs)
%   +1 and +1 because the center of mvec is 0, but matrix-ind starts with 1
    kfsp = kfsp( mlim(1)+szkf+1:mlim(2)+szkf+1 , 1:fend_ind);

end