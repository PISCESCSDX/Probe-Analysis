function [freq kfsp] = kfmean(A, tt, wl, ol)
%==========================================================================
%function [freq kfsp] = kfmean(A, tt, wl, ol)
%--------------------------------------------------------------------------
% Loads the complete couronne matrix (2eSamplesx64, Samples=16,17) and
% calculates the average fft2d. The frequency resolution decreases with the
% amount of averaging steps.
%
% IN: A  matrix [64 x Samples]
%     tt time vector
%     wl length of the window (% of original data length)
%     ol overlap of the windows
%OUT: freq frequency scale
%     kfsp mean kf values (not complex! real!)
%--------------------------------------------------------------------------
% EX:  [freq kfsp] = kfmean(A, tt, 80, 0.5);
%==========================================================================

dt = tt(2)-tt(1);

% overlap of windows (0<ol<1)
if nargin < 4; ol = 0.5; end
if nargin<3; wl = 30; end


% calculate steps and stepsize in dependency on wl and ol
[~, T] = size(A);
W = floor( (T/100)*wl );
% amount of steps
steps = ceil( (T-W)/((1-ol)*W) );
% calculate stepwidth
if steps == 0
  stepwidth = 0;
else
  stepwidth = floor( (T-W)/steps )-1;
end

kfsp = 0; freq = 0;
for i=0:steps
  inda = 1 + i*stepwidth;
  indb = 1 + i*stepwidth + W;
  sig = A(:, inda:indb);
  % calculate the fft2d;  kfsp1 is the complex kf-spectrum
  [freq kfsp1] = fft2d(sig, 1/dt, 1);
  % calculate the mean !! take the absolute value of the windowed
  % kf-spectrum - phase information gets lost
  kfsp = kfsp + abs(kfsp1);
end

% divide by (steps+1) to get the mean; +1 because it starts at 0
kfsp = kfsp/(steps+1);
    
end