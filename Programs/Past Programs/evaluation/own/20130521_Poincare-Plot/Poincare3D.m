function [S, dtau, acf] = Poincare3D(s, fs, tau_delay, itint)
%==========================================================================
%function [S, dtau, acf] = Poincare3D(s, fs, tau_delay, itint)
%--------------------------------------------------------------------------
% May-21-2013, C. Brandt, San Diego
% Poincare3D computes the 3D Poincare plot
%--------------------------------------------------------------------------
%INPUT:
% s: data vector
% fs: sampling frequency
% tau_delay: value of time delay (in time units)
%         -1: find time lag automatically
% itint: time index intervall from which signals are generated
%OUTPUT:
% S: structure array: S.s1 = tint, S.s2 = tint +1*tau, S.s3 = tint +2*tau
% dtau: structure: time delay .points, .time
% acf: structure: auto-correlation function
%      .cf auto-correlation function
%      .tau tau vector
%--------------------------------------------------------------------------
% EXAMPLE: 
%==========================================================================

% Calculate auto-correlation function acf
[acf.cf, acf.tau] = ccf(s, s, fs);

%----------------------------------------------------- Calculate time delay
if tau_delay == -1
  % Find number indices between maximum and first null of corr-fct
  dtau.points = ccf_timedelay(acf.cf);
  dtau.time = dtau.points ./ fs;
else
  % Set the time delay manually
  dtau.points = round( tau_delay .* fs );
  dtau.time = tau_delay;
end

%------------------------------------- Define the three time lagged signals
S.s1 = s( itint + 0*dtau.points );
S.s2 = s( itint + 1*dtau.points );
S.s3 = s( itint + 2*dtau.points );

end