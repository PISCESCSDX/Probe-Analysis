function [dtau,pa] =Velocimetry2D_Extract_TimeDelay(s1,s2,...
  fs,winrat,olap,fend,k,calcphstd)
%==========================================================================
%function [dtau,pa]
%=Velocimetry2D_Extract_TimeDelay(s1,s2,fs,winrat,olap,fend,k,calcphstd)
%--------------------------------------------------------------------------
% May-15-2013, C. Brandt, San Diego
% Velocimetry2D_Extract_TimeDelay calculates the time lag between two
% signals 's1' and 's2'.
%--------------------------------------------------------------------------
%IN:
% s1: signal 1
% s2: signal 2
% calcphstd: if 0 much faster!
%OUT:
% dtau: time lag vector (frequency dependent)
% pa: structure: .cpsd: cross power spectrum
%                .freq: frequency vector
%                .phstd: standard deviation for phase (and time lags)
%--------------------------------------------------------------------------
%EXAMPLE:
%==========================================================================

% Calculate cross power spectral density
[A,ph,frq,phstd] = cpsdmean(s1,s2,fs,winrat,olap,fend,k,calcphstd);

% Calculate time lag dtau
%   omega = d theta / d t
%      dt = d theta / (2*pi*f)
dtau = ph*pi ./ (2*pi*frq);

pa.cpsd = A;
pa.ph = ph;
pa.freq = frq;
pa.phstd = phstd;

end