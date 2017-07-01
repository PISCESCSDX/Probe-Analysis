function [A_out, ph_out] = calibration_cm(coilno, f0, A_in, ph_in)
%20080224-So-23:25 Brandt
%function [A_out, ph_out] = calibration_cm(coilno, f0, A_in, ph_in)
% Calculate current amplitude [A] and phase [rad] for a signal measured 
% with coil #coilno of the current monitor box (8 Rogowski coils).
%
% APPLY E and Phase to Measurement:
%==================================
% USAGE of CALIBRATION: Get E and phase for a special f0 [Hz]
% E(f0,i)  = polyval(cmcal.E{i}, f0);
% ph(f0,i) = spline(cmcal.f{i}, cmcal.ph{i}, f0);
% => A_real = A_cm / E
% => phi_real = phi_cm - ph(f0)
%
% IN: coilno: no. of the Rogowski coil
%     f0: frequency [Hz]
%     A_in:  input amplitude [V]
%     ph_in: input phase [rad]
%OUT: A_out: output amplitude [A]
%     ph_in: output phase [rad]
%
% EX: [A1, p1] = calibration_cm(1, 2000, 0.5, 1)

if nargin<4; error('More input needed.'); end


% LOAD CALIBRATION FILE
pth_cmcal = '/home/saturn/calibration/current-monitor-box/';
load([pth_cmcal 'current-monitor-cal.mat']);


% APPLY CALIBRATION
  E = polyval(cmcal.E{coilno}, f0);
  A_out = A_in/E;

% APPLY CALIBRATION
  ph = spline(cmcal.f{coilno}, cmcal.ph{coilno}, f0);
  ph_out = ph_in - ph;

end