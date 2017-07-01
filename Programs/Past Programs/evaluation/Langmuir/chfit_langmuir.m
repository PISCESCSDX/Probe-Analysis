function chfit_langmuir(fn, pl, pr, B, tr, n, f, ho)
%==========================================================================
%function chfit_langmuir(fn, pl, pr, B, tr, n, f, ho)
%--------------------------------------------------------------------------
% CHFIT_KINETIC_MAG fits Langmuir probe characteristics in magnetized
% plasmas. The fit procedure is done by a software tool 'chfit' written by
% Stefan Zegenhagen. The following assumptions have to hold for the plasma:
% - a cylindrical probe is used
% - plasma collisionless (probe radius << mean free path of electrons)
% - plasma not magnetized
%--------------------------------------------------------------------------
% IN: fn: filename of raw data (rawdata: two columns: voltage, current)
%    fit: 
%     pl: probe length [m]
%     pr: probe radius [m]
%      B: B-field      [T]
%     tr: [string] truncate, ex: '-40:15'
%      n: smooth over n points
%      f: filter strength, 0-10
%     ho: offset plasma potential, -10 - +10 (useful: 0.5V steps)
%--------------------------------------------------------------------------
% EX1: kenfit('meas_*.dat', 0.003, 0.0001, 90, '-80:60', 5, 1, 0.0)
%==========================================================================

error('This function is not ready yet.')

end