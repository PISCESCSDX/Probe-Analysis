function [iufit, fitpar] = ...
  chfit_kinetic_mag(fn, pl, pr, par, B, tr, lo, ho, n, f, ...
  epsabs, epsrel, fsaveYN)
%==========================================================================
%function [iufit, fitpar] = ...
%  chfit_kinetic_mag(fn, pl, pr, par, B, tr, lo, ho, n, f, ...
%  epsabs, epsrel, fsaveYN)
%--------------------------------------------------------------------------
% CHFIT_KINETIC_MAG fits Langmuir probe characteristics in magnetized
% plasmas. The fit procedure is done by a software tool 'chfit' written by
% Stefan Zegenhagen. The following assumptions have to hold for the plasma:
% - a cylindrical probe is used
% - probe radius >> electron Larmor radius
% - electron-neutral collisions are dominant (i.e., plasma weakly ionized)
%--------------------------------------------------------------------------
% IN: fn: filename of raw data (rawdata in ascii: voltage, current)
%     pl: probe length (m)
%     pr: probe radius (m)
%    par: [0, 1] probe orientation: perpendicular (0), parallel (1)
%      B: B-field (T)
%     tr: [string] truncate, ex: '-40:15' (V)
%     lo: lower offset plasma potential, -10 - +10 (useful: 0.5V steps)
%     ho: offset plasma potential, -10 - +10 (useful: 0.5V steps)
%      n: smooth over n points
%      f: filter strength, 0-10
% epsabs: look chfit --help=fits for details
% epsrel: look chfit --help=fits for details
%fsaveYN: (0) no files are saved; (1) files are saved; 
%OUT: iufit: structure array containing ne,Te, Vp, Vf, R, RDer
%     fitpar: structure array containing fitparameters
%--------------------------------------------------------------------------
% EX1: [iufit, fitpar] = chfit_kinetic_mag('data01_ui.dat', 0.002, ...
% 0.0005, 0, 0.090, '-80:60', 0.0, 0.0, 5, 1, 1e-6, 1e-6, 1)
% C. Brandt 09.02.2012, San Diego
%==========================================================================

% Use kinetic theory for magnetized plasmas
fittype = 'kinetic_mag';

% Fit procedure
fitvec = ['!chfit'          ...
' -l='     num2str(pl)      ...
' -fit='   fittype          ...
' -r='     num2str(pr)      ...
' -par='   num2str(par)     ...
' -b='     num2str(B)       ...
' -T 1'                     ...
' -t='     tr               ...
' -n='     num2str(n)       ...
' -f='     num2str(f)       ...
' -lo='     num2str(lo)      ...
' -ho='     num2str(ho)      ...
' -epsabs=' num2str(epsabs)  ...
' -epsrel=' num2str(epsrel)  ...
' '       fn];

eval(fitvec);

% Load fit parameters

fitfn = strcat( fn(1:end-4), '.mod');
check = which(fitfn);

if ~isempty(check)
  [~, pa, Rsqr, d_Rsqr] = loadken(fn(1:end-4), 0);
  iufit.ne = pa.n1;     % electron density
  iufit.te = pa.t1;     % electron temperature
  iufit.pp = pa.pp;     % plasma potential
  iufit.fp = pa.fp;     % floating potential
  iufit.R  = Rsqr;      % I-V: quality of fit measure: 1-Chi-squared
  iufit.RDer = d_Rsqr;  % dI/dV: quality of fit measure: 1-Chi-squared
else
  iufit = [];
end

% Store fit parameters
fitpar.pl  = pl;
fitpar.pr  = pr;
fitpar.par = par;
fitpar.B   = B;
fitpar.tr  = tr;
fitpar.nn  = n;
fitpar.ff  = f;
fitpar.ho  = ho;
fitpar.lo  = lo;


if fsaveYN
  savefn = [fn(1:end-4) '_fit.mat'];
  save(savefn, 'iufit', 'fitpar');
end

end