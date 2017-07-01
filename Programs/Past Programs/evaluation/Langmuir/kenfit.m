function kenfit(fn, pl, pr, B, tr, n, f, ho)
%==========================================================================
%function kenfit(fn, pl, pr, B, tr, n, f, ho)
%--------------------------------------------------------------------------
% KENFIT evaluates all probe measurement files in one directory.
% For kinetic theory.
%--------------------------------------------------------------------------
% IN: fn: filename string of characteristics, ex: 'meas*.dat'
%     pl: probe length [m]
%     pr: probe radius [m]
%      B: B-field      [T]
%     tr: [string] truncate, ex: '-40:15'
%      n: smooth over n points
%      f: filter strength, 0-10
%     ho: offset plasma potential, -10 - +10 (useful: 0.5V steps)
%--------------------------------------------------------------------------
% EX1: kenfit('meas_*.dat', 0.003, 0.0001, 90, '-80:60', 5, 1, 0.0)
% EX2: kenfit('meas_*.dat', 0.003, 0.0001, 90, '-40:15', 5, 1, 0.0)
% chfit -fit=kinetic_mag -l=0.003 -r=0.0001 -b=0.102 -T 1 -t=-40:15 
%       -n=5 -f=1 -ho=0 meas_000000.dat
% PLOT with:
%   load iu-data
%   kenprofiles(iudat.ne, iudat.pp, iudat.te)
%
% Allgemeines Vorgehen zur IU-Auswertung:
% (0) irgendeine meas_*.dat laden und Spannungsgrenzen ansehen
% (1) alle Kennlinien mit gleichen Parametern mit "kenfit" auswerten
%     mit etwas verkleinerten Spannungsgrenzen aus (1) ermittelt
%     auf die Ausgabe des Wertes "mean-rsqr" achten -> moeglichst > 0.9
% (2) plotten: load iudata; kenprofiles(iudat.ne, iudat.pp, iudat.te)
% (3) checkken anwenden und vergleichen ob die fits passen
%     - ggf. den Spannungsbereich ca. um -40:15 anpassen
%     - ggf. Parameter n, f, ho ändern
% (4) Erneut von (1) beginnen
%==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Characteristic file name list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a  = dir(fn);
al = length(a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go through all characteristic files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:al
% fittype = 'langmuir';
fittype = 'kinetic_mag';
eval(['!chfit -l=' num2str(pl) ...
          ' -fit=' fittype ...
          '   -r=' num2str(pr) ...
          '   -b=' num2str(B) ...
          '   -T 1' ...
          '   -t=' tr ...
          '   -n=' num2str(n) ...
          '   -f=' num2str(f) ...
          '  -ho=' num2str(ho) ...
          ' ' a(i).name]);

  [fit,par,Rsqr,d_Rsqr] = loadken(a(i).name(1:end-4),0);
  iudat.ne(i) = par.n1;
  iudat.te(i) = par.t1;
  iudat.pp(i) = par.pp;
  iudat.fp(i) = par.fp;
  iudat.rsqr(i)  = Rsqr;
  iudat.drsqr(i) = d_Rsqr;
end

fitpar.pl = pl;
fitpar.pr = pr;
fitpar.B  = B;
fitpar.tr = tr;
fitpar.nn = n;
fitpar.ff = f;
fitpar.ho = ho;

save('iu-data.mat', 'iudat', 'fitpar')

disp([ 'mean-rsqr: ' num2str(mean(iudat.rsqr)) ]);

end