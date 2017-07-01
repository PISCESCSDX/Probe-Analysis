function allkenfit(pl, pr, curr)
%function allkenfit(pl, pr, curr)
% Fits all meas*.dat files with chfit.
% (function AUTOKENFIT is more advanced!)
% IN: pl: probe length [m]
%     pr: probe radius [m]
%     curr: coil current [A]
%OUT: evaluation files for all fits
% EX: allkenfit(0.004, 0.0001, 160);
% look also for "autokenfit", "loadprofiles"


a = dir('meas*.dat');
al = length(a);

B = bfield(curr);
  
for i=1:al
eval(['!chfit -l=' num2str(pl) ...
      ' -fit=kinetic_mag' ...
      ' -r=' num2str(pr) ...
      ' -b=' num2str(B) ...
      ' -t=-198:198' ...
      ' -n=1' ...
      ' -f=2' ...
      ' -ho=2' ...
      ' ' a(i).name]);
end;

% get the radial profile
for i=1:al
  [fit,par] = loadken(a(i).name(1:end-4),1);
  dens(i) = par.n1;
end;

end