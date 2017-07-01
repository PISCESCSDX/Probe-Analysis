function loadprofiles
%LOADPROFILES loads the evaluation files of a directory containing
% the data of a radial characteristic profile.

filetree = fnametree(pwd);
folname = [filetree{3}(1:8)];

ad  = dir('meas*dat');
al = length(ad);

for i=1:al
  fname = ad(i).name(1:end-4);
  fn=my_filename(i-1,6,'meas_','');
  [fit,par] = loadken( fname );
  n(i) = par.n1;
  t(i) = par.t1;
  pp(i) = par.pp;
  fp(i) = par.fp;
  
  a=load([fname '.dat']);
  is(i) = mean(a(1:100,2));
end

figeps(16,25,1)

subplot(321); plot(n, 'bo-');
title('density profile');

subplot(322); plot(-is, 'bo-');
title('I_{sat} profile');

subplot(323); plot(t, 'ro-'); ylabel('T [eV]');
title('T-profile');

subplot(324); plot(pp, 'ko'); ylabel('V_{p} [V]');
title('V_p profile');

subplot(325); plot(fp, 'ko'); ylabel('V_{f} [V]');
title('V_f profile');

print('-depsc2',  ['iueva' folname filetree{2} filetree{1} '.eps']);

save profile_data.mat n t pp

end