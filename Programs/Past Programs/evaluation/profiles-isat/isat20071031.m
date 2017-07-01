% ION-SATURATION profile
% FLUCTUATION-profile

% LOAD OFFSETS
  load('../offsets.mat');

% LOAD RAW DATA
  fls=dir('B2*');
  num = length(fls);

  for i=1:num
    [A tt] = readmdf(fls(i).name);
    b = A(:,1)-mean(A(:,1));
    fluc(i) = mean(abs(b));
    IsDC(i) = mean(A(:,2));
    int_i = mean(A(:,6));
  end;

  int = mean(int_i);
% REMOVE Isat OFFSET
  IsDC = IsDC - isoff;
% REMOVE INTERFEROMETER OFFSET
  int = int - intoff;

% <<INPUT>> CREATE r-VECTOR
  r = (-6:0.5:6)/100;

% CALCULATE n-PROFILE
  [nprof.n_x nprof.n_y nprof.is_x nprof.is_y nprof.exp10 alpha] = isprof2nprof(r, IsDC, int);

% PLOT ISAT
  fonts = 14;
figeps(12,8,1);
  hold on
  plot(nprof.is_x, nprof.is_y, 'ko');
  plot(nprof.n_x, nprof.n_y, 'r');
  hold off
  xlabel('radial position [m]');
  ylabel(['n [10^{' num2str(nprof.exp10) '} m^{-3}]']);
  box on;
  set(gca, 'fontsize', fonts);
  % PRINT ISAT
  ftree = fnametree(pwd);
  print('-depsc2', [ftree{2} ftree{1} '_isatprof.eps']);

% PLOT FLUCTUATION LEVEL
figeps(12,8,2);
  flucprof.r = r;
  flucprof.y = fluc;
  plot(flucprof.r, flucprof.y, '-o');
  xlabel('radius [cm]');
  ylabel('\delta n [a.u.]');
  % PRINT fluc-profile
  ftree = fnametree(pwd);
  print('-depsc2', [ftree{2} ftree{1} '_flucprof.eps']);

% SAVE TO mat-file
fn = ['isat_' ftree{2} ftree{1} '.mat'];
save(fn, 'nprof','fluc');