% MAT-FILE NAME
 fname = 'isat_mex_dw_z2_off.mat';
 ylb = 'i_{sat} [mA]';

% STEPWIDTH of r-Scan [m]
  dx = 0.005;

% LOAD U-I-Probe-Characteristics in VAR 'data' and 'I_sat'
  a = dir('meas*.dat');
  aend = length(a);
  for i=1:aend
    data{i} = load_raw_data(a(i).name, '%--', ' ');
    I_sat(i) = mean( data{i}(10:100,2) );
  end
  
% PLOT I_sat
  fig(12,8,1);
  % stepwidth in [cm]
  dx = dx*100;
  rv = ((1:aend)-21)*dx;
  % I_sat in mA
  ImA = I_sat*1e3;
  subplotprof(rv, ImA, ylb);

save(fname, 'rv', 'ImA', 'ylb');