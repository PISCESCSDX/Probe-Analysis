% EXAMPLE for cross-correlation

% measured data
  [A tt] = readmdf('cou0002.MDF');
  tt = ((1:2^16)'-1)*800e-9;
  dt = tt(2)-tt(1);
  fs = 1/dt;
% artificial noisy signal
  a=0.1*randn(size(tt,1), 1);
  b=0.1*randn(size(tt,1), 1);
% phase shifted signal
  phi = 180;
  phsh = (phi/360)*2*pi;
  aint = 30000:35000;
  a(aint) = sin(2*pi*5e3*tt(aint));  
  aint = 31000:36000;  
  b(aint) = sin(2*pi*5e3*tt(aint) + phsh);    
  %
  nA = 0.1;
  p1 = sin(2*pi*5e3*tt) + nA*a;
  p2 = sin(2*pi*5e3*tt+phsh) + nA*b;

  p1 = A(:,1);
  p2 = A(:,1);

% cross correlation function of p1 and p2
 [cpp, tau] = ccf(p1, p2, fs);
%   [cpp, tau] = ccf(a, b, fs);

% PLOT
  plot(tau, cpp);
  
  figure
  plot(tau(1:400), p1(1:400), 'b')
  hold on
  plot(tau(1:400), p2(1:400), 'r')
% save corrtest.mat corr  