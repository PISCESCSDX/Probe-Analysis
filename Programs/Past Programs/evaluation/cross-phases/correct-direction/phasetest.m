fsample = 2.0000e+05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tn = (0:2^16-1)/fsample;
nois = 0.1;
dphi = -2*pi/3;

x = sin(2*pi*8000*tn);
Tn = x + nois*randn(length(tn),1)';

x = sin(2*pi*8000*tn + dphi);
T1 = x + nois*randn(length(tn),1)';

ind=(1:100);

figeps(8,6,1)
hold on
  plot(tn(ind), Tn(ind), 'b')
  plot(tn(ind), T1(ind), 'r')
hold off
mkplotnice('time','signals',12);
ht = puttextonplot(gca, 5, 90, 'blue: Tn~sin(a)' ); set(ht, 'fontsize',10);
ht = puttextonplot(gca, 5, 80, 'red:  T1~sin(a - pi/2)' ); set(ht, 'fontsize',10);
print('-depsc2', '-r300', 'signals.eps')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the cross-phase
n = 100;      % n=100 is ok
olap = 0.5;
[A, ph, f, phstd] = cpsdmean(Tn, T1, fsample, n, olap, 60e3);

[h1 h2] = plotfftphase(f, 20*log10( abs(A) ), ph, phstd, [-170 -60]);
ht = puttextonplot(gca, 5, 90, 'cpsdmean(Tn,T1)'); set(ht, 'fontsize',10);
print('-depsc2', '-r300', 'cpsdmean.eps');