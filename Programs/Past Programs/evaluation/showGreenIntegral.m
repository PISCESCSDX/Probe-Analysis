t = (0:10000)/10;

% Noise level
noise = 0.3;

% Background constant
co = 0;

% Uncorrelated signals
x1 = co + noise*randn(numel(t),1);
y1 = co + noise*randn(numel(t),1);
[~, xyint1] = integratexy(x1,y1);

% Correlated signals
phi = +1.0*pi/2;
x2 = co + noise*randn(numel(t),1) + sin(t');
y2 = co + noise*randn(numel(t),1) + sin(t'+phi);
[~, xyint2] = integratexy(x2,y2);

figeps(18,22,1); clf; 
subplot(3,2,1)
  hold on
  plot(x1,'k');
  plot(y1,'b'); set(gca,'xlim',[0 1000])
  hold off
  title('uncorrelated signals')
subplot(3,2,2)
  hold on
  plot(x2,'k');
  plot(y2,'r'); set(gca,'xlim',[0 1000])
  hold off
  title('correlated signals')
subplot(3,2,3)
  plot(x1,y1,'x')
  title('phase space: x vs y')
subplot(3,2,4)
  plot(x2,y2,'rx')
  title('phase space: x vs y')
subplot(3,2,5)
  plot(cumsum(xyint1))
  set(gca,'xlim',[t(1) t(end)])
  title('cumulated Green integral')
subplot(3,2,6)
  hold on
  plot(cumsum(xyint1),'b')
  plot(cumsum(xyint2),'r')
  hold off
  set(gca,'xlim',[t(1) t(end)])
  title('cumulated Green integral')