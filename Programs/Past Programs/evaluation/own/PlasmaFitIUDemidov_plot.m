figure(2); clf
axes('Position',[0.13 0.55 0.80 0.37])
ax{1} = get(gca,'Position');
  % Dummy axis
  set (gca, 'Box', 'on', 'Color', 'white', 'XTick', [], 'YTick', []);
ax2 = axes('position', ax{1});
  set(ax2, 'XAxisLocation','top')
  set(ax2, 'YAxisLocation','right')
  hold on
  hand.line = plot(xind, 1e3*ye, 'o-', 'Color', rgb('NavajoWhite'));
  hold off
  set(ax2,'xlim',[xind(1) xind(end)])
  xlabel('sample #','FontSize',12)
  set(gca,'FontSize',12)
ax1 = axes('position', ax{1});
  hold on
  plot(x          ,1e3*ye   ,'o-','Color', rgb('NavajoWhite'))
  plot(x(xifitres),1e3*Iefit,'r-','LineWidth', 2)
  hold off
  set(ax1,'xlim',[x(1) x(end)])
  set(ax1,'color','none')
  set(hand.line, 'visible', 'off')

lx = 10; ly=120; dy=-20;
d = max(ye)-min(ye);
set(ax1, 'ylim', 1e3*[min(ye)-0.02*d max(ye)+0.02*d])
set(ax2, 'ylim', 1e3*[min(ye)-0.02*d max(ye)+0.02*d])
mkplotnice('probe voltage (V)', 'fit I_e (mA)', 12, '-25', '-40');
ilb=-1;
ilb=ilb+1; puttextonplot(gca, [0 0], lx, ly+ilb*dy, ...
                ['R=' sprintf('%4.2f', iudata.RsqrIe)], 0, 12, 'k');
ilb=ilb+1; puttextonplot(gca, [0 0], lx, ly+ilb*dy, ...
                ['T_e=' sprintf('%4.2f', iudata.Te) 'eV'], 0, 12, 'k');
ilb=ilb+1; puttextonplot(gca, [0 0], lx, ly+ilb*dy, ...
                ['n_e=' sprintf('%4.2g', iudata.ne) 'm^{-3}'], 0, 12, 'k');
ilb=ilb+1; puttextonplot(gca, [0 0], lx, ly+ilb*dy, ...
                ['V_p=' sprintf('%4.2f', iudata.Vp) 'V'], 0, 12, 'k');
ilb=ilb+1; puttextonplot(gca, [0 0], lx, ly+ilb*dy, ...
                ['V_f=' sprintf('%4.2f', iudata.Vf) 'V'], 0, 12, 'k');
ilb=ilb+1; puttextonplot(gca, [0 0], lx, ly+ilb*dy, ...
          ['\Gamma=' sprintf('%4.2g',iudata.Fli) 'm^{-2}s^{-1}'],0,12,'k');

axes('Position',[0.13 0.09 0.80 0.37])
hold on
  plot(x,1e3*Dye, 'o-', 'Color', rgb('NavajoWhite'))
  yf = 1e3*ffun([fit.ne fit.Te], x(xifitres));
  plot(x(xifitres), yf, 'r-', 'LineWidth', 2)
hold off
d = max(yf)-min(yf);
set(gca, 'xlim', [x(1) x(end)], 'ylim', [min(yf)-0.5*d max(yf)+0.5*d])
mkplotnice('probe voltage (V)', 'dI_e/dV (mA/V)', 12, '-25', '-40');
puttextonplot(gca, [0 1], 10, -35, 'fit to dI_e/dV', 0, 12, 'r');
puttextonplot(gca, [0 0], 10, 90, ...
                ['R=' sprintf('%4.2f', iudata.RsqrDIe)], 0, 12, 'k');