function ShowFFTDec(basedir, ipic, dt)
%==========================================================================
%function ShowFFTDec(basedir, ipic, dt)
%--------------------------------------------------------------------------
% Formerly known as 'mkfftdecpic.m'
%
% ipic: e.g. 1:10
% dt: e.g. 0.2 OR -1 for manual pressing button
%------------------------------------------------------------------ EXAMPLE
% mkfftdecpic('18680_f650_ap1.4', 1:400, 0.5)
%==========================================================================


%==========================================================================
% PLOTS
%==========================================================================
figx = 24;
figy = 9.375;
figeps(figx,figy,1,0,5);
clf
round(6.75*figx); % in mm
fonts = 12;

% distance constants
x0 = 0.070;
dx = 0.116;
xw = 0.100;

yw0= 0.30; yw1= 0.30;
y0 = 0.58; y1 = 0.14;

% Define maximal mode number plotted
m_max = 6;

% Plot variables (depends on normalization of En_m and En_mr)
% Determine from mkfftdecpic_help
ylim{1} = [  0  30];
ylim{2} = [  0  20];
ylim{3} = [  0  30];
ylim{4} = [  0  10];
ylim{5} = [  0   3];
ylim{6} = [  0   2];
ylim{7} = [  0   2];
ylim{8} = [  0   2];
cl = 1.0*[-1 1];
cl0 = 1.0*cl.*[1 1];

%==========================================================================
% MODE PLOTS
%==========================================================================
% Define color limits (light fluctuation in percent/100)
xytick = [-50 0 50];

info = cineInfo([basedir '.cine']);

load([basedir '_statistics.mat'])
% Normalize Modes to average amplitude in percent (='AAP')
pavg = matmean(movstat.avg) / ...
  mean( (movstat.lightfluc.max-movstat.lightfluc.min)/2 );


for i = ipic
%   disp_num(i,la)
  % Load plot data
  fn = mkstring([basedir '\mkfftdec_'], '0', i, info.NumFrames, '.mat');
  load(fn);
  % Plot camera measurement (removed offset)
  axes('position', [x0+0*dx y0 xw yw0])
  pcolor(mode.cam.x, mode.cam.y, mode.cam.pic/pavg);
  axis square
  shading interp
  colormap(pastell(64))
  %colormap(fireice(64))
  %set(gca, 'xlim', 5.1*[-1 1], 'ylim', 5.1*[-1 1], 'clim', cl);
  set(gca, 'clim', cl0);
  set(gca, 'xtick', xytick, 'ytick',  xytick)
  mkplotnice('x (mm)', 'y (mm)', fonts, '-20', '-30');
  str = ['t=' sprintf('%0.2f',(i-1)/mode.fs*1000) 'ms'];
  puttextonplot(gca, [0 1], 5, -15, str, 0, 12, 'k');
  mknicecolorbar('NorthOutside', ...
    '\Delta\Gamma_\lambda/\Gamma_{\lambda{0}}',fonts,0.15,0.1,3);

  % PLOT mode decomposed plots
  for im=1:m_max+1
    axes('position',[x0+im*dx y0 xw yw0])
    if im<m_max+1
      pcolor(mode.m{im}.x, mode.m{im}.y, mode.m{im}.pic/pavg) % ??? WHy *AAP
    else
      if im>1
        % inverse FFT
        Q = 0;
        for l=1:m_max-1
          Q = Q + mode.m{l}.pic/pavg;
        end
        pcolor(mode.m{im}.x, mode.m{im}.y, Q)
      else
        pcolor(mode.m{im}.x, mode.m{im}.y, mode.m{im}.pic/pavg)
      end
    end
    axis square
    shading interp
    colormap(pastell(64))
    % colormap(fireice(64))
    set(gca, 'xtick', xytick, 'ytick', xytick)
    if im>1
      if im==m_max
        set(gca, 'clim', cl0);
      else
        set(gca, 'clim', cl);
      end
    else
      set(gca, 'clim', cl);
    end
    mkplotnice('x (cm)', '-1', fonts, -20);
    if im<m_max+1
      puttextonplot(gca,[0 1],35,35,['m=' num2str(mode.mvec(im))],0,fonts,'k');
    else
    % Superimposed plot (inverse FFT)
      puttextonplot(gca,[0 1],5,15,['ifft: m0 to ' num2str(m_max-1)],0,fonts,'k');
    end
  end



%==========================================================================
% RADIAL LIGHT INTENSITY PROFILES
%==========================================================================
axes('position', [x0+0*dx y1 xw yw1])
plot(radprof.r, radprof.En_m, 'k', 'Linewidth', 1.5);
set(gca, 'xlim', [0 radprof.r(end)])
set(gca, 'xtick', 10*[2 4 6 8])
set(gca, 'ylim', ylim{1});
mkplotnice('r (mm)', 'A^2 (arb.u.)', fonts, '-20', '-30');
%puttextonplot(gca, [0 0], 10, 92, str, 0, fonts-2, 'k');
%

col = hsv(m_max);
for im=1:m_max
  axes('position', [x0+im*dx y1 xw yw1])
  if im<length(mode.mvec)+1
    plot(radprof.r, radprof.En_mr{im}, 'Color', col(im,:), 'Linewidth', 2);
  else
    hold on
    plot(radprof.r, radprof.En_m,   'Color', 'r', 'Linewidth', 2);
    plot(radprof.r, radprof.En_cam, 'Color', 'k', 'Linewidth', 2);
    hold off
  end
  % Next row: radprof.m<im>{i} = En_mr{im}*rprof_norm;
  %eval(['radprof.m' num2str(im) '{' num2str(i) '} = En_mr{' num2str(im) '}*rprof_norm;'])
  set(gca, 'ylim', ylim{im+1});
  set(gca, 'xlim', [0 radprof.r(end)])
  set(gca, 'xtick', 10*[2 4 6 8])
  mkplotnice('r (mm)', '-1', fonts, '-22');
  if im==length(mode.mvec)+1
    puttextonplot(gca, [0 0], 5, 80, str, 0, fonts-2, 'k');
    puttextonplot(gca, [1 1], -87, -8, 'ifft: m{\rightarrow}6', 0, fonts-2, 'r');
    puttextonplot(gca, [1 1], -30, -8, '{\rightarrow}32', 0, fonts-2, 'k');
  end
end


%save i info cl xytick basedir x0 dx y0 xw yw0 fac AAP fonts m_max y1 yw1 ylim dt
clear mode sm lightfluc
%load i info cl xytick basedir x0 dx y0 xw yw0 fac AAP fonts m_max y1 yw1 ylim dt

drawnow

%==========================================================================
% PRINT
%==========================================================================
% fname = ['fftdec_' num2str(i) '.eps'];
% print_adv([0 0 0 0 0 0 0 0  1 1 1 1 1 1 1 0 1], '-r300', fname, 50, 4);

% return
%fname = ['fftdec_' num2str(i) '.jpg'];
% print('-djpeg', '-r300', fname);


if dt==-1
  input('press any button to show next picture')
else
  pause(dt)
end
clf


%==========================================================================

end

% clear all; clf

end