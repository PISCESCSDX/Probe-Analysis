%==========================================================================
% Calculate a Poincaré Plot
%--------------------------------------------------------------------------
N = 2^14;
dt = 5e-6;
t = (0:N-1)'*dt;
fs = 1/dt;

f1 = 3000; a1=0.5; ph1 = 0.5*pi;
f2 = 3400; a2=0.4; ph2 = -0.5*pi;
f3 = 2900; a3=0.7; ph3 = +0.5*pi;
noise = 0.1;

S1 = a1*sin(2*pi*f1*t + ph1) + noise*randn(N,1);
S2 = a2*sin(2*pi*f2*t + ph2) + noise*randn(N,1);
S3 = a3*sin(2*pi*f3*t + ph3) + noise*randn(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ### INPUT: fixed time length for CCF shown [s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ccf_lt = 50/1e3;
int_long  = ceil( 1.2*ccf_lt/dt );
int_short = ceil( 1.0*ccf_lt/dt );

% ### INPUT
ind = 1:int_long;

s1 = S1(ind,1);
s2 = S2(ind,1);
s3 = S3(ind,1);
 t = t(ind,1);

%==========================================================================
% Calculate Poincaré 3D and 2D plots
%--------------------------------------------------------------------------
figeps(10,6,1);
plot(t, s1, 'b-o')

% Calculate time delay automatically
tau_delay = -1;
% Define time interval
itint = 1:int_short;

% Calculate Variables for the 3D Poincaré plots
[S1, dtau1, acf1] = Poincare3D(s1, fs, tau_delay, itint);
[S2, dtau2, acf2] = Poincare3D(s2, fs, tau_delay, itint);
[S3, dtau3, acf3] = Poincare3D(s3, fs, tau_delay, itint);
%==========================================================================


%==========================================================================
% PLOTS
%----------------------------------------------- Auto-correlation functions
figeps(10,6,2);
plot(acf1.tau, acf1.cf, 'k-o')
figeps(10,6,3);
plot(acf2.tau, acf2.cf, 'r-o')
figeps(10,6,4);
plot(acf3.tau, acf3.cf, 'b-o')
%--------------------------------------------- Calculate the Poincaré plane
zeroplane = 3;
[px1, py1] = Poincare2D(S1.s1, S1.s2, S1.s3, zeroplane);
[px2, py2] = Poincare2D(S2.s1, S2.s2, S2.s3, zeroplane);
[px3, py3, Poincare2D_info] = Poincare2D(S3.s1, S3.s2, S3.s3, zeroplane);
%-------------------------------------------------------- 3D Poincaré plots
fonts = 12;
figx = 15; epsx = round(figx*7.5);
figeps(figx,6,5); clf;

y0 = 0.20;
ax{1} = [0.10 y0 0.28 0.70];
ax{2} = [0.63 y0 0.30 0.70];

alim = 1.2*[-1 1];
dist = 1;
ind = 1:dist:size(S1.s1,1);
msz = 4; % MarkerSize
%-------------------------------------------------------- Plot the signal 1
i=1;
axes('position', ax{i})
hp = plot3(S1.s1(ind), S1.s2(ind), S1.s3(ind), '.k', 'MarkerSize', msz);
set(hp, 'LineWidth', 0.1)
set(gca, 'xlim', alim , 'ylim', alim,'zlim', alim)
set(gca,'CameraPosition',  [+1 +2 +2])
[hxl hyl] = mkplotnice('s(t)', 's(t+\tau)', fonts, '-10', '0');
set(hxl, 'rotation', 10); set(hyl, 'rotation', -55);
zlabel('s(t+2\tau)', 'fontsize', fonts)
box off
set(gca, 'xgrid', 'on', 'ygrid', 'on')
puttextonplot(gca, [0 1], -40, -150, '(a)', 0, fonts, 'k');

i=2;
axes('position', ax{i})
plot(px1, py1, 'ko', 'markerfacecolor', 'k', 'markersize', 0.2);
mkplotnice('s(t)', 's(t+\tau)', fonts);
set(gca, 'xlim', alim, 'ylim', alim)
puttextonplot(gca, [0 1],  5, -15, '(b)', 0, fonts, 'k');
puttextonplot(gca, [0 1], 57, 8, Poincare2D_info, 0, fonts, 'k');
%-------------------------------------------------------- Plot the signal 2
figeps(figx,6,6); clf;
i=1;
axes('position', ax{i})
hp = plot3(S2.s1(ind), S2.s2(ind), S2.s3(ind), '.k', 'MarkerSize', msz);
set(hp, 'LineWidth', 0.1)
set(gca, 'xlim', alim , 'ylim', alim,'zlim', alim)
set(gca,'CameraPosition',  [+1 +2 +2])
[hxl hyl] = mkplotnice('s(t)', 's(t+\tau)', fonts, '-10', '0');
set(hxl, 'rotation', 10); set(hyl, 'rotation', -55);
zlabel('s(t+2\tau)', 'fontsize', 12)
box off
set(gca, 'xgrid', 'on', 'ygrid', 'on')
puttextonplot(gca, [0 1], -40, -150, '(a)', 0, fonts, 'k');

i=2;
axes('position', ax{i})
plot(px2, py2, 'ko', 'markerfacecolor', 'k', 'markersize', 0.2);
mkplotnice('s(t)', 's(t+\tau)', 12);
set(gca, 'xlim', alim, 'ylim', alim)
puttextonplot(gca, [0 1],  5, -15, '(b)', 0, fonts, 'k');
puttextonplot(gca, [0 1], 57, 8, Poincare2D_info, 0, fonts, 'k');
%-------------------------------------------------------- Plot the signal 3
figeps(figx,6,7); clf;
i=1;
axes('position', ax{i})
hp = plot3(S3.s1(ind), S3.s2(ind), S3.s3(ind), '.k', 'MarkerSize', msz);
set(hp, 'LineWidth', 0.1)
set(gca, 'xlim', alim , 'ylim', alim,'zlim', alim)
set(gca,'CameraPosition',  [+1 +2 +2])
[hxl hyl] = mkplotnice('s(t)', 's(t+\tau)', fonts, '-10', '0');
set(hxl, 'rotation', 10); set(hyl, 'rotation', -55);
zlabel('s(t+2\tau)', 'fontsize', 12)
box off
set(gca, 'xgrid', 'on', 'ygrid', 'on')
puttextonplot(gca, [0 1], -40, -150, '(a)', 0, fonts, 'k');

i=2;
axes('position', ax{i})
plot(px3, py3, 'ko', 'markerfacecolor', 'k', 'markersize', 0.2);
mkplotnice('s(t)', 's(t+\tau)', 12);
set(gca, 'xlim', alim, 'ylim', alim)
puttextonplot(gca, [0 1],  5, -15, '(b)', 0, fonts, 'k');
puttextonplot(gca, [0 1], 57, 8, Poincare2D_info, 0, fonts, 'k');
%================================================================ End PLOTS