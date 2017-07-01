%==========================================================================
% piscesprobe_v2012-04-13 14:15
% Evaluation of retracted probe measurements at PISCES-A
%--------------------------------------------------------------------------
%==========================================================================


%==========================================================================
% Parameters: Probe, Probe Circuit, DAQ & Plasma
%--------------------------------------------------------------------------
% Probe Paramters
%--------------------------------------------------------------------------
% Total A (m^2) is given, but we need r, l
  probe.A = 6.48e-6;
% assume A is calculated as: A = 2*pi*r*l + pi*(r^2)
% assume r (m)
  probe.r = 0.040/2 * 2.54/100;
% calculate length (m)
  probe.l = (probe.A - pi*(probe.r^2)) / (2*pi*probe.r);
%--------------------------------------------------------------------------
% Probe Circuit & DAQ Parameters
%--------------------------------------------------------------------------
iu.vol.fac = 1/100; % Voltage attenuation
iu.cur.Rsh = 101.0; % Probe shunt resistor
iu.cur.fac = 1/1;   % Attenuation at LeCroy Differential Amplifier
iu.pos.fac = 50.8;  % Position Factor [priv. comm. Ben Hudon]
fs = 100e3;         % Sample frequency
%--------------------------------------------------------------------------
% Plasma parameters
%--------------------------------------------------------------------------
% Magnetic field on the axis
  I_Sour =  50;
  I_Trim =   0;
  I_Main = 400;
  B = B_PISCES(I_Sour, I_Trim, I_Main);
% Ion Mass [u0]
% m_ion = 2.01410178; % Deuterium
  m_ion = 4.002602; % Helium
%==========================================================================




%==========================================================================
% load file list of raw data files
%--------------------------------------------------------------------------
a = dir('*.raw');
la= length(a);
%==========================================================================




%==========================================================================
% Type of radial fits: 1: gaussian fits, 2: average smoothing
% # Input file number
num =    input('Shot Number                          : ');
%num = 7;

% InPlot = input('Calculate(Empty)   Just Plot(1)      : '); % ***
InPlot = [];
if isempty(InPlot)                                              % IF InPlot

% Ask if Check Plots shall be shown
chk =    input('Show Check Plots?   No(empty) Yes(1) : ');

%--------------------------------------------------------------------------
% Display evaluation summary
%--------------------------------------------------------------------------
disp('==================================')
disp('PISCES probe evaluation summary')
disp('----------------------------------')
disp(['rawdata file: ' a(num).name])
disp(['I_Sour = ' sprintf('%.0f', I_Sour) ' A'])
disp(['I_Trim = ' sprintf('%.0f', I_Trim) ' A'])
disp(['I_Main = ' sprintf('%.0f', I_Main) ' A'])
disp(['     B = ' sprintf('%.0f', B*1e3) ' mT'])
disp(['   L_p = ' sprintf('%.1f', probe.l*1e3) ' mm'])
disp(['   R_p = ' sprintf('%.1f', probe.r*1e3) ' mm'])
disp(['   A_p = ' sprintf('%.2f', probe.A*1e6) ' mm^2'])
disp(['   f_s = ' sprintf('%.0f', fs/1e3) ' kHz'])
disp('==================================')
%--------------------------------------------------------------------------


%==========================================================================
% 1. Extraction of voltage, current, and radial position
%--------------------------------------------------------------------------
% Using LabVIEW vi karl_probe one mesaurement yields
% two files: <filename>.raw and <filename>.hea
% The header file lines:
% # 1: number of points (e.g. 100000)
% # 2: scan rate (Hz)   (e.g. 100000)
% The raw file contains pure raw data, i.e. not changed by LabVIEW
%  col 1: 1/100 probe voltage (V)
%  col 2: voltage at probe shunt (V)
%         - measurement of shunt voltage at 101 Ohm shunt
%         - differential amplifier: Read Effective Gain
%  col 3: N/A (Nov 2011 - April 2012)
%  col 4: Position
  fname = a(num).name;
  dat = load(fname);
% Calculate probe voltage (V)
  iu.vol =  dat(:,1) / iu.vol.fac;
% Remove Offset
% Position of probe outside of vessel: over 2.5 V
  ind = dat(:,4)>2.5;
  h = dat(ind,2);
  offset = mean(h(100:end-100));
% calculate probe current (A)
  iu.cur = ((dat(:,2)-offset) / iu.cur.fac) / iu.cur.Rsh;
% calculate radial position
  iu.pos = dat(:,4) * iu.pos.fac;
% define the time vector
  iu.tvec    = (0:length(iu.vol)-1) * (1/fs);
%==========================================================================




%==========================================================================
% 2. Extraction of single characteristics
%--------------------------------------------------------------------------
% Procedure to get the single ramps:
% 1 vsm = smooth voltage trace
% 2 diff = differentiate vsm
% 3 diffsm = smooth diff
% 4 find zero crossings of diffsm
% # Input: smooth value probe voltage
  smooth_vol = 100;
  vd = diff_discrete(iu.tvec', filtmooth(iu.vol,smooth_vol));
% # Input: smooth value: derivation of probe voltage
  smooth_vdiff = 10;
  vdsm = filtmooth(vd, smooth_vdiff);
  vdsmsig = sign(vdsm);
% Find zero crossings: they are positions between two ramps
  vec1 = vdsmsig;
  vec2 = [0; vdsmsig(1:end-1)];
  vecdiff = vec2-vec1;
% Find positions where vecdiff is > 0 -> ramp transitions
  tpos = find(abs(vecdiff) > 0);
%==========================================================================



%==========================================================================
% Extraction of single characteristics and of average position
%--------------------------------------------------------------------------
tpos(1) = 1;
rvec = zeros(1,length(tpos)-1);
ctr = 0;
for i=2:length(tpos)
  ctr = ctr + 1;
  ind = tpos(i-1):tpos(i);
  % Flip to positive ramp direction
  if vecdiff(tpos(i)) < 0
    ind = fliplr(ind);
  end
  
  % Sort function + Averaging of multiple data points
  [v1 v2] = XYAVGMultiples(iu.vol(ind), iu.cur(ind));
  pr.vol{i-1} = v1';
  pr.cur{i-1} = v2';
  
  % Average position of current characteristic
  rvec(ctr) = mean(iu.pos(ind));
end
% Save r-vector to file
fn  = a(num).name;
savefn = [fn(1:end-4) '_rvec.mat'];
save(savefn, 'rvec');
%==========================================================================




%>>>===============================
% Video of  characteristics
%----------------------------------
% figeps(12,8,1,30,30);
% xlim = [-150 50];
% ylim = 1e3*[-0.3 0.05];
% for i=1:length(pr.cur)
%   clf
%   axes('position', [0.18 0.18 0.75 0.75])  
%   plot(pr.vol{i}, 1e3*pr.cur{i}, 'x')
%   set(gca, 'xlim', xlim, 'ylim', ylim)
%   mkplotnice('voltage (V)', 'current (mA)', 12, '-30', '-30');
%   puttextonplot(gca, [0 1], 5, -15, num2str(i), 0, 12, 'k');
%   input('<<< press any key >>>')
% end
%<<<===============================





%>>>=============================================
% Save extracted single characteristics to files
%------------------------------------------------
fn  = a(num).name;
% Delete Old Saved Characteristics:
  delete( [fn(1:end-4) '_iu*.dat'] );
for i=1:length(pr.cur)
  % Filename of each characteristic
  savefn = mkstring([fn(1:end-4) '_iu'],'0',i,length(pr.cur),'.dat');
  % Save characteristic
  A = [pr.vol{i} pr.cur{i}];
  %save(savefn, '-ascii', '-double', '-tabs', 'A')
  dlmwrite(savefn, A, 'delimiter', '\t', 'precision', 6);
end
%<<<=============================================





disp('Calculating plasma paramters ...')
%>>>=======================================================================
% Evaluation of single characteristics
%--------------------------------------------------------------------------
% This m-file is for radial IV-evaluation using 'PlasmaFitIULangmuir'

% find Pisces-A rawdata file
rawbase = a(num).name(1:end-4);
f = dir([rawbase '_iu*.dat']);
lf= length(f);

% Preallocating radial plasma parameters
ne.n = zeros(lf,1);   ne.e = zeros(lf,2);
Te.n = zeros(lf,1);   Te.e = zeros(lf,2);
Vp.n = zeros(lf,1);   Vp.e = zeros(lf,2);
Vf.n = zeros(lf,1);   Vf.e = zeros(lf,2);
R = zeros(lf,1);

%==========================================================================
% Calculate plasma parameters from I-V-characteristic
%--------------------------------------------------------------------------
% num = input('Num #: ');
for i=1:length(f)
  disp(['# ' num2str(i)])
  A = load(f(i).name);
% Assign input variables
  data.voltage = A(:,1);
  data.current = A(:,2);
  data.r = probe.r;
  data.l = probe.l;
  data.m_ion = m_ion;
  data.B = B;
% Show plots of the fits
%   chk = 0;
% Set the start parameters and the fit boundaries
  fixpar.Te    = 5;
  fixpar.Telim = [0.5 50];
  fixpar.ne    = [];
  fixpar.nelim = [];
  fixpar.Vp    = [];
% Use whole characteristic (up to V_p)
%  Vint = [-Inf, NaN];
  Vint = [-Inf, NaN];
% # Input: Smooth over 'vol_avg' Volts
  vol_avg = 1;
% # INPUT: Level-Range of Slope (parameter with strong influence on T_e!!)
  p_IexpLevel = 0.95;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iudata = PlasmaFitIUDemidov_v10(data, chk, fixpar, Vint, vol_avg, ...
  p_IexpLevel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================


%==========================================================================
% Store variables to radial positions
%--------------------------------------------------------------------------
if isempty(iudata)
  ne.n(i,1) = NaN;
  ni.n(i,1) = NaN;
  Te.n(i,1) = NaN; Te.e(i,:) = NaN;
  Vp.n(i,1) = NaN; Vp.e(i,:) = NaN;
  Vf.n(i,1) = NaN; Vf.e(i,:) = NaN;
  Vf2.n(i,1)= NaN; Vf2.e(i,:) = NaN;
  Fli.n(i,1)= NaN;
  R(i,1)    = NaN;
else
  ne.n(i,1) = iudata.ne;
  ni.n(i,1) = iudata.ni;
  Te.n(i,1) = iudata.Te; Te.e(i,:) = iudata.Teerr;
  Vp.n(i,1) = iudata.Vp; Vp.e(i,:) = iudata.Vperr;
  Vf.n(i,1) = iudata.Vf;
  Vf2.n(i,1)= iudata.Vf2;
  Fli.n(i,1)= iudata.Fli;
  R(i)      = iudata.Rsqr;
end

end
%==========================================================================



% # Input: Manually set the center of the radial vector
rcenter = 70;
rvec = (rvec-rcenter)/10;


%==========================================================================
% Calculate Averages of the radial profiles
%--------------------------------------------------------------------------
% Electron density profile
%--------------------------
[~, ir] = sort(rvec);
x = rvec(ir); y = ne.n(ir);
ind = ~isnan(y);
x = x(ind); y = y(ind);
%
nefit.x = min(x):0.1:max(x);
nefit.y = csaps(x, y, 0.90, nefit.x);
%--------------------------
% Ion density profile
%--------------------------
[~, ir] = sort(rvec);
x = rvec(ir); y = ni.n(ir);
ind = ~isnan(y);
x = x(ind); y = y(ind);
%
nifit.x = min(x):0.1:max(x);
nifit.y = csaps(x, y, 0.90, nifit.x);
%--------------------------
% Floating potential
%--------------------------
[~, ir] = sort(rvec);
x = rvec(ir); y = Vf.n(ir);
ind = ~isnan(y);
x = x(ind); y = y(ind);
%
%[x,y,~] = FctDeleteOutlier(x,y,0.15);
Vffit.x = min(x):0.1:max(x);
Vffit.y = csaps(x, y, 0.9, Vffit.x);
%------ Vf2 (from fit) -------
x = rvec(ir); y = Vf2.n(ir);
ind = ~isnan(y);
x = x(ind); y = y(ind);
%
Vf2fit.x = min(x):0.1:max(x);
Vf2fit.y = csaps(x, y, 0.9, Vf2fit.x);
%--------------------------
% Plasma potential
%--------------------------
[~, ir] = sort(rvec);
x = rvec(ir); y = Vp.n(ir);
ind = ~isnan(y);
x = x(ind); y = y(ind);
%
%[x,y,~] = FctDeleteOutlier(x,y,0.15);
Vpfit.x = min(x):0.1:max(x);
Vpfit.y = csaps(x, y, 0.9, Vpfit.x);
%--------------------------
% Temperature profile
%--------------------------
x = rvec(ir); y = Te.n(ir); ey = Te.e(ir);
ind = ~isnan(y);
x = x(ind); y = y(ind); ey = ey(ind);
%
Tefit.x = min(x):0.1:max(x);
Tefit.y = csaps(x, y, 0.98, Tefit.x, 1./ey');
%--------------------------
% Ion Flux Profile
%--------------------------
x = rvec(ir); y = Fli.n(ir);
ind = ~isnan(y);
x = x(ind); y = y(ind);
%
Flifit.x = min(x):0.1:max(x);
Flifit.y = csaps(x, y, 0.90, Flifit.x);
%==========================================================================




%==========================================================================
% Plot n_e, V_p, V_f, T_e
%--------------------------------------------------------------------------
% Save radial profiles
save([rawbase '_radprof.mat'],'rvec','ne','ni','Vf','Vp','Te','Fli',...
  'R', 'nefit','nifit','Vffit','Vpfit','Tefit','Flifit')

end % End of InPlot


%==========================================================================
rawbase = a(num).name(1:end-4);
load( [rawbase '_radprof.mat'] );
%--------------------------------------------------------------------------

fonts = 12;
%--------------------------------------------------------------------------
% Use these lines for: 2 x 2 PLOTS
%----------------------------------
zoom = 1.2;
figx = zoom*12; figy = zoom*12;
figeps(figx,figy,3,47,100); clf;

x0 = 0.12; y0 = 0.10;
dx = 0.50; dy = 0.33;
%xw = 0.75;
yw = 0.22;
% Golden Cut:
xw = 1.618*yw*figy/figx;
ax{1} = [x0+0*dx y0+2*dy xw yw];
ax{2} = [x0+1*dx y0+2*dy xw yw];
ax{3} = [x0+0*dx y0+1*dy xw yw];
ax{4} = [x0+1*dx y0+1*dy xw yw];
ax{5} = [x0+0*dx y0+0*dy xw yw];
ax{6} = [x0+1*dx y0+0*dy xw yw];
% for l=1:6; axes('position', ax{l}); end
%----------------------------------
% Use these lines for: 2 x 2 PLOTS
%----------------------------------
% zoom = 1.0;
% figx = zoom*12; figy = zoom*20;
% figeps(figx,figy,3,47,100); clf;
% 
% x0 = 0.18; y0 = 0.08;
% dx = 0.50; dy = 0.31;
% %xw = 0.75;
% yw = 0.28;
% % Golden Cut:
% xw = 1.618*yw*figy/figx;
% ax{1} = [x0+0*dx y0+2*dy xw yw];
% ax{2} = [x0+0*dx y0+1*dy xw yw];
% ax{3} = [x0+0*dx y0+0*dy xw yw];
%========================================================================

% MarkerSize
MSz = 3;
% LineWidth of Error bars
errwdth = 0.1;

xlim = 4.1*[-1 1]; xtick = -6:1:6;
j = 0;

% --------- ni ---------
j = j+1; axes('position', ax{j})
mcol = [0 0 0];
log10max = floor( log10( max(ni.n) ) )-1;
hold on
  plot(rvec, ni.n/(10^(log10max)), 'o', 'MarkerSize', MSz, ...
  'MarkerFaceColor', mcol, 'MarkerEdgeColor', 'w', 'LineWidth', 0.3)
  dx = [];
  dy = (10^(-log10max))* 1*abs([ni.n.*(1-R) ni.n.*(1-R)]);
  my_errorbar(rvec,ni.n/(10^(log10max)),dx,dy,errwdth, mcol)
  plot(nifit.x, nifit.y/(10^log10max), 'k')
hold off
ymax = 1.1*max(ni.n)/(10^log10max); [dy, y1] = plotticks(0, ymax, 4);
ylim = [0 ymax]; ytick = y1:dy:ymax;
set(gca,'xlim',xlim,'xtick',xtick,'ylim',ylim,'ytick',ytick)
mkplotnice('-1', ['n_i (10^{' num2str(log10max) '}m^{-3})'], ...
  fonts, '', '-26');

% --------- ne ---------
j = j+1; axes('position', ax{j})
mcol = [0 0 0];
log10max = floor( log10( max(ne.n) ) )-1;
hold on
  plot(rvec, ne.n/(10^log10max), 'o', 'MarkerSize', MSz, ...
  'MarkerFaceColor', mcol, 'MarkerEdgeColor', 'w', 'LineWidth', 0.3)
  dx = [];
  dy = (10^(-log10max))* 1*abs([ne.n.*(1-R) ne.n.*(1-R)]);
  my_errorbar(rvec,ne.n/(10^log10max),dx,dy,errwdth, mcol)
  plot(nefit.x, nefit.y/(10^log10max), 'k')
hold off
ymax = 1.1*max(ne.n)/(10^log10max); [dy, y1] = plotticks(0, ymax, 4);
ylim = [0 ymax]; ytick = y1:dy:ymax;
set(gca,'xlim',xlim,'xtick',xtick,'ylim',ylim,'ytick',ytick)
mkplotnice('-1', ['n_e (10^{' num2str(log10max) '}m^{-3})'], ...
  fonts, '', '-26');

% --------- Te ---------
j = j+1; axes('position', ax{j})
mcol = [1 0 0];
hold on
 plot(rvec, Te.n, 'ro', 'MarkerSize', MSz, ...
  'MarkerFaceColor', mcol, 'MarkerEdgeColor', 'w', 'LineWidth', 0.6)
  dx = [];
  dy = abs(Te.e-[Te.n Te.n]);
  my_errorbar(rvec, Te.n, dx, dy, errwdth, mcol)
  plot(Tefit.x, Tefit.y, 'r')
hold off
ymax = 1.1*max(Te.n); [dy, y1] = plotticks(0, ymax, 4);
ylim = [0 ymax]; ytick = y1:dy:ymax;
set(gca,'xlim',xlim,'xtick',xtick,'ylim',ylim,'ytick',ytick)
mkplotnice('r (cm)', 'T_e (eV)', fonts, '-25', '-26');

% --------- Ion Flux ---------
j = j+1; axes('position', ax{j})
mcol = [0 0 0];
log10max = floor( log10( max(Fli.n) ) )-1;
hold on
  plot(rvec, Fli.n/(10^(log10max)), 'o', 'MarkerSize', MSz, ...
  'MarkerFaceColor', mcol, 'MarkerEdgeColor', 'w', 'LineWidth', 0.3)
  dx = [];
  dy = (10^(-log10max))* 1*abs([Fli.n .*(1-R) Fli.n .*(1-R)]);
  my_errorbar(rvec,Fli.n/(10^(log10max)),dx,dy,errwdth, mcol)
  plot(Flifit.x, Flifit.y/(10^(log10max)), 'k')
hold off
ymin = 0; ymax = 1.1*max(Fli.n/(10^log10max));
[dy, y1] = plotticks(ymin, ymax, 4);
ylim = [ymin ymax]; ytick = y1:dy:ymax;
set(gca,'xlim',xlim,'xtick',xtick,'ylim',ylim,'ytick',ytick)
mkplotnice('-1', ['\Gamma_i (10^{' num2str(log10max) '}m^{-2}s^{-1})'], ...
  fonts, '', '-26');

% --------- Vf ---------
j = j+1; axes('position', ax{j})
hold on
plot(rvec, Vf.n, 'o', 'MarkerSize', MSz, ...
  'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'w', 'LineWidth', 0.3)
plot(Vffit.x, Vffit, 'Color', 'c')
hold off
Mmin = matmin([Vf.n Vp.n]); Mmax = matmax([Vf.n Vp.n]); dM = Mmax-Mmin;
ymin = Mmin - 0.05*dM; ymax = Mmax + 0.05*dM;
[dy, y1] = plotticks(ymin, ymax, 4);
ylim = [ymin ymax]; ytick = y1:dy:ymax;
set(gca,'xlim',xlim,'xtick',xtick,'ylim',ylim,'ytick',ytick)
mkplotnice('r (cm)', 'V_f (V)', fonts, '-25', '-28');
% --------- Vp ---------
j = j+1; axes('position', ax{j})
mcol = [0 0 1];
hold on
plot(rvec, Vp.n, 'o', 'MarkerSize', MSz, ...
  'MarkerFaceColor', mcol, 'MarkerEdgeColor', 'w', 'LineWidth', 0.3)
  dx = [];
  dy = 10*abs(Vp.e-[Vp.n Vp.n]);
  my_errorbar(rvec,Vp.n,dx,dy,errwdth, mcol)
plot(Vpfit.x, Vpfit, 'Color', mcol)
hold off
set(gca,'xlim',xlim,'xtick',xtick,'ylim',ylim,'ytick',ytick)
mkplotnice('r (cm)', 'V_p (V)', fonts, '-25', '-28');

%--------------------------------------------------------------------------
% Number String
num_str = mkstring('','0',num,1000,'');
% Print as *.eps
fn = ['d10_' num_str '_' num2str(vol_avg) '.eps'];
print('-depsc2', fn)
%==========================================================================