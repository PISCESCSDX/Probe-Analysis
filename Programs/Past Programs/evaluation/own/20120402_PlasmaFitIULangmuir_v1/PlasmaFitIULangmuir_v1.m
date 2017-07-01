function [iudata] = PlasmaFitIULangmuir(data, checkonoff, fixpar)
%==========================================================================
% *** INFO: This m-file is a backup of the development. The version to use
%           is PlasmaFitIULangmuir.m!
%function [iudata] = PlasmaFitIULangmuir(data, checkonoff, fixpar)
%--------------------------------------------------------------------------
% PlasmaFitIULangmuir(data, checkonoff) calculates the plasma parameters 
% n_e (electron density, T_e (electron temperature), V_p (plasma potential)
% and the V_f (floating potential) from a current voltage characteristic
% using a cylindrical probe.
% Sub-m-files: subPlasmaIsatPoly, 
%--------------------------------------------------------------------------
% Assumptions for the applicability of the Langmuir theory.
% (1) The plasma is collisionless. This holds usually, if the pressure and
%     the denisity is low and the electron temperature is relatively high.
% (2) If the plasma is magnetized, the mean free path of the electrons is
%     replaced by the gyration radius and should be very large.
% Rule of thumb: R/lambda_e << 1 (R: probe radius, lambda_e: mean free 
% path of electrons or the gyration radius)
%--------------------------------------------------------------------------
% IN: data.voltage, data.current, data.l, data.r
%     checkonoff: 0: no plots, 1: check plots
%     fixpar.Te: fixed electron temperature
%     fixpar.ne: fixed electron density
%     fixpar.Vp: fixed plasma potential
%OUT: iudata.ne: electron density (m^-3)
%     iudata.Te: electron temperature (eV)
%     iudata.Vp: plasma potential (V)
%     iudata.Vf: floating potential (V)
%--------------------------------------------------------------------------
% EX: iudata = PlasmaFitIULangmuir(data, checkonoff, fixpar);
%==========================================================================


% Extract the characteristic and the probe parameters
x = data.voltage;
y = data.current;


% 1. Estimate the ion saturation current (use filter)
%--------------------------------------------------------------------------
data.voltlim = [];
[PolyIsat, voltlim] = subPlasmaIsatPoly(data);

if isempty(PolyIsat)
  disp('No results can be obtained from the IV-characteristic.')
  return
end


% XXX --->>> Check ion saturation current detection:
if checkonoff
  figeps(12,8,1,0,90);
  axes('position', [0.15 0.15 0.80 0.80])
  hold on
    % Plot Rawdata
    plot(x,y,'ko-')
      xIs = min(x):max(x);
      yIs = polyval(PolyIsat, xIs);
    % Plot Ion saturation current fit
    plot(xIs,yIs,'r-')
    % Plot limits for estimation of ion saturation current
    line(voltlim(1)*[1 1], [min(y) max(y)])
    line(voltlim(2)*[1 1], [min(y) max(y)])
  hold off
  mkplotnice('probe voltage (V)', 'probe current (A)', 12, '-25', '-50');
end
% XXX  <<<--- Check ion saturation current


% 2. Calculate the electron current
% This means removing the ion current from the measured characteristic.
% Below the plasma potential the ion current can be estimated by 
% calculating a linear regression for the ion saturation current.
%--------------------------------------------------------------------------
vol_avg = 6;
% = average volts * (points per Volt)
par_ysm = ceil( vol_avg * (numel(x) / (x(end)-x(1))) );

ysm = smooth(y, par_ysm);
Iisat = polyval(PolyIsat, x);
ye = ysm - Iisat;


% 3. Find the plasma potential
% Calculate the derivative (of the electron current) and search its 
% minimum, i.e. the plasma potential can be estimated by the position 
% of the maximum slope of the characteristic. Or the position, where the
% second derivative crosses zero min(I') --> I" = 0
%--------------------------------------------------------------------------
Dye = diff_discrete(x, ye);
% 3.1 smooth the derivative
Dyesm = smooth(Dye, par_ysm/2);

% Find plasma potential: minimum of first derivative (upper boundary)
%---------------------------------------------------------------------
[~, iVp] = min(Dyesm);
if isempty(fixpar.Vp)
  iudata.Vp = x(iVp);
else
  iudata.Vp = fixpar.Vp;
  ind = find(x<iudata.Vp);
  iVp = ind(end);
end
% Find zero crossing - floating potential (lower boundary)
%----------------------------------------------------------

% A: Take the whole characteristic until the plasma potential
posleft = 1;
% % B: Cut from floating potential
% ii = find(ysm<0);
% xvf = ii(1);
% posleft = xvf;

% Cut fit x vector
xcut  =  x(posleft:iVp);
yecut = ye(posleft:iVp);

% 4. Calculate the fit to the characteristic
% The fit parameters are the n_e and T_e.
%--------------------------------------------------------------------------
natconst
Ap = 2*pi*data.r*data.l;
%--------------------------
% 4.1 Define fit-functions
%--------------------------
% Between the floating potential and the plasma potential an exponential 
% increase of the electron current is assumed:
% fit formula: I(V) = I_e,sat * exp(-e (phi_p - V)/T_e)
%              I_e,sat = -e_0 * n_e * A_probe * <v+>
%              -->     = -e_0 * n_e * A_probe * sqrt(e_0 T_e/(2*pi*m_e))
% a: n_e
% b: T_e
% x: V probe potential
% y: I_e electron current
%     y = -e0*A_p*sqrt(e0/(2*pi*m_e))*ne*sqrt(T_e)*exp(-(phi_p - x)/T_e);
% --> y = s1*a*sqrt(b)*exp(-(s2-x)/b);
s1 = sprintf( '%5.5e', -e0*Ap*sqrt(e0/(2*pi*m_e)) );
s2 = sprintf( '%5.5e', iudata.Vp);
fitfct = fittype([s1 '*a*sqrt(b)*exp(-(' s2 ' - x)/b)']);
opts = fitoptions(fitfct);
%--------------------------------
% 4.2 Determine start parameters
%--------------------------------
% Start parameter T_e: T_e = -I_e(V_p) / I_e'(V_p)
Te_start = ye(iVp) / Dyesm(iVp);
% Start parameter n_e: n_e = -I_e(V_p)/[Ap*e0*sqrt{e0*Te_start/(2*pi*m_e)}]
ne_start = -ye(iVp)/(Ap*e0*sqrt(e0*Te_start/(2*pi*m_e)));
set(opts,'start',[ne_start    Te_start]);
opts.Upper =     [1e23        100];
opts.Lower =     [1e10          0];
%--------------------------------
% 4.3 Fit
%--------------------------------
ffit = fit(xcut, yecut, fitfct, opts);
Iefit = ffit(x);

iudata.ne = ffit.a;
iudata.Te = ffit.b;

%--------------------------------------
% 4.4 Calculate the floating potential
%--------------------------------------
% With smoothing, the potential, where the current vanishes is looked for.
% Sum fit of electron current and linear fit of ion saturation current:
yfitei = Iefit + Iisat;
  il = find(yfitei> 0);
  ir = find(yfitei<=0);
  % 2nd order polynomial at zero crossing
  ind = il(end-2):ir(2);                  % incices at zero crossing
  pf = polyfit(x(ind), yfitei(ind), 2);   % polynomial
  % Find y=zero crossings of parabula:
  hp = pf(2)/pf(1); hq = pf(3)/pf(1);
  x01 = -hp/2 + sqrt((hp^2)/4 - hq);
  x02 = -hp/2 - sqrt((hp^2)/4 - hq);
  % Find closest to I_e+I_i fit function
  if abs(x01-x(il(end))) < abs(x02-x(il(end)))
    % x01 is the zero crossing
    iudata.Vf = x01;
  else
    % x02 is the zero crossing
    iudata.Vf = x02;
  end


% XXX  <<<--- Check fit of electron current between Vf and Vp
figeps(12,8,2,30,90);
axes('position', [0.15 0.15 0.80 0.80])
hold on
plot(data.voltage, ye, 'ko')
plot(x,Iefit, 'r')
set(gca, 'xlim', [data.voltage(1) data.voltage(end)], ...
  'ylim', [min(data.current) max(data.current)])
hold off
mkplotnice('probe voltage (V)', 'fit I_e (A)', 12, '-25', '-50');
% XXX  <<<--- Check ion saturation current


% 4.5 *** Calculate the deviation of the fit from the rawdata
%---------------------------------------------------------
% ***Calculate the Normalized Integral Deviation
% NID = int_V |I_mod(V) - I_raw(V)| dV / int_V I_raw(V) dV

end