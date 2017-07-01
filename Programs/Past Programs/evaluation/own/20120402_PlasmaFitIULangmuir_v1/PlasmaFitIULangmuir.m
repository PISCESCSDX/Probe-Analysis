function [iudata] = PlasmaFitIULangmuir(data, checkonoff, fixpar, Vint, ...
  vol_avg)
%==========================================================================
%function [iudata] = PlasmaFitIULangmuir(data, checkonoff, fixpar, ...
% Vint, vol_avg)
%--------------------------------------------------------------------------
% PlasmaFitIULangmuir calculates the plasma parameters 
% n_e (electron density), T_e (electron temperature), 
% V_p (plasma potential) and the V_f (floating potential) from a 
% current-voltage characteristic using a cylindrical electrode 
% (Langmuir probe).
% === External m-files ===
% EXTREMUM, FIGEPS, PUTTEXTONPLOT, MKPLOTNICE
% === Problems ===
% If ion saturation current is found to be increasing no floating
% potential may be found.
%
% If the plasma potential is not found look to procedure
% "Check derivative for minima positions"
%
% Sub-m-files: subPlasmaIsatPoly
% Todo: Optimize Finding Procedure of the plasma potential
% Todo: Include error calculation for T_e, n_e, V_p and V_f
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
%     fixpar.Te: input electron temperature
%     fixpar.Telim: lower, upper T_e limits [Te1 Te2], e.g. [3 4]
%     fixpar.ne: input electron density AND/OR
%     fixpar.nelim: lower, upper n_e limit [ne1 ne2], e.g. [1e15 1e17]
%     fixpar.Vp: input plasma potential AND/OR
%     fixpar.Vplim: (!)Delta low, upp V_p limits [DVp1 DVp2], e.g. [-2 1]
%                   if empty: no variation of Vp, just the extremum of I'
%     Vint: voltage interval: as below or mixed possible
%           [] OR [NaN, NaN]: [V_f, V_p]
%           [NaN, 20]: [Vf, 20]
%           [-40, NaN]: [-40, Vp]
%           [-40, 20]: [-40, 20]
%           [-Inf, 20]: [left end, 20]
%           [-40, Inf]: [-40, right end]
%           [-Inf, Inf]: [-40, right end]
%     vol_avg: smooth characteristic over voltage range (default 3V)
%              (used for smoothing the electron current)
%OUT: iudata.ne: electron density (m^-3)
%     iudata.Te: electron temperature (eV)
%     iudata.Vp: plasma potential (V)
%     iudata.Vf: floating potential (V)
%--------------------------------------------------------------------------
% EX:
% checkonoff=1; fixpar.Te=[];fixpar.Telim=[];fixpar.ne=[];...
% fixpar.nelim=[];fixpar.Vp=[];fixpar.Vplim=[-2 2];Vint=[-Inf, NaN];...
% vol_avg=6;
% iudata=PlasmaFitIULangmuir(data,checkonoff,fixpar,Vint,vol_avg);
%==========================================================================

% Set default voltage smoothing range
if nargin <5
  vol_avg = 3;
end

% Set default voltage cut interval
if nargin <4
  Vint = [-Inf, NaN];
end

% Set default fixed parameters
if nargin <3
  fixpar.Te = []; fixpar.Telim = [];
  fixpar.ne = []; fixpar.nelim = [];
  fixpar.Vp = []; fixpar.Vplim = [-1 1];
end



% Extract the characteristic and the probe parameters
% Cut margin: 5% both sides
L5 = round(0.05*length(data.voltage));
ind = L5:length(data.voltage)-L5;
data.voltage = data.voltage(ind);
data.current = data.current(ind);

x = data.voltage;
y = data.current;

% Some Tests: I-V data has the shape of a typical characteristic
%--------------------------------------------------------------------------
% First simple test whether data may be Langmuir characteristic
% (average of left side must be larger than average of right side)
if mean(y(1:round(length(x)/2))) <= mean(y(round(length(x)/2):end))
  disp(['I-V characteristic is not asymmetric enough, i.e. most ' ...
      'probable not evaluable.'])
  iudata = [];
  return
end


% 1. Estimate the ion saturation current (use filter)
%--------------------------------------------------------------------------
data.voltlim = [];
[PolyIsat, voltlim, ~] = subPlasmaIsatPoly(data);

if isempty(PolyIsat)
  disp('No results can be obtained from the IV-characteristic.')
  iudata = [];
  return
end


% XXX --->>> Check ion saturation current detection:
if checkonoff
  figeps(12,8,1,0,100); clf;
  axes('position', [0.20 0.20 0.78 0.78])
  title('detection of I_{sat} voltage-limits');
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
  set(gca, 'xlim', [data.voltage(1) data.voltage(end)], ...
  'ylim', [min(data.current) max(data.current)])
  mkplotnice('probe voltage (V)', 'probe current (A)', 12, '-25', '-50');
end
% XXX  <<<--- Check ion saturation current


% 2. Calculate the electron current
% This means removing the ion current from the measured characteristic.
% Below the plasma potential the ion current can be estimated by 
% calculating a linear regression for the ion saturation current.
%--------------------------------------------------------------------------
% = average volts * (points per Volt)
par_ysm = ceil( vol_avg * (numel(x) / (x(end)-x(1))) );

ysm = filtmooth(y, par_ysm);
Iisat = polyval(PolyIsat, x);
ye = ysm - Iisat;


% 3. Find the plasma potential
% Calculate the derivative (of the electron current) and search its 
% minimum, i.e. the plasma potential can be estimated by the position 
% of the maximum slope of the characteristic. Or the position, where the
% second derivative crosses zero min(I') --> I" = 0
%--------------------------------------------------------------------------
Dye = diff_discrete(x, ye);
% 3.1 smooth the derivative (empiric amount of smoothing points)
Dyesm = filtmooth(Dye, 3*par_ysm);

% Procedure: Finding the plasma potential: minimum of first derivative
%----------------------------------------------------------------------
% (i) Find all minima of I'
% (ii) Look which belongs to the plasma potential
%   a) Smooth I' strongly (ca. 30% smoothing length)
%   b) Find minima smaller than the "very smooth" function
%   c) Find minima most close to the ion saturation current regime
%------------------------------------------------------------------------
% (i)
[~,~,pmin,imin] = extrema(Dyesm);
  % Remove first of 60% and last element
  ii = imin < round(0.60*length(Dyesm));
  jj = imin==length(Dyesm);
% (ii.a) Smooth over 20% of the characteristic length
Dye_verysmooth = filtmooth(Dyesm, round(0.20*numel(Dyesm)) );
  % Remove minima positions ii and jj from the following criterion ii.b
  pmin(ii) = Dye_verysmooth(imin(ii));
  pmin(jj) = Dye_verysmooth(imin(jj));
% (ii.b)
  % Define the level below the "very smoothed" I' (= Dyesm)
  % Factor 2.2 found by testing. It seems to work reliable at small signal 
  % to noise data (but there may still be data that does not work.)
  % Define 'diff_min_smooth_avg': level below 'very smooth' I'
  diff_min_smooth_avg = 2.2*mean( pmin - Dye_verysmooth(imin) );
  % Find positions of minima below this level:
  i_absmin = pmin - Dye_verysmooth(imin) < diff_min_smooth_avg;
% (ii.c)
[iVp, ~] = min( imin(i_absmin) );

if isempty(fixpar.Vp)
  Vp_start = x(iVp);
else
  Vp_start = fixpar.Vp;
  ind = find(x<Vp_start);
  iVp = ind(end);
end



% XXX  >>>>--- Check derivative for minima positions
if checkonoff
figeps(12,13,2,30,100); clf;
axes('position', [0.2 0.60 0.75 0.35])
title('detection of V_p');
hold on
  plot(x,ysm)
  if ~isempty(iVp)
    line(x(iVp)*[1 1], [min(ysm) max(ysm)])
  end
hold off
set(gca, 'xlim', [data.voltage(1) data.voltage(end)], ...
'ylim', [min(data.current) max(data.current)])
mkplotnice('probe voltage (V)', 'I (A)', 12, '-30', '-30');
%
% Green: Very Smooth Derivative
% Magenta: Level (below "very smooth DIe") above local minima are neglected
axes('position', [0.2 0.10 0.75 0.35])
title('detection of V_p');
hold on
  plot(Dyesm,'r')
  plot(imin, Dyesm(imin), 'o')
  plot(Dye_verysmooth, 'g')
  plot(Dye_verysmooth + diff_min_smooth_avg*ones(length(Dyesm),1), 'm')
  if ~isempty(iVp)
    line(iVp*[1 1], [min(Dyesm) max(Dyesm)])
  end
hold off
set(gca, 'xlim', [1 length(Dyesm)], ...
'ylim', [min(Dyesm) max(Dyesm)])
mkplotnice('# samples', 'dI/dV (A/V)', 12, '-30', '-30');
end
% XXX  <<<--- Check derivative for minima positions


% If no Plasma Potential could be detected, finish this m-file.
if isempty(iVp)
  disp('No plasma potential could be found! Return no Output!')
  iudata = [];
  return
end
 
  
% Find zero crossing - floating potential (lower boundary)
%----------------------------------------------------------
% Cut fit x vector according to Vint
ii = find(ysm<0);
    % If no Floating Potential could be detected, finish this m-file.
    if isempty(ii)
      disp('No floating potential could be found! Return no Output!')
      iudata = [];
      return
    end
iVf = ii(1);
if isempty(Vint) || ( isnan(Vint(1)) && isnan(Vint(2)) )
  posleft = iVf;     posright = iVp;
else
  if isnan(Vint(1)); posleft  = iVf;      end
  if Vint(1)==-Inf;  posleft  = 1;        end
  if isnan(Vint(2)); posright = iVp;      end
  if Vint(2)==Inf;   posright = numel(x); end
  if Vint(1)>-Inf && Vint(1)<Inf
    posleft  = findind(x,Vint(1));
  end
  if Vint(2)>-Inf && Vint(2)<Inf
    posright = findind(x,Vint(2));
  end
end
indcut = posleft:posright;
xcut  =  x(indcut);
yecut = ye(indcut);
  

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
% c: phi_p, plasma potential
% x: V probe potential
% y: I_e electron current
%     y = -e0*A_p*sqrt(e0/(2*pi*m_e))*ne*sqrt(T_e)*exp(-(phi_p - x)/T_e);
% --> y = s1*a*sqrt(b)*exp(-(c-x)/b);
s1 = sprintf( '%5.5e', -e0*Ap*sqrt(e0/(2*pi*m_e)) );
fitfct = fittype([s1 '*a*sqrt(b)*exp(-(c - x)/b)']);
opts = fitoptions(fitfct);
%--------------------------------
% 4.2 Determine start parameters
%--------------------------------

% Start parameter T_e: T_e = -I_e(V_p) / I_e'(V_p)
if isempty(fixpar.Te)
  Te_start = ye(iVp) / Dyesm(iVp);
else
  Te_start = fixpar.Te;
end

% Start parameter n_e: n_e = -I_e(V_p)/[Ap*e0*sqrt{e0*Te_start/(2*pi*m_e)}]
if isempty(fixpar.ne)
  ne_start = -ye(iVp)/(Ap*e0*sqrt(e0*Te_start/(2*pi*m_e)));
else
  ne_start = fixpar.ne;
end
set(opts,'start',[ne_start  Te_start  Vp_start]);

% Set fit limits [ne, Te, Vp]
  if isempty(fixpar.nelim)
    n0 = 0.1*ne_start;      n1 = 10*ne_start;
  else
    n0 = fixpar.nelim(1);   n1 = fixpar.nelim(2);
  end

  if isempty(fixpar.Telim)
    T0 = 0.1*Te_start;      T1 = 10*Te_start;
  else
    T0 = fixpar.Telim(1);   T1 = fixpar.Telim(2);
  end

  if isempty(fixpar.Vplim)
    % Vp will be not varied
    Vp0 = Vp_start;         Vp1 = Vp_start;
  else
    Vp0 = Vp_start + fixpar.Vplim(1);
    Vp1 = Vp_start + fixpar.Vplim(2);
  end

%--------------------------------
% 4.3 Boundaries & Fit
%--------------------------------
% In case upper boundaries are smaller than lower boundaries or numbers are
% not real, do not fit, and output NaN.
if isreal([n0 n1 T0 T1 Vp0 Vp1])
  if n1>=n0 && T1>=T0 && Vp1>=Vp0
    opts.Upper = [n1 T1 Vp1];
    opts.Lower = [n0 T0 Vp0];
    
    opts.weights = ones(1,iVf-1);
    opts.weights(iVf:iVp) = 100;
    [ffit, gof] = fit(xcut, yecut, fitfct, opts);
      % % Chi Square Deviation: gof.rsquare is exactly:
      % sumerrfit = sum( (yecut-Iefit(indcut)).^2 );
      % sumdiffmn = sum( (yecut - mean(yecut)).^2 );
      % iudata.Rsqr = 1 - sumerrfit/sumdiffmn;
    iudata.Rsqr = gof.rsquare;
    
    % Calculate fit function along original x-vector (here e-current)
    Iefit = ffit(x);
      % Extract the fit parameters
      iudata.ne = ffit.a;
      iudata.Te = ffit.b;
      iudata.Vp = ffit.c;
      % Extract the confidence intervals
      c = confint(ffit);
      
      iudata.neerr = [c(1,1) c(2,1)];
      iudata.Teerr = [c(1,2) c(2,2)];
      iudata.Vperr = [c(1,3) c(2,3)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF CONTINUES BELOW AT "NO-OUTPUT" %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------
% 4.4 Calculate the floating potential
%--------------------------------------
% With smoothing, the potential, where the current vanishes is looked for.
% Before adding the linear ion saturation current fct. check for problems:
% Problem Case 1: Iisat has negative values -> make them to NaN
  ind = Iisat<0;
  Iisat(ind) = NaN;
% Sum fit of electron current and linear fit of ion saturation current:
yfitei = Iefit + Iisat;
  il = find(yfitei> 0);
  ir = find(yfitei<=0);
% Test whether points around zero crossing exist (may be NaN)
if ~isempty(il) && ~isempty(ir)
  % 2nd order polynomial at zero crossing
  ind = il(end-2):ir(2);                  % indices at zero crossing
  pf = polyfit(x(ind), yfitei(ind), 2);   % polynomial
  % Find y=zero crossings of parabula:
  hp = pf(2)/pf(1); hq = pf(3)/pf(1);
  x01 = -hp/2 + sqrt((hp^2)/4 - hq);
  x02 = -hp/2 - sqrt((hp^2)/4 - hq);
  % Find closest to I_e+I_i fit function
  if abs(x01-x(il(end))) < abs(x02-x(il(end)))
    % x01 is the zero crossing
    iudata.Vffit = x01;
  else
    % x02 is the zero crossing
    iudata.Vffit = x02;
  end
  
else
  iudata.Vffit = NaN;
end


% 5. Calculate the deviation of the fit from the rawdata
%--------------------------------------------------------
% Calculate the Normalized Integral Deviation from smooth Ie
integral1 = int_discrete(xcut,  abs(yecut-Iefit(indcut)) );
integral2 = int_discrete(xcut,  abs(yecut) );
%
% NID = int_V |I_mod(V) - I_raw(V)| dV / int_V I_raw(V) dV
iudata.NID = integral1/integral2;


% XXX  <<<--- Check fit of electron current between Vf and Vp
if checkonoff
figeps(12,8,3,60,100); clf;
axes('position', [0.20 0.20 0.78 0.78])
hold on
plot(data.voltage, ye, 'ko')
plot(x(indcut),Iefit(indcut), 'r')
set(gca, 'xlim', [data.voltage(1) data.voltage(end)], ...
  'ylim', [min(data.current) max(data.current)])
hold off
mkplotnice('probe voltage (V)', 'fit I_e (A)', 12, '-25', '-50');
puttextonplot(gca, [0 1], 10, -40, 'fit to I_e', 0, 14, 'k');
puttextonplot(gca, [0 0], 10, 90, ['R=' sprintf('%4.2f', iudata.Rsqr)], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 70, ['T_e=' sprintf('%4.2f', iudata.Te) 'eV'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 50, ['n_e=' sprintf('%4.2g', iudata.ne) 'm^{-3}'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 30, ['V_p=' sprintf('%4.2f', iudata.Vp) 'V'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 10, ['V_f=' sprintf('%4.2f', iudata.Vffit) 'V'], 0, 12, 'k');
end
% XXX  <<<--- Check ion saturation current

%%%%%%%%%%%%%%%
% "NO-OUTPUT" %
%%%%%%%%%%%%%%%

  else   % Found at least one lower boundary larger than a upper boundary.
   disp('No fit: Lower and upper boundaries are not logic.')
    iudata.ne = NaN; iudata.neerr = NaN;
    iudata.Te = NaN; iudata.Teerr = NaN;
    iudata.Vp = NaN; iudata.Vperr = NaN;
    iudata.Vffit = NaN;
    iudata.Rsqr = NaN;
    iudata.NID = NaN;
    return    
  end
else  % Found at least one imaginary number, i.e. a fit makes no sense.
  disp('No fit: Boundary values of ne, Te and Vp are not real numbers.')
    iudata.ne = NaN; iudata.neerr = NaN;
    iudata.Te = NaN; iudata.Teerr = NaN;
    iudata.Vp = NaN; iudata.Vperr = NaN;
  iudata.Vffit = NaN;
  iudata.Rsqr = NaN;
  iudata.NID = NaN;
  return
end

end