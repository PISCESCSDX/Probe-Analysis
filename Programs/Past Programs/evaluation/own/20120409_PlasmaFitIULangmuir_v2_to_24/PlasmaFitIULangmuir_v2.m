function [iudata, raw] = PlasmaFitIULangmuir_v2(data, checkonoff, fixpar, Vint, ...
  vol_avg)
%==========================================================================
%function [iudata] = PlasmaFitIULangmuir_v2(data, checkonoff, fixpar, ...
% Vint, vol_avg)
% Last Change: 2012-04-03 C. Brandt, San Diego
%--------------------------------------------------------------------------
% PlasmaFitIULangmuir_v2 calculates the plasma parameters 
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


% 1. Cut Characteristic (both sides 5%)
%--------------------------------------------------------------------------
% Extract the characteristic and the probe parameters
% Cut margin: 5% both sides
L5 = round(0.05*length(data.voltage));
ind = L5:length(data.voltage)-L5;
data.voltage = data.voltage(ind);
data.current = data.current(ind);

x = data.voltage;
y = data.current;


% 2. Exclusion Test: I-V characteristic is asymmetric
%--------------------------------------------------------------------------
% First simple test whether data may be Langmuir characteristic
% (average of left side must be larger than average of right side)
if mean(y(1:round(length(x)/2))) <= mean(y(round(length(x)/2):end))
  disp(['I-V characteristic is not asymmetric enough, i.e. most ' ...
      'probable not evaluable.'])
  iudata = []; raw = [];
  return
end


% 3. Calculate the ion saturation current
%--------------------------------------------------------------------------
%data.voltlim = []; % *** Is this line needed?
[PolyIsat, voltlim, ~, raw.iVf] = subPlasmaIsatPoly_v2(data);


% 4. Exclusion Test: Ion saturation current was found
%--------------------------------------------------------------------------
if isempty(PolyIsat)
  disp('No results can be obtained from the IV-characteristic.')
  iudata = []; raw = [];
  return
end
% XXX --->>> Check ion saturation current detection:
if checkonoff
  figeps(12,8,1,0,100); clf;
  axes('position', [0.20 0.16 0.78 0.78])
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


% 5. Calculate the electron current (by subtraction of ion sat. current)
% This means removing the ion current from the measured characteristic.
% Below the plasma potential the ion current can be estimated by 
% calculating a linear regression for the ion saturation current.
%--------------------------------------------------------------------------
% Smoothing points: = (AVG over Volt) * (points per Volt)
par_ysm = ceil(vol_avg * (numel(x)/(x(end)-x(1))) );
% Filter and smooth characteristic
ysm = filtmooth(y, par_ysm);
% Calculate the ion saturation current (linear polynom)
Iisat = polyval(PolyIsat, x);
% Subtract the ion saturation current to obtain basically the e-current
ye = ysm - Iisat;
% Output: rawdata (variable raw) of logarithmic electron current
yelog = real(log(ye));
raw.V = x;
raw.Ielog = yelog;


% 6. Exclusion Test: Maximum of first derivative of electron current
%    (approx plasma potential) must not be at the right margin
%    (because then the plasma potential as well as the electron saturation
%    current is most probably beyond the measurement)
%--------------------------------------------------------------------------
%  Smooth logarithm of Ie over 18 Volts (empirically found to work)
L_sm  = ceil(18 * (numel(x)/(x(end)-x(1))) );
% First derivative of electron current
Dyelog  = diff_discrete(x, yelog);
% Strong Smoothing & Filtering of D(yelog)
Dyelog_cut_sm = filtspline(x, Dyelog, L_sm, 0.5);
% Exclusion if maximum of 1st derivative of yelog is at the right margin
[~, imax] = max(Dyelog_cut_sm);
if imax==length(yelog)
  disp('The maximum of smoothed DI is at the right margin! No evaluation!')
  iudata = []; raw = [];
  return
end


% 7. Cut out Ion Saturation Regime
% The absolute value of the 2nd derivative of the logarithm of the
% electron current [i.e. abs(Diff2(yelog))] shows strong peaks in the ion
% saturation regime, since it corresponds to noise.
%--------------------------------------------------------------------------
% Calculate 2nd derivative of logarithm of electron current
D2yelog = diff_discrete(x, Dyelog);
% Take the absolute value and smooth
AbsD2yelog = abs(D2yelog);
AbsD2yelog_sm = smooth(AbsD2yelog,1);
% Calculate total average value
% (everything larger belongs to noise fluct. in the ion saturation regime)
avg = mean(AbsD2yelog_sm);
ind = AbsD2yelog_sm<avg;
AbsD2yelog_sm(ind) = 0;
[~,imax,~,~] = extrema(AbsD2yelog_sm);
% Cut out exponential and saturation region of electron current
i_Ieexp_start = max(imax);

% Average and Maximum Value of ion saturation region
yelog_noise_avg = mean(yelog(1:i_Ieexp_start));
yelog_noise_max = max(yelog(1:i_Ieexp_start));

% Cut everything below the noise level: Find first index from i_Ieexp_start
% which is larger than the noise level
ind = yelog > yelog_noise_max;
inum = 1:length(yelog); inum = inum(ind);
% Redefine the exponential start index
i_Ieexp_start = inum(1);
% New cut electron current  regime
ind_Ie = i_Ieexp_start:length(yelog);
xelog_cut =     x(ind_Ie);
yelog_cut = yelog(ind_Ie);


% 8. Extract exponential electron current region
%--------------------------------------------------------------------------
% Local maxima of D(Ie) correspond to lowest Te estimation.
% Smooth 1st derivative of Ie-region  and  look at local maxima:
% (Smooth over 2 Volts)
pts = ceil(2 * (numel(x)/(x(end)-x(1))) );
Dyelog_sm = filtspline(x, Dyelog, pts, 0.5);
Dyelog_cut_sm = filtspline(xelog_cut, Dyelog(ind_Ie), pts, 0.5);
[~,imax,~,~] = extrema(Dyelog_cut_sm);
% Define: Max of D(Ie) must be in the first half of the Ie-region
while imax(1) > length(Dyelog_cut_sm)/2
  imax = imax(2:end);
end

% 8.1 Look at points within a certain level-range around the local maximum:
% # INPUT: Level-Range of Slope (parameter with strong influence on T_e!!)
DIeexp_level = 0.2;
i_start = imax(1) + i_Ieexp_start - 1;
value = Dyelog_sm(i_start);
% Go to the left
status = 0; ctr=i_start;
while status==0
  ctr=ctr-1;
  if abs(value-Dyelog_sm(ctr)) > DIeexp_level*value
    ctr_le = ctr+1;
    status = 1;
  end
end
% Go to the right
status = 0; ctr=i_start;
while status==0
  ctr=ctr+1;
  if abs(value-Dyelog_sm(ctr)) > DIeexp_level*value
    ctr_ri = ctr-1;
    status = 1;
  end
end
% Define indices of exponential electron current region
ind_Ieexp = ctr_le:ctr_ri;
% Index where maximum of D(Ie) is detected (= minimum of Te)
ind_Ieexp_max = imax(1) + i_Ieexp_start - 1;
%----------------- !!! imax(1) is also used below !!! ---------------------
% Find the Start Index of the electron saturation region
% (local minimum of 2nd Derivative of Ie):
D2yelog = diff_discrete(x, Dyelog);
D2yelog_sm = diff_discrete(xelog_cut, Dyelog_cut_sm);
[~,~,~,imin] = extrema(D2yelog);
  ind = sort(imin);
  iminnew = ind>imax(1);
  iminnew = ind(iminnew);
  j = iminnew(1);
i_Iesat_start = j + i_Ieexp_start - 1;

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXX --->>> Check Plots: (1) log Ie, (2) D(log(Ie)), (3) D2(log(Ie))
% if checkonoff
figeps(12,15,3,30,100); clf; 
subplot(3,1,1)
hold on
plot(x, yelog, '-o')
plot(xelog_cut, yelog_cut, 'r-o')
hold off

subplot(3,1,2)
hold on
plot(x, Dyelog, '-o')
plot(xelog_cut, Dyelog_cut_sm, 'r-o')
plot(xelog_cut(ind_Ieexp_max-i_Ieexp_start+1), ...
  Dyelog_cut_sm(ind_Ieexp_max-i_Ieexp_start+1), 'rx')
hold off

subplot(3,1,3)
hold on;
plot(xelog_cut, D2yelog_sm, 'r-o')
plot(x, D2yelog, 'b-o')
plot(x(i_Iesat_start), D2yelog(i_Iesat_start), 'rx')
hold off

figeps(12,11,4,30,0);clf;
hold on
plot(x,yelog, 'bo-')
plot(x(ind_Ie),yelog(ind_Ie), 'g*-')
plot(x(ind_Ieexp),yelog(ind_Ieexp), 'r*')
line([x(1) x(end)], yelog_noise_avg*[1 1], 'Color', 'k')
line([x(1) x(end)], yelog_noise_max*[1 1], 'Color', 'b')
hold off

% end
% XXX --->>> Check ion saturation current detection:
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% Calculate the electron temperature from estimates of the slope
raw.Te2 = 1 /       Dyelog_cut_sm(imax(1));
raw.Te3 = 1 / mean( Dyelog(ind_Ieexp) );
disp(num2str([raw.Te2 raw.Te3]))



% 8. Find the plasma potential
% Start at maximum of 1st Derivative and go to the left, find zero crossing
i_start = ind_Ieexp_max;
value = D2yelog(i_start);
% Go to the left
status = 0; ctr=i_start;
while status==0
  ctr=ctr-1;
  if sign(D2yelog(ctr))>sign(D2yelog(ctr+1))
    ctr_le = ctr;
    status = 1;
  end
end
disp( num2str( [ctr_le x(ctr_le)]) )


[~,imax,~,imin] = extrema( filtspline(x,D2yelog,2,0.9) );
le = imax<i_start; le = imax(le); le = sort(le); le=le(end);
ri = imin>i_start; ri = imin(ri); ri = sort(ri); ri=ri(1);
D2yelog_sm = filtspline(x,D2yelog,2,0.9);
% Find Zero Crossing (Inflection Point) = Plasma Potential
% Go from ri to left and find zero crossing
ctr = ri;
status=0;
while status==0
  ctr=ctr-1;
  if sign(D2yelog_sm(ctr))>sign(D2yelog_sm(ctr+1))
    ctr_le2 = ctr;
    status = 1;
  end
end
disp( num2str( [ctr_le2 x(ctr_le2)] ) )

figeps(12,14,1,0,100);clf;
subplot(3,1,1); 
hold on
plot(  yelog,'k-o')
ind = i_start-10:length(yelog);
line([ctr_le  ctr_le] , [min(yelog(ind)) max(yelog(ind))], 'Color', 'b')
line([ctr_le2 ctr_le2], [min(yelog(ind)) max(yelog(ind))], 'Color', 'r')
hold off
set(gca, 'xlim', [ind(1) ind(end)])

subplot(3,1,2); 
plot(Dyelog_sm,'b-o')
set(gca, 'xlim', [i_start-10 length(yelog)])

subplot(3,1,3);
hold on
plot(D2yelog,'k-.')
plot(D2yelog_sm, 'b-o')

plot(le, D2yelog_sm(le) , 'r-*')
plot(ri, D2yelog_sm(ri) , 'c-*')

hold off
set(gca, 'xlim', [i_start-10 length(yelog)])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


iVp = ctr_le2;



%if checkonoff
  input('>>> Press Enter To Proceed <<<')
%end




Dye = diff_discrete(x, ye);
% 3.1 smooth the derivative (empiric amount of smoothing points)
Dyesm = filtmooth(Dye, 3*par_ysm);

if isempty(fixpar.Vp)
  Vp_start = x(iVp);
else
  Vp_start = fixpar.Vp;
  ind = find(x<Vp_start);
  iVp = ind(end);
end


% If no Plasma Potential could be detected, finish this m-file.
if isempty(iVp)
  disp('No plasma potential could be found! Return no Output!')
  iudata = [];
  return
end
 
ii=find(ye>0.01*ye(iVp));
iVf = ii(end);
    % If no Floating Potential could be detected, finish this m-file.
    if isempty(ii)
      disp('No floating potential could be found! Return no Output!')
      iudata = [];
      return
    end
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
%ye = ye-max(ye); %*** added 04.04.2012
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
  ne_start = -ye(iVp) / (Ap*e0*sqrt(e0*Te_start/(2*pi*m_e)));
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
    Vp1 = Vp_start + 0*fixpar.Vplim(2); % *** 4.4.12
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
    [ffit, gof] = fit(xcut, yecut, fitfct, opts);
      % % Chi Square Deviation: gof.rsquare is exactly:
      % sumerrfit = sum( (yecut-Iefit(indcut)).^2 );
      % sumdiffmn = sum( (yecut - mean(yecut)).^2 );
      % iudata.Rsqr = 1 - sumerrfit/sumdiffmn;
    iudata.Rsqr = gof.rsquare;
    
% ****     % Calculate fit function along original x-vector (here e-current)
     Iefit = ffit(x);
%     Iefit2= ffit(xcut);
%     figure; hold on
%     plot(xcut, yecut, 'o')
%     ne = ffit.a; Te = ffit.b; Vp = ffit.c; d = ffit.d;
%     plot(xcut, -e0*Ap*sqrt(e0/(2*pi*m_e)) *ne*sqrt(Te).*exp(-(Vp - xcut)./Te)+d, 'g');
%     
% y1 = filtmooth(real(log(ye)),30);
% y2 = csaps(x, y1, 0.5, x);
% y3 = smooth(diff_discrete(x,y2),20);
% yfit = csaps(x, y3, 0.5, x);
% figure
% hold on
% plot(x, y1, 'o')
% plot(x, y3, 'r')

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
    iudata.Vf = x01;
  else
    % x02 is the zero crossing
    iudata.Vf = x02;
  end
  
else
  iudata.Vf = NaN;
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
puttextonplot(gca, [0 0], 10, 10, ['V_f=' sprintf('%4.2f', iudata.Vf) 'V'], 0, 12, 'k');
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
    iudata.Vf = NaN;
    iudata.Rsqr = NaN;
    iudata.NID = NaN;
    return    
  end
else  % Found at least one imaginary number, i.e. a fit makes no sense.
  disp('No fit: Boundary values of ne, Te and Vp are not real numbers.')
    iudata.ne = NaN; iudata.neerr = NaN;
    iudata.Te = NaN; iudata.Teerr = NaN;
    iudata.Vp = NaN; iudata.Vperr = NaN;
  iudata.Vf = NaN;
  iudata.Rsqr = NaN;
  iudata.NID = NaN;
  return
end

end