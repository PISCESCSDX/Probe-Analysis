function iudata = PlasmaFitIULangmuir_a11(data, p_chk, p_fix, p_Vint, ...
  p_Vavg, p_IexpLevel, p_ni)
%==========================================================================
%function iudata = PlasmaFitIULangmuir_a11(data, p_chk, p_fix, p_Vint, ...
% p_Vavg, p_IexpLevel, p_ni)
% Last Change: 2012-04-21 05:41 C. Brandt, San Diego
% This is a beta Version: a11. (last working version before: a10)
% Like version a10 this version provides the option 'p_ni' of calculating 
% ne from I_i,sat (p_ni='Iisat', ne=ni using Bohm criterion) or from 
% electron current fit (p_ni='Ie').
% But in this version another fit routine is used ('lsqcurvefit') since it
% works more reliable than the normal 'fit' method.
% (See details in the fit section below.)
%--------------------------------------------------------------------------
% PlasmaFitIULangmuir_a11 calculates the plasma parameters 
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
% IN: data.voltage, data.current, data.l, data.r, data.m_ion (in u0)
%     p_chk: 0: no plots, 1: check plots
%     p_fix.Te: input electron temperature
%     p_fix.Telim: lower, upper T_e limits [Te1 Te2], e.g. [3 4]
%     p_fix.ne: input electron density AND/OR
%     p_fix.nelim: lower, upper n_e limit [ne1 ne2], e.g. [1e15 1e17]
%     p_fix.Vp: input plasma potential AND/OR
%     p_fix.Vplim: (!)Delta low, upp V_p limits [DVp1 DVp2], e.g. [-2 1]
%                   if empty: no variation of Vp, just the extremum of I'
%     p_Vint: voltage interval: as below or mixed possible
%           [] OR [NaN, NaN]: [V_f, V_p]
%           [NaN, 20]: [Vf, 20]
%           [-40, NaN]: [-40, Vp]
%           [-40, 20]: [-40, 20]
%           [-Inf, 20]: [left end, 20]
%           [-40, Inf]: [-40, right end]
%           [-Inf, Inf]: [-40, right end]
%     p_Vavg: smooth characteristic over voltage range (default 3V)
%              (used for smoothing the electron current)
% p_IexpLevel: level range of slope (strong influence on T_e!!)
% p_ni: 'Iisat': calculate ne from ion saturation current
%          'Ie': calculate ne from fit to electron current
%OUT: iudata.ne: electron density (m^-3)
%     iudata.Te: electron temperature (eV)
%     iudata.Vp: plasma potential (V)
%     iudata.Vf: floating potential (V)
%--------------------------------------------------------------------------
% EX:
% p_chk=1; p_fix.Te=[];p_fix.Telim=[];p_fix.ne=[];...
% p_fix.nelim=[];p_fix.Vp=[];p_fix.Vplim=[-2 2];p_Vint=[-Inf, NaN];...
% p_Vavg=1;
% iudata=PlasmaFitIULangmuir_a11(data,p_chk,p_fix,p_Vint,p_Vavg);
%==========================================================================

% # INPUT: Level-Range of Slope (parameter with strong influence on T_e!!)
if nargin <6
  p_IexpLevel = 0.95;
end

% Set default voltage smoothing range
if nargin <5
  p_Vavg = 3;
end

% Set default voltage cut interval
if nargin <4
  p_Vint = [-Inf, NaN];
end

% Set default fixed parameters
if nargin <3
  p_fix.Te = []; p_fix.Telim = [];
  p_fix.ne = []; p_fix.nelim = [];
  p_fix.Vp = []; p_fix.Vplim = [-1 1];
end


%--------------------------------------------------------------------------
% Load natural physical constants
%--------------------------------------------------------------------------
natconst
% Probe area
A_p = 2*pi*data.r*data.l;
% Ion mass
m_i = data.m_ion*u0;


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
  ll = numel(y);
  dx = (x(end)-x(1)) / numel(x);
  indx = 1:ll;


% 2. Exclusion Test: I-V characteristic is asymmetric
%--------------------------------------------------------------------------
% First simple test whether data may be Langmuir characteristic
% (average of left side must be larger than average of right side)
if mean(y(1:round(length(x)/2))) <= mean(y(round(length(x)/2):end))
  disp(['I-V characteristic is not asymmetric enough, i.e. most ' ...
      'probable not evaluable.'])
  iudata = [];
  return
end


% 3. Calculate the ion saturation current
%--------------------------------------------------------------------------
[PolyIsat, ~, i_isat, Vf] = subPlasmaIsatPoly(data);


% 4. Exclusion Test: Ion saturation current was found & ne from Iisat
%--------------------------------------------------------------------------
if isempty(PolyIsat)
  disp('Could not detect I_i,sat.')
  iudata = [];
  return
else
  % Determine maximum Ion saturation current
  Iisat_max = polyval(PolyIsat, x(i_isat(1)) );
end
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXX --->>> Check ion saturation current detection:
if p_chk
  figeps(12,8,1,0,100); clf;
  axes('position', [0.20 0.16 0.78 0.72])
  title('I_{sat} voltage limits');
  hold on
    % Plot Rawdata
    plot(x,y*1e3,'ko-')
    ind_isat = i_isat(1):i_isat(2);
    plot(x(ind_isat),y(ind_isat)*1e3,'bo-')
    xIs = min(x):max(x);
    yIs = polyval(PolyIsat, xIs);
    % Plot Ion saturation current fit
    plot(xIs,yIs*1e3,'r-')
  hold off
  set(gca, 'xlim', [data.voltage(1) data.voltage(end)], ...
  'ylim', 1e3*[min(data.current) max(data.current)])
  mkplotnice('probe voltage (V)', 'probe current (mA)', 12, '-25', '-50');
end
% XXX  <<<--- Check ion saturation current
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


% 5. Calculate the electron current (by subtraction of ion sat. current)
% This means removing the ion current from the measured characteristic.
% Below the plasma potential the ion current can be estimated by 
% calculating a linear regression for the ion saturation current.
%--------------------------------------------------------------------------
% Smoothing points: = (AVG over Volt) * (points per Volt)
par_ysm = ceil(p_Vavg * (numel(x)/(x(end)-x(1))) );
% Filter and smooth characteristic
ysm = filtmooth(y, par_ysm);
% Calculate the ion saturation current (linear polynom)
Iisat = polyval(PolyIsat, x);
% Subtract the ion saturation current to obtain basically the e-current
ye = ysm - Iisat;
% Output: rawdata (variable raw) of logarithmic electron current
yelog = real(log(ye));


% 6. Exclusion Test: Maximum of first derivative of electron current
%    (approx plasma potential) must not be at the right or left margin
%    (because then the plasma potential as well as the electron saturation
%    current is most probably beyond the measurement)
%--------------------------------------------------------------------------
%  Smooth logarithm of Ie over 18 Volts (empirically found to work)
L_sm  = ceil(5/dx);
% First derivative of electron current
Dyelog  = diff_discrete(x, yelog);
% Strong Smoothing & Filtering of D(yelog)
Dyelog_sm = filtspline(x, Dyelog, L_sm, 0.5);
% Exclusion if maximum of 1st derivative of yelog is at the right margin
% (= region containing the exponential increase!)
[~, imax] = max(Dyelog_sm(i_isat(2):end));
if imax==length(yelog(i_isat(2):end)) || imax==1
 disp('The maximum of smoothed DI is at the left or right margin! No Ev!')
  iudata = [];
  return
end
 
%==========================================================================
% Most Difficult Automatic Procedures: Define limits of the fit range
%==========================================================================
% 7. Search for Characteristic Positions of the electron current
% 7.1 Start of exponential electron region (end of noise level)
%    - Start of electron saturation region
%    - Inflection point of characteristic (Plasma Potential)
% -------------------------------------------------------------------------
% 7.1 Rough estimation of the Start of exponential growth:
%     The absolute value of the 2nd derivative of the logarithm of the
%     electron current [i.e. abs(Diff2(yelog))] shows strong peaks in the 
%     ion saturation regime, since it corresponds to noise.
%     i_Ie_endofnoise: rough estimate
%     i_Ieexp_start: start exp. Ie (more to right hand side than i_Iee*1)
%     ind_Ie: indices of electron current (beginning at exp. increase)
%--------------------------------------------------------------------------
% Calculate 2nd derivative of logarithm of electron current
  D2yelog = diff_discrete(x, Dyelog);
% Take the absolute value and smooth
  AbsD2yelog = abs(D2yelog);
% Calculate total average value
% (everything larger belongs to noise fluct. in the ion saturation regime)
  avg = mean(AbsD2yelog);
  ind = AbsD2yelog<avg;
  AbsD2yelog(ind) = 0;
  %-----------------------------------------------------
  % Help procedure LONGEST ZERO: Find longest zero part:
  ind = AbsD2yelog<eps;
  % This help function peaks positively when a zero part starts
  % and negatively when it ends:
  zero_peaks = diff_discrete(indx',sign(ind-0.1));
  zero_peaks(1) = 0;
  zmin = min(zero_peaks);        zmax = max(zero_peaks);
  pkpos = zero_peaks > zmax/10;   pkpos = indx(pkpos);
  pkneg = zero_peaks < zmin/10;   pkneg = indx(pkneg);
  % Delete all negative peaks below the first positive peak:
  if pkneg(1)<=pkpos(1)
    pkneg = pkneg(2:end);
  end
  % Negative Peaks may be missing one peak at the end: if so add one at end
  if numel(pkneg)-numel(pkpos) == -1
    pkneg(end+1) = ll;  % ll = length(y); Defined at Start
  end
  % Calculate the difference between the start and end peaks:
    % Control Plot: figure; hold on; plot(pkneg, 'o'); plot(pkpos, 'r*')
  pk_diff = (pkneg-1) - pkpos;
  % Find length of maximum zero area
  [pkmax, ~] = max(pk_diff);
  % If more equal-sized maximum areas exist take the first one!
  ind = pk_diff==pkmax; helpx = 1:numel(pkpos); ind = min(helpx(ind));
  i_Ie_endofnoise = pkpos(ind);
  % --- Check Plots: Detection of the Noise level End Positon ----
  %   close all; figure(1); clf;
  %   hold on
  %   plot(AbsD2yelog);
  %   plot(i_Ie_endofnoise,AbsD2yelog(i_Ie_endofnoise), 'r*')
  %   hold off
  %   input('Press')
  % --
  % END Help procedure LONGEST ZERO
  %------------------------------------------------------------------------
  
if isempty(i_Ie_endofnoise)
 disp('Did not find the end of the noise! No Evaluation!')
  iudata = [];
  return
end


% Rough Estimation of electron current regime from exponential I_e
%
% Average and Maximum of the noisy ion saturation region
  yelog_noise_avg = mean(yelog(1:i_Ie_endofnoise));
  yelog_noise_max = max(yelog(1:i_Ie_endofnoise));
% Cut everything below the noise level
  ind = yelog > yelog_noise_max;
  inum = indx(ind);
% Redefine the exponential start index
i_Ieexp_start = inum(1);
% New cut electron current regime
ind_Ie = i_Ieexp_start:length(yelog);


% 7.2 Find exponential electron current region
%--------------------------------------------------------------------------
% Local maxima of D(Ie) correspond to lowest Te estimation.
% Smooth 1st derivative of Ie-region  and  look at local maxima:
% (Smooth over 2 Volts)
pts = ceil( 2/dx );
% Strong filtering:
Dyelog_sm = filtspline(x, Dyelog, pts, 0.01);
[~,imax,~,~] = extrema(Dyelog_sm);
% Define: Max of D(Ie) must be in the first half of the Ie-region (!)
% ii: Find Maximum of Slope
  ii = imax >= i_Ie_endofnoise-2*pts; ii = imax(ii);
  ii2= ii < (numel(yelog)-ceil(numel(ind_Ie)/2)); ii2 = ii(ii2);
if  numel(ii2)==0
 disp('Error: ii2 empty!')
  iudata = [];
  return
  
  else
    imax = ii2(1);
end
  

% Look at points within a certain level-range around the local maximum
% of the 1st derivative (find the exponential slope region)
% Left and Right Indices
ind = Dyelog_sm < p_IexpLevel*Dyelog_sm(imax);
  i_lower = indx(ind);
  ile = i_lower<imax; ile = i_lower(ile); ile = ile(end);
  iri = i_lower>imax; iri = i_lower(iri); iri = iri(1);
% Define indices of exponential electron current region
ind_Ieexp = ile:iri;
% Index where maximum of D(Ie) is detected (= minimum of Te)
i_dIe_max = imax;


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Check Plots: (1) rawdata log(I_e), (2) smoothed log(I_e)
if p_chk
figeps(10,15,1); clf;
subplot(2,1,1)
yelog_raw = real(log(-y+polyval(PolyIsat,x)));
hold on
plot(yelog_raw, 'bo-')
plot(ind_Ie,yelog_raw(ind_Ie), 'g*-')
plot(ind_Ieexp,yelog_raw(ind_Ieexp), 'r*')
hold off
set(gca,'xlim',[1 numel(yelog)])
mkplotnice('', '', 12, '0', '0');
title('log(rawdata ye)')

subplot(2,1,2)
hold on
plot(yelog, 'bo-')
plot(ind_Ie,yelog(ind_Ie), 'g*-')
plot(ind_Ieexp,yelog(ind_Ieexp), 'r*')
line([1 numel(yelog)], yelog_noise_avg*[1 1], 'Color', 'k')
line([1 numel(yelog)], yelog_noise_max*[1 1], 'Color', 'b')
hold off
set(gca,'xlim',[1 numel(yelog)])
mkplotnice('', '', 12, '0', '0');
title('log(smoothed ye)')
end
% XXX <<<---
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



% 7.2.1 Histogram Statistics of the Characteristic (Double Check 7.1)
%--------------------------------------------------------------------------
% If 7.1/7.2 fails finding the electron slope, this procedure may grip.
%
% It looks at the histogram of log(I_e). A usable characteristic reveals 
% two maxima, one at the ion saturation regime and the other at the 
% electron current regime.
% The minimum between both serves here for the determination of another
% index, which will be used instead of the end found for the
% exponential slope: ind_Ieexp (This procedure relies on the determination
% of the noise level - for low noise levels, the above procedure [7.1]
% yields no correct values).
%--------------------------------------------------------------------------
[hy,hx] = hist(yelog,round(ll/2));
  yout = filtspline(hx, hy, 2, 0.9);
  [~,~,~,histmin] = extrema(yout);
  histmin = sort(histmin);
  yemin = hx(histmin);
  % Convert from histogram transform to x Indices
  imin = zeros(length(yemin),1);
  for j=1:length(yemin)
    ind = yelog > yemin(j);
    imin(j) = min(indx(ind));
  end
  % Check if minima are present anyway
  if ~isempty(imin)
    % If the desired minimum is at the end:
    if histmin(end)==length(hx)
        % If more than one minimum is detected:
        if length(histmin)>1
          % Next minimum larger than ind_Ieexp end index?
          if imin(end-1) > ind_Ieexp(end)
            i_fit_ie_end = imin(end-1);
          else
            i_fit_ie_end = ind_Ieexp(end);
          end
        else
          i_fit_ie_end = ind_Ieexp(end);
        end
    else
      if imin(end) > ind_Ieexp(end)
        i_fit_ie_end = imin(end);
      else
        i_fit_ie_end = ind_Ieexp(end);
      end
    end
  else
    i_fit_ie_end = ind_Ieexp(end);
  end
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
if p_chk
figeps(22,10,2,0,0); clf;
subplot(2,2,2)
    hist(y,25);     title('histogram of lin. rawdata')
subplot(2,2,4)
    hold on
      plot(hx,hy);
      plot(hx,yout, 'r')
    hold off;       title('histogram log(I_e)')
subplot(2,2,3)
    plot(yelog)
    line(i_fit_ie_end*[1 1], [min(yelog) max(yelog)], 'color', 'r')
    title('log(I_e)')
subplot(2,2,1)
    plot(y)
    line(i_fit_ie_end*[1 1], [min(y) max(y)], 'color', 'r')
    title('I, rawdata')
end
% XXX <<<---
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


% 7.3 Inflection point of characteristic (Plasma Potential)
%     Start at maximum of 1st Derivative, go to left, find zero crossing
i_start = i_dIe_max;
[~,~,~,imin] = extrema( filtspline(x,D2yelog,2,0.9) );
ri = imin>i_start; ri = imin(ri); ri = sort(ri); ri=ri(1);
D2yelog_sm = filtspline(x,D2yelog,2,0.9);
% Find Zero Crossing (Inflection Point) = Plasma Potential
% Go from ri to left and find zero crossing
ctr = ri;
status=0;
while status==0
  ctr=ctr-1;
  if sign(D2yelog_sm(ctr))>sign(D2yelog_sm(ctr+1))
    i_Ie_inflection = ctr;
    status = 1;
  end
end



% 8. Define rough estimate of plasma potential
%--------------------------------------------------------------------------
% i_Ieexp_start: is the first index above the noise level of log(I_e)
% i_dIe_max:     is the index of the maximum slope of log(I_e)
% (also a possible option:  iVp = i_Ie_inflection)
% For a rough estimation of the Plasma Potential use the largest of both:
if i_Ieexp_start > i_dIe_max
  iVp = i_Ieexp_start;
else
  iVp = i_dIe_max;
end




% 9. Define rough estimate of plasma potential (for fit-start-parameters)
%--------------------------------------------------------------------------
Dye = diff_discrete(x, ye);
% 3.1 smooth the derivative (empiric amount of smoothing points)
Dyesm = filtmooth(Dye, 3*par_ysm);

[~,~,~,imin] = extrema(Dyesm);
imin2 = sort(imin);
dimin = round( (imin2(end)-imin2(1)) / numel(imin2) );
Dyesm2= filtspline(x, Dye, round(1*dimin), 0.1);  % Best Found by Trying
[~,~,~,imin] = extrema(Dyesm2);
while imin(1)==1 || imin(1)==ll
 imin = imin(2:end);
end
i_fit_ie_end_2 = imin(1);



%TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
% Figures showing the different important indices:
if p_chk
close all; figeps(20,15,1,0,100); clf;
subplot(2,2,1)
hold on
plot(ye,'bo-')
line(i_fit_ie_end*[1 1], [min(ye) max(ye)])
line(i_fit_ie_end_2*[1 1], [min(ye) max(ye)], 'Color', 'r')
hold off
dy=20;
puttextonplot(gca, [0 0], 5, 100-0*dy, ...
  ['i-Ie-endofnoise: ' num2str(i_Ie_endofnoise)], 0, 8, 'k');
puttextonplot(gca, [0 0], 5, 100-1*dy, ...
  ['i-Ieexp-start: ' num2str(i_Ieexp_start)], 0, 8, 'k');
puttextonplot(gca, [0 0], 5, 100-2*dy, ...
  ['i-Ie-inflection: ' num2str(i_Ie_inflection)], 0, 8, 'k');
puttextonplot(gca, [0 0], 5, 100-3*dy, ...
  ['i-fit-ie-end ' num2str(i_fit_ie_end)], 0, 8, 'b');
puttextonplot(gca, [0 0], 5, 100-4*dy, ...
  ['i-fit-ie-end-2 ' num2str(i_fit_ie_end_2)], 0, 8, 'r');
title('Comparison Fit Ends')

subplot(2,2,3)
hold on
plot(Dyesm,'bo-')
line(i_fit_ie_end*[1 1], [min(Dyesm) max(Dyesm)])
hold off

subplot(2,2,2)
hold on
plot(yelog,'ko-')
line(i_fit_ie_end*[1 1], [min(yelog) max(yelog)])
hold off

subplot(2,2,4)
hold on
plot(Dyelog,'ko-')
line(i_fit_ie_end*[1 1], [min(Dyelog_sm) max(Dyelog_sm)])
line(i_fit_ie_end*[1 1], [min(Dyelog_sm) max(Dyelog_sm)], 'Color', 'r')
hold off
end
%TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

% *** Delete ???
% if isempty(p_fix.Vp)
% %   Vp_start = x(iVp);
%   Vp_start = x(round((i_fit_ie_end_2+i_fit_ie_end)/2));
% else
%   Vp_start = p_fix.Vp;
%   ind = find(x<Vp_start);
%   iVp = ind(end);
% end


% 9.1 Excusion test: If no Plasma Potential detected, finish this m-file.
if isempty(iVp)
  disp('No plasma potential could be found! No Output!')
  iudata = [];
  return
end


% 10. If necessary correct floating potential to left (*** bad method)
%--------------------------------------------------------------------------
if Vf.i>iVp
  iVp = Vf.i +1;
end


%==========================================================================
% 12. Fit the characteristic
%     (fit parameters:  n_e,  T_e, V_p)
%--------------------------------------------------------------------------
% 12.1 Define fit range (ind_fit)
%--------------------------------------------------------------------------
if isempty(p_Vint) || ( isnan(p_Vint(1)) && isnan(p_Vint(2)) )
  i_fit0 = Vf.i;
%  i_fit1 = i_fit_ie_end;
  i_fit1 = round((i_fit_ie_end + i_fit_ie_end_2)/2);
else
  if isnan(p_Vint(1))  % Left: Vf
    i_fit0 = Vf.i;
  end
  if p_Vint(1)==-Inf  % Left:  1
    i_fit0 = 1;
  end
  if isnan(p_Vint(2)) % Right: Vp
%    i_fit1 = i_fit_ie_end;
    i_fit1 = round((i_fit_ie_end + i_fit_ie_end_2)/2);
  end
  if p_Vint(2)==Inf   % Right: 1
    i_fit1 = numel(x);
  end
  if p_Vint(1)>-Inf && p_Vint(1)<Inf  % Left: manual
    i_fit0  = findind(x,p_Vint(1));
  end
  if p_Vint(2)>-Inf && p_Vint(2)<Inf  % Right: manual
    i_fit1 = findind(x,p_Vint(2));
  end
end
% Define fit and weight indices:
ind_fit = i_fit0:i_fit1;


%--------------------------------------------------------------------------
% 12.2 Define Fit Function
%--------------------------------------------------------------------------
% Between the floating potential and the plasma potential an exponential 
% increase of the electron current is assumed:
% e0: elementary charge
% A_p: effective probe surface
% V: probe potential
% Phi: plasma potential
% Te_V: electron temperature in Volt (not in eV)
% I_e(V)   = I_esat * exp(-(Phi-V)/Te_V)
%   I_esat = -e0 *n_e *A_p *<v+>              %<v+> avg velocity to surface
%   -->    = -e0 *n_e *A_p *sqrt(e0*Te_V/(2*pi*m_e))
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% (A) n_e from fit of electron current and n_i from ion saturation current
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% I_e(V) = -e0*n_e*A_p*sqrt(e0*Te_V/(2*pi*m_e))*exp(-(Phi-V)/Te_V)
% Rearrange for fit function f(x):
%  I_e(V) = -e0*A_p*sqrt(e0/(2*pi*m_e))*n_e*sqrt(Te_V)*exp(-(Phi-V)/Te_V)
% Fit function with 3 fit parameters: a:n_e, b:Te_V, c:Phi
%    f(x) = v1*a*sqrt(b)*exp(-(c-x)/b)

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% (B) n_e=n_i of ion saturation current
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% I_e(V) = -e0*n_e*A_p*sqrt(e0*Te_V/(2*pi*m_e))*exp(-(Phi-V)/Te_V)
% From the ion saturation current with n_e=n_i
% and n_i = I_i,sat/(0.61*e0*A_p*sqrt(e0*Te_V/m_i))
% I_e(V) = -(1/0.61)*I_i,sat*sqrt(m_i/(2*pi*m_e))*exp(-(Phi-V)/Te_V)
% Rearrange for fit function f(x):
% Fit function with 2 fit parameters: a:Te_V, b:Phi
%     f(x) = v1*exp(-(a-x)/b)

%--------------------------------------------------------------------------
% 4.2.1 Determine Fit Start Parameters
%--------------------------------------------------------------------------
% Start parameter Te_V at Plasma Potential: Te_V = -I_e(V_p) / I_e'(V_p)
if isempty(p_fix.Te)
  Te_start = ye(iVp) / Dyesm(iVp);
else
  Te_start = p_fix.Te;
end
% Start parameter f(n_e)=const*n_e at plasma potential:
%   const*n_e = -I_e(V_p)/sqrt(Te_start)
if isempty(p_fix.ne)
  ne_start = ye(iVp)/sqrt(Te_start);
else
  ne_start = p_fix.ne;
end

% Switch between fit procedure (A) and (B)
switch p_ni
  case 'Ie'
    % OLD: fit_start = [ne_start  Te_start];
    fit_start = [ne_start  Te_start 0];
  case 'Iisat'
    fit_start = Te_start;
end
%--------------------------------------------------------------------------
% 4.2.2 Determine Fit Boundaries
%--------------------------------------------------------------------------
% Set fit limits [n_e, Te_V, Vp]
if isempty(p_fix.nelim)
  n0 =   0.1*ne_start;
  n1 = 100.0*ne_start;
else
  n0 = p_fix.nelim(1);
  n1 = p_fix.nelim(2);
end
  
if isempty(p_fix.Telim)
  T0 = 0.3*Te_start;
  T1 = 3.0*Te_start;
else
  T0 = p_fix.Telim(1);
  T1 = p_fix.Telim(2);
end

%--------------------------------
% 4.3 Fit
%--------------------------------
% Define Fit Function
% Rearrange for fit function f(x):
% I_e(V) = -e0*n_e*A_p*sqrt(e0*Te_V/(2*pi*m_e))*exp(-(Phi-V)/Te_V)
% I_e(V) = -e0*A_p*sqrt(e0/(2*pi*m_e))*n_e*sqrt(Te_V)*exp(-(Phi-V)/Te_V)
% Fit function with 3 fit parameters: a:f(n_e), b:Te_V, c:Phi
% f(x) = c(1)*sqrt(c(2))*exp(-(c(3)-x)/c(2))
% c(1): -e0*A_p*sqrt(e0/(2*pi*m_e))*n_e
% c(2): electron temperature
% Set the lsqcurvefit fitoptions
options = optimset('Display','off','TolFun',1e-8);

% In case upper boundaries are smaller than lower boundaries or numbers are
% not real, do not fit, and output NaN.
if isreal([n0 n1 T0 T1])
  if n1<n0 && T1>=T0
    % Find the best fit depending on the fit range and thus on V_plasma
    % Half-Length of Range variation (full range = 2*N2, see below):
    N2 = 220;
    indloop = [ind_fit (ind_fit(end)+1:ind_fit(end)+N2)];
    Rsqr = zeros(2*N2,1); Rstd = Rsqr; Rrat = Rsqr;  % Initialize variables
    for j=1:2*N2
      % Indices for the loop
      iloopfit = indloop(1:end-2*N2+j);
      % If loopfit is longer than x break up:
      if iloopfit(end)>numel(x); break; end
      % Define plasma potential (at end of fit region) an fit function
      vp = x(iloopfit(end));
      % Switch between fit procedure (A) and (B)
      switch p_ni
        case 'Ie'
          % ub = [n0 T1];   % ! n is negative (because electron current)
          % lb = [n1 T0];   % ! n is negative (because electron current)
          % ffun = @(c,x) c(1)*sqrt(c(2)).*exp(-(vp-x)./c(2));
          ub = [n0 T1 1];   % ! n is negative (because electron current)
          lb = [n1 T0 0];   % ! n is negative (because electron current)
          ffun = @(c,x) c(1)*sqrt(c(2)).*exp(-(vp-x)./c(2)) + c(3);
        case 'Iisat'
          ub = T1;
          lb = T0;
          ffun = @(c,x) c(1)*sqrt(c(2)).*exp(-(vp-x)./c(2)); % ***
      end
      % Fit whole electron current up to the plasma potential region
      % OLD: xx = x(iloopfit); yy = ye(iloopfit);
      xx = x(iloopfit); yy = y(iloopfit);
      % [ffit,~] = lsqcurvefit(@myfun_PlasmaFitIULangmuir_Ie, ...
      %   fit_start, xx, yy);
      ffit = lsqcurvefit(ffun, fit_start, xx, yy, lb, ub, options);
      yfit = ffun(ffit,xx);
      
      
      % % Chi Square Deviation:
        sef = sum( (yy-yfit).^2 );
        sdm = sum( (yy-mean(yy)).^2 );
      Rsqr(j) = 1 - sef/sdm;     % The larger the better
      Rstd(j) = std(yy-yfit);    % The smaller the better
      Rrat(j) = Rsqr(j)/Rstd(j); % The larger the better

      %%%%%%%%%%%%%%%%%%%%%
      % Control Plot
      %%%%%%%%%%%%%%%%%%%%%
      clf
      hold on
      plot(xx,  yy,  'bo')
      plot(xx, yfit, 'r')
      %%%%%%%%%%%%%%%%%%%%%

    end

    % Find fit with the lowest deviation (best plasma potential)
    [~, iRmax] = max(Rrat);
    % Define indices of fit region:
    iloopfit = indloop(1:end-2*N2+iRmax);
    % Define plasma potential (at end of region) an fit function
    vp = x(iloopfit(end));
    ffun = @(c,x) c(1)*sqrt(c(2)).*exp(-(vp-x)./c(2));
    switch p_ni
      case 'Ie'
        ub = [n0 T1];   % ! n is negative (because electron current)
        lb = [n1 T0];   % ! n is negative (because electron current)
      case 'Iisat'
        ub = [T1 Vp1];
        lb = [T0 Vp0];
    end
    % Fit whole electron current up to the plasma potential
      xx = x(iloopfit); yy = ye(iloopfit);
      ffit = lsqcurvefit(ffun,fit_start, xx, yy, lb, ub, options);
      yfit = ffun(ffit,xx);
      % Calculate Chi Square Deviation:
      sef = sum( (yy-yfit).^2 );
      sdm = sum( (yy-mean(yy)).^2 );
      iudata.Rsqr = 1 - sef/sdm;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      figeps(10,20,1,0,100); clf
      subplot(4,1,1);
      hold on
      plot(ye, 'ro'); plot(yfit, 'b');
        set(gca, 'ylim', [min(ye) max(ye)])
      ix = (0:numel(Rsqr)-1)+indloop(end-2*N2+1);
      subplot(4,1,2); plot(ix, Rsqr, 'k-');
      subplot(4,1,3); plot(ix, Rstd, 'b-');
      subplot(4,1,4); plot(ix, Rrat, 'm-');
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate fit function along the whole original x-vector
    Iefit = ffun(ffit,x);
    % Extract the confidence intervals
      a0 = [ffit(1) ffit(2)];
      [a,r,~,COVB] = nlinfit(xx,yy,ffun,a0);
    c = nlparci(a,r,'covar',COVB);
    
    % Calculate Output parameters (depending on fit method p_ni)
    switch p_ni
      case 'Ie'
        iudata.Te = ffit(2);
        iudata.Vp = vp;
        iudata.ni = Iisat_max/(0.61*e0*A_p*sqrt(e0*iudata.Te/m_i));
        iudata.ne = ffit(1)/(-e0*A_p*sqrt(e0/(2*pi*m_e)));
        iudata.neerr = iudata.ne*[1 1]; % ***
        iudata.Teerr = [c(2,1) c(2,2)];
        iudata.Vperr = iudata.Vp*[1 1]; % *** wdth Rrat-peak below 95% max
      case 'Iisat'
        iudata.Te = ffit(1);
        iudata.Vp = ffit(2);
        iudata.ni = Iisat_max/(0.61*e0*A_p*sqrt(e0*iudata.Te/m_i));
        iudata.ne = iudata.ni;
        iudata.neerr = iudata.ne*[1 1]; % ***
        iudata.Teerr = [c(1,1) c(2,1)];
        iudata.Vperr = iudata.Vp*[1 1]; % ***
    end
    iudata.Vf= Vf.V;
    % Calculate the Ion Flux: Gamma_i = n_i*<v_th,i>
    iudata.Fli = iudata.ni*sqrt(2*e0*iudata.Te/m_i);
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF CONTINUES BELOW AT "NO-OUTPUT" %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------
% 4.4 Calculate the floating potential from fit an linear I_i,sat
%-----------------------------------------------------------------
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
    iudata.Vf2 = x01;
  else
    % x02 is the zero crossing
    iudata.Vf2 = x02;
  end
  
else
  iudata.Vf2 = NaN;
end


% 5. Calculate the deviation of the fit from the rawdata
%--------------------------------------------------------
% Calculate the Normalized Integral Deviation from smooth Ie
integral1 = int_discrete(x(1:i_fit_ie_end), ...
  abs(ye(1:i_fit_ie_end)-Iefit(1:i_fit_ie_end)));
integral2 = int_discrete(x(1:i_fit_ie_end),  abs(ye(1:i_fit_ie_end)) );
%
% NID = int_V |I_mod(V) - I_raw(V)| dV / int_V I_raw(V) dV
iudata.NID = integral1/integral2;


% XXX  <<<--- Check fit of electron current between Vf and Vp
if p_chk
figeps(16,16,3,30,100); clf;
subplot(2,1,1)
hold on
  plot(1e3*ye, 'co-', 'Color', 0.8*[1 1 1])
  Ie = ffun(ffit,x(iloopfit));
  DIe = diff_discrete(x(iloopfit),Ie);
      yy = Dye(iloopfit);
      sef = sum( (yy-DIe).^2 );
      sdm = sum( (yy-mean(yy)).^2 );
      iudata.Rsqr2 = 1 - sef/sdm;
  plot(iloopfit, 1e3*Ie, 'r-', 'LineWidth', 2)
hold off
set(gca, 'ylim', 1e3*[min(data.current) max(data.current)])
hold off
set(gca, 'xlim', [1 ll])
mkplotnice('probe voltage (V)', 'fit I_e (mA)', 12, '-25', '-40');
puttextonplot(gca, [0 1], 10, -40, 'fit to I_e', 0, 14, 'r');
puttextonplot(gca, [0 0], 10, 90, ...
  ['R=' sprintf('%4.2f', iudata.Rsqr)], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 70, ...
  ['T_e=' sprintf('%4.2f', iudata.Te) 'eV'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 50, ...
  ['n_e=' sprintf('%4.2g', iudata.ne) 'm^{-3}'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 30, ...
  ['V_p=' sprintf('%4.2f', iudata.Vp) 'V'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 10, ...
  ['V_f=' sprintf('%4.2f', iudata.Vf2) 'V'], 0, 12, 'k');

subplot(2,1,2)
  hold on
  plot(x,1e3*Dyesm, 'bo-', 'Color', 0.8*[1 1 1])
  plot(x(iloopfit),1e3*DIe, 'r-', 'Linewidth', 2)
  hold off
  set(gca, 'xlim', [x(1) x(end)])
  mkplotnice('probe voltage (V)', 'dI_e/dV (mA/V)', 12, '-25', '-40');
  puttextonplot(gca, [0 0], 10, 90, ...
    ['R=' sprintf('%4.2f', iudata.Rsqr2)], 0, 12, 'k');

save fitdata_a11_0001-48.mat ye iloopfit Iefit data iudata Ie DIe Dye par_ysm x
  
% figeps(17,12,10,40,0); plot(x,y);
% input('>>> Press Any Key <<<')

end
% XXX  <<<--- Check ion saturation current

%%%%%%%%%%%%%%%
% "NO-OUTPUT" %
%%%%%%%%%%%%%%%
  else   % Found at least one lower boundary larger than a upper boundary.
   disp('No fit: Lower and upper boundaries are not logic.')
    iudata.ne  = NaN; iudata.neerr = NaN;
    iudata.Te  = NaN; iudata.Teerr = NaN;
    iudata.Vp  = NaN; iudata.Vperr = NaN;
    iudata.Vf  = NaN;
    iudata.Vf2 = NaN;
    iudata.Fli = NaN;
    iudata.Rsqr= NaN;
    iudata.NID = NaN;
    return
  end
else  % Found at least one imaginary number, i.e. a fit makes no sense.
  disp('No fit: Boundary values of ne, Te and Vp are not real numbers.')
    iudata.ne  = NaN; iudata.neerr = NaN;
    iudata.Te  = NaN; iudata.Teerr = NaN;
    iudata.Vp  = NaN; iudata.Vperr = NaN;
    iudata.Vf  = NaN;
    iudata.Vf2 = NaN;
    iudata.Fli = NaN;
    iudata.Rsqr= NaN;
    iudata.NID = NaN;
  return
end

end