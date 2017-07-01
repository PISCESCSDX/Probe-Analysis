function [iudata, raw] = PlasmaFitIULangmuir_v24(data, checkonoff, ...
  fixpar, Vint, vol_avg, DIeexp_level)
%==========================================================================
%function [iudata] = PlasmaFitIULangmuir_v24(data, checkonoff, fixpar, ...
% Vint, vol_avg)
% Last Change: 2012-04-13 10:15 C. Brandt, San Diego
% - Test version: try to improve fit range detection
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
% DIeexp_level: level range of slope (strong influence on T_e!!)
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

% # INPUT: Level-Range of Slope (parameter with strong influence on T_e!!)
if nargin <6
  DIeexp_level = 0.95;
end

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
  iudata = []; raw = [];
  return
end


% 3. Calculate the ion saturation current
%--------------------------------------------------------------------------
[PolyIsat, ~, i_isat, raw.iVf] = subPlasmaIsatPoly_v2(data);


% 4. Exclusion Test: Ion saturation current was found
%--------------------------------------------------------------------------
if isempty(PolyIsat)
  disp('Could not detect I_i,sat.')
  iudata = []; raw = [];
  return
end
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXX --->>> Check ion saturation current detection:
if checkonoff
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
  iudata = []; raw = [];
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
  zmin = min(zero_peaks);        zmax = max(zero_peaks);
  pkpos = zero_peaks > zmax/2;   pkpos = indx(pkpos);
  pkneg = zero_peaks < zmin/2;   pkneg = indx(pkneg);
  % Delete all negative peaks below the first positive peak:
  if pkneg(1)<=pkpos(1)
    pkneg = pkneg(2:end);
  end
  % Negative Peaks may be missing one peak at the end: if so add one at end
  if numel(pkneg)-numel(pkpos) == -1
    pkneg(end+1) = ll;  % ll = length(y); Defined at Start
  end
  % Calculate the difference between the start and end peaks:
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
 disp('Did not find the start of the exponential growth! No Ev!')
  iudata = []; raw = [];
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
%     (
  ii = imax >= i_Ie_endofnoise-2*pts; ii = imax(ii);
  ii2= ii < (numel(yelog)-ceil(numel(ind_Ie)/2)); ii2 = ii(ii2);
  imax = ii2(1);

% Look at points within a certain level-range around the local maximum
% of the 1st derivative (find the exponential slope region)
% Left and Right Indices
ind = Dyelog_sm < DIeexp_level*Dyelog_sm(imax);
  i_lower = indx(ind);
  ile = i_lower<imax; ile = i_lower(ile); ile = ile(end);
  iri = i_lower>imax; iri = i_lower(iri); iri = iri(1);
% Define indices of exponential electron current region
ind_Ieexp = ile:iri;
% Index where maximum of D(Ie) is detected (= minimum of Te)
i_dIe_max = imax;


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Check Plots: (1) rawdata log(I_e), (2) smoothed log(I_e)
if checkonoff
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


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Histogram Plots
%--------------------------------------------------------------------------
[hy,hx] = hist(yelog,round(ll/2));
  yout = filtspline(hx, hy, 2, 0.9);
  [~,~,~,imin] = extrema(yout);
  imin = sort(imin);
  % Check if minima are present anyway
  if ~isempty(imin)
    % If the desired minimum is at the end:
    if imin(end)==length(hx)
        % If more than one minimum is detected:
        if length(imin)>1
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
  
if checkonoff
figeps(22,10,2,0,0); clf;
subplot(2,2,2)
hist(y,25)

subplot(2,2,4)
hold on
plot(hx,hy);
  plot(hx,yout, 'r')
hold off

subplot(2,2,3)
plot(yelog)
line(i_fit_ie_end*[1 1], [min(yelog) max(yelog)], 'color', 'r')

subplot(2,2,1)
plot(y)
line(i_fit_ie_end*[1 1], [min(y) max(y)], 'color', 'r')
end
% XXX <<<---
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


% 7.3 Inflection point of characteristic (Plasma Potential)
%     Start at maximum of 1st Derivative, go to left, find zero crossing
i_start = i_dIe_max;
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
    i_Ie_inflection = ctr;
    status = 1;
  end
end


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% (1) log(Ie), (2) D(log(Ie)), (3) D2(log(Ie))
if checkonoff
ind1 = 1:i_start-1;
ind2 = i_start:length(yelog);

figeps(10,14,3,20,100);clf;
subplot(3,1,1); 
hold on
plot(ind1,yelog(ind1),'k-o')
plot(ind2,yelog(ind2),'b-v')
line(i_dIe_max*[1 1] , 1.1*[min(yelog(ind)) max(yelog(ind))], 'Color', 'b')
line(i_Ie_inflection*[1 1],[min(yelog(ind)) max(yelog(ind))], 'Color', 'r')
hold off
set(gca, 'xlim', [ind2(1)-10 ind2(end)])
set(gca, 'ylim', ...
  [min(yelog(ind2(1)-10:ind2(end))) max(yelog(ind2(1)-10:ind2(end)))])
mkplotnice('', '', 12, '0', '0');

subplot(3,1,2);
hold on
plot(ind1,Dyelog_sm(ind1),'k-o')
plot(ind2,Dyelog_sm(ind2),'b-v')
hold off
set(gca, 'xlim', [ind2(1)-10 ind2(end)])
set(gca, 'ylim', [min(Dyelog_sm(ind2(1)-10:ind2(end))) ...
  max(Dyelog_sm(ind2(1)-10:ind2(end)))])
mkplotnice('', '', 12, '0', '0');

subplot(3,1,3);
hold on
plot(ind1,D2yelog_sm(ind1), 'k-o')
plot(ind2,D2yelog_sm(ind2), 'b-v')
plot(le, D2yelog_sm(le) , 'r-*')
plot(ri, D2yelog_sm(ri) , 'c-*')
hold off
set(gca, 'xlim', [ind2(1)-10 ind2(end)])
set(gca, 'ylim', [min(D2yelog_sm(ind2(1)-10:ind2(end))) ...
  max(D2yelog_sm(ind2(1)-10:ind2(end)))])
mkplotnice('', '', 12, '0', '0');
end
% XXX <<<---
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


% 8. Define rough estimate of plasma potential
%--------------------------------------------------------------------------
iVp = i_dIe_max;
%iVp = i_Ie_inflection;


% 9. Define rough estimate of plasma potential (for fit-start-parameters)
%--------------------------------------------------------------------------
Dye = diff_discrete(x, ye);
% 3.1 smooth the derivative (empiric amount of smoothing points)
Dyesm = filtmooth(Dye, 3*par_ysm);



[~,~,~,imin] = extrema(Dyesm);
imin2 = sort(imin);
dimin = round( (imin2(end)-imin2(1)) / numel(imin2) );
Dyesm2= filtspline(x, Dye, 2*dimin, 0.1);
[~,~,~,imin] = extrema(Dyesm2);
while imin(1)==1 || imin(1)==ll
 imin = imin(2:end);
end
i_fit_ie_end_2 = imin(1);

%TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
% Figures showing the different important indices:
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
  ['i-Ieexp-start: ' num2str(i_Ieexp_start)], 0, 8, 'b');
puttextonplot(gca, [0 0], 5, 100-2*dy, ...
  ['i-Ie-inflection: ' num2str(i_Ie_inflection)], 0, 8, 'r');
puttextonplot(gca, [0 0], 5, 100-3*dy, ...
  ['i-fit-ie-end ' num2str(i_fit_ie_end)], 0, 8, 'm');

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

figure;
plot(imin2, 'ro')

figure;
subplot(2,1,1)
hold on; plot(Dyesm, 'ro-'); plot(filtspline(x, Dye, 2*dimin, 0.1), 'b-');
hold off
subplot(2,1,2)
plot(diff_discrete(x, filtspline(x, Dye, 2*dimin, 0.1)), 'mo-')
%TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
%   
% 
% i0 = ind_Ieexp(1);
% i_start  = i0+1;
% 
% dev3= zeros(1,ll-i_start+1);
% for i1=i_start:ll
%   [p, s] = polyfit((i0:i1)', yelog(i0:i1), 1);
%   dev3(i1-i_start+1) = sum(abs( yelog(i0:i1) - polyval(p,i0:i1)' ))/(i1-i0);
% end
% 
% figure
% hold on
% plot(i_start:ll, filtspline(i_start:ll, dev3,dimin,0.9),'bO-')
% hold off
% %i_fit_ie_end = :
% 
% 
% 
% i_start  = i0-30;
% L=ll-i_start;
% d1 = 3;
% for d = d1:floor(L/4)
%   for j=i_start:i_start+50
%     ind = j:j+d-1;
%     testf = yelog;
% %    testf = filtspline(x,yelog_raw,1,0.7);
%     [p, s] = polyfit(ind', testf(ind), 1);
%     mat(d-d1+1, j-i_start+1) = ...
%       sign(p(1))*sum(abs( testf(ind) - polyval(p,ind)' )) / d;
%   end
% end
% 
% % Make column of mat NaN if mat(1,i)<0 
% for j=i_start:i_start+50
%   if mat(1, j-i_start+1)<0
%     mat(:,j-i_start+1)=NaN;
%   end
% end
% figeps(20,10,4,40,30);
% subplot(1,2,1)
% surf(mat);
% 
% mat = mat./matmax(mat);
% xo = i_start:i_start+50;
% yo = d1:floor(L/4);
% [xn,yn,An] =interp_matrix(xo,yo,mat,[2*numel(xo) 2*numel(yo)]);
% 
% subplot(1,2,2)
% %plotcontour(i_start:i_start+10, d1:floor(L/4), mat, 20, '', '', 12)
% plotcontour(xn, yn, An, 40, 'i', 'd', 12)
% colorbar
% colormap pastell


if isempty(fixpar.Vp)
  Vp_start = x(iVp);
else
  Vp_start = fixpar.Vp;
  ind = find(x<Vp_start);
  iVp = ind(end);
end

% 9.1 Excusion test: If no Plasma Potential detected, finish this m-file.
if isempty(iVp)
  disp('No plasma potential could be found! No Output!')
  iudata = [];
  return
end


% 10. If necessary correct floating potential to left
%--------------------------------------------------------------------------
if raw.iVf<iVp
  iVf = raw.iVf;
else
  iVf = iVp-1;
end


%==========================================================================
% 12. Fit the characteristic
%     (fit parameters:  n_e,  T_e, V_p)
%--------------------------------------------------------------------------
% 12.1 Define fit range (ind_fit)
%--------------------------------------------------------------------------
if isempty(Vint) || ( isnan(Vint(1)) && isnan(Vint(2)) )
  i_fit0 = iVf;
  i_fit1 = round((i_fit_ie_end + i_fit_ie_end_2)/2);
  i_w0 = iVf;
  i_w1 = i_fit_ie_end - iVf + i_w0;
else
  if isnan(Vint(1))  % Left: Vf
    i_fit0 = iVf;
      i_w0 = 1; %----------------- Start of Weights (in exponential region)
  end
  if Vint(1)==-Inf  % Left:  1
    i_fit0 = 1;
      i_w0 = iVf; %--------------- Start of Weights (in exponential region)
  end
  if isnan(Vint(2)) % Right: Vp
%     i_fit1 = i_fit_ie_end;
    i_fit1 = round((i_fit_ie_end + i_fit_ie_end_2)/2);
      i_w1 = i_fit1-iVf+i_w0; %End of Weights (in exponential region)
  end
  if Vint(2)==Inf   % Right: 1
    i_fit1 = numel(x);
      i_w1 = i_fit_ie_end-iVf+i_w0; %End of Weights (in exponential region)
  end
  if Vint(1)>-Inf && Vint(1)<Inf  % Left: manual
    i_fit0  = findind(x,Vint(1));
    i_w0    = [];
  end
  if Vint(2)>-Inf && Vint(2)<Inf  % Right: manual
    i_fit1 = findind(x,Vint(2));
    i_w1   = [];
  end
end
% Define fit and weight indices:
ind_fit = i_fit0:i_fit1;
ind_w   = i_w0:i_w1;

%--------------------------------------------------------------------------
% Load natural physical constants, probe area
%--------------------------------------------------------------------------
natconst
Ap = 2*pi*data.r*data.l;
%--------------------------------------------------------------------------
% 12.2 Define Fit Function
%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
% 4.2 Determine Fit Start Parameters
%--------------------------------------------------------------------------
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
  n0 =  0.1*ne_start;
  n1 = 10.0*ne_start;
else
  n0 = fixpar.nelim(1);
  n1 = fixpar.nelim(2);
end

if isempty(fixpar.Telim)
  T0 =  0.1*Te_start;
  T1 = 10.0*Te_start;
else
  T0 = fixpar.Telim(1);
  T1 = fixpar.Telim(2);
end

if isempty(fixpar.Vplim)
  % Vp will be not varied
  Vp0 = Vp_start;
  Vp1 = Vp_start;
else
  Vp0 = Vp_start + fixpar.Vplim(1);
  % *** Control the influence of an increase of the upper Vp-boundary
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
    % Increase Weights of Fit in the exponential region:
    opts.weights = ones(1,numel(ind_fit));
    opts.weights(ind_w) = 1000;
    % Fit whole electron current up to the plasma potential
    [ffit, gof] = fit(x(ind_fit), ye(ind_fit), fitfct, opts);
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
integral1 = int_discrete(x(1:i_fit_ie_end), ...
  abs(ye(1:i_fit_ie_end)-Iefit(1:i_fit_ie_end)));
integral2 = int_discrete(x(1:i_fit_ie_end),  abs(ye(1:i_fit_ie_end)) );
%
% NID = int_V |I_mod(V) - I_raw(V)| dV / int_V I_raw(V) dV
iudata.NID = integral1/integral2;


% XXX  <<<--- Check fit of electron current between Vf and Vp
if checkonoff
figeps(16,16,3,30,100); clf;
subplot(2,1,1)
hold on
  plot(1:i_fit1, 1e3*ye(1:i_fit1), 'co-')
  plot(1:i_fit1, 1e3*ffit(x(1:i_fit1)), 'r')
  %
  plot(i_fit1+1:numel(x), 1e3*ye(i_fit1+1:end), 'bo-')
hold off
set(gca, 'ylim', 1e3*[min(data.current) max(data.current)])
hold off
set(gca, 'xlim', [1 ll])
mkplotnice('probe voltage (V)', 'fit I_e (mA)', 12, '-25', '-40');
puttextonplot(gca, [0 1], 10, -40, 'fit to I_e', 0, 14, 'k');
puttextonplot(gca, [0 0], 10, 90, ...
  ['R=' sprintf('%4.2f', iudata.Rsqr)], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 70, ...
  ['T_e=' sprintf('%4.2f', iudata.Te) 'eV'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 50, ...
  ['n_e=' sprintf('%4.2g', iudata.ne) 'm^{-3}'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 30, ...
  ['V_p=' sprintf('%4.2f', iudata.Vp) 'V'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 10, ...
  ['V_f=' sprintf('%4.2f', iudata.Vffit) 'V'], 0, 12, 'k');

subplot(2,1,2)
  plot(1e3*Dyesm, 'bo-')
  set(gca, 'xlim', [1 ll])
  mkplotnice('probe voltage (V)', 'dI/dV (mA/V)', 12, '-25', '-40');
  
% input('>>> Press Any Key <<<')
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