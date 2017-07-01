function iudata = PlasmaFitIUDemidov(data, p_chk, p_fix, p_Vint, ...
  p_Vavg, p_IsRg, p_fitgoal, p_scanstepsize)
%==========================================================================
% function iudata = PlasmaFitIUDemidov(data, p_chk, p_fix, p_Vint, ...
%   p_Vavg, p_IsRg, p_fitgoal, p_stepsize)
%
% Last Change: 2013-08-09 16:19 C. Brandt, San Diego
% This is alpha version: v2.0.
% This version includes a manual fitting option asked after every fit.
% For Details look into the file: 2012-04-19_IUFit.pdf
%--------------------------------------------------------------------------
% PlasmaFitIUDemidov calculates the plasma parameters n_e (electron 
% density), T_e (electron temperature), V_p (plasma potential) and the V_f 
% (floating potential) from a current-voltage characteristic using a 
% cylindrical electrode (Langmuir probe).
% === External m-files ===
% EXTREMUM, FIGEPS, PUTTEXTONPLOT, MKPLOTNICE
% === Problems ===
% If ion saturation current is found to be increasing no floating
% potential may be found.
%
% Sub-m-files: subPlasmaIsatPoly.m
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
%           [-40, 20]: [-40, 20]
%           [-Inf, 20]: [left end, 20]
%           [-40, Inf]: [-40, right end]
%           [-Inf, Inf]: [left end, right end]
%     p_Vavg: smooth characteristic over voltage range (default 3V)
%              (used for smoothing the electron current)
%     p_IsRg: ion saturation current region (in Volt)
%     p_fitgoal: (1) best fit for I_e or (2) dI_e/dV
%     p_scanstepsize: stepsize of fit range scan (in sample steps)
%OUT: iudata.ne: electron density (m^-3)
%     iudata.ni:      ion density (m^-3)
%     iudata.Te: electron temperature (eV)
%     iudata.Vp: plasma potential (V)
%     iudata.Vf: floating potential (V)
%--------------------------------------------------------------------------
% EX:
% p_chk=1; p_fix.Te=[];p_fix.Telim=[];p_fix.ne=[];...
% p_fix.nelim=[];p_fix.Vp=[];p_fix.Vplim=[-2 2];p_Vint=[-Inf, NaN];...
% p_Vavg=1;
% iudata=PlasmaFitIUDemidov_v11(data,p_chk,p_fix,p_Vint,p_Vavg);
%==========================================================================


% # INPUT: Ion saturation current range in Volt
if nargin < 7 || isempty(p_fitgoal)
  p_fitgoal = input('Best fit for (1) I_e (2) dI_e/dV (3) sum? ');
end

% # INPUT: Ion saturation current range in Volt
if nargin < 6
  p_IsRg = input('ion saturation current range (in Volt): ');
end

% Set default voltage smoothing range
if nargin <5
  p_Vavg = 3;
end

% Set default voltage cut interval
if nargin <4
  p_Vint = [-Inf, +Inf];
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
  r_p = data.r;
  L_p = data.l;
A_p = 2*pi *r_p *L_p;
% Ion mass
m_i = data.m_ion*u0;
% pre-factor for fitting
ne_prefac = -(3*sqrt(pi)/4) *( r_p*log(pi*L_p/(4*r_p))/A_p ) *(data.B/e0);

iudata.probe.L = L_p;
iudata.probe.R = r_p;
iudata.probe.A = A_p;
iudata.plasma.B= data.B;
    
% 1. Cut Characteristic (both sides 5%)
%--------------------------------------------------------------------------
% Extract the characteristic and the probe parameters
% Cut margin: 2% both sides
L5 = round(0.02*length(data.voltage));
if L5==0
  disp('I-V characteristic is too short -> no automatic evaluation')
  iudata = [];
  return
end
ind = L5:length(data.voltage)-L5;
data.voltage = data.voltage(ind);
data.current = data.current(ind);

x = data.voltage;
xind = 1:numel(x);
y = data.current;
  dx = (x(end)-x(1)) / numel(x);


% 2. Exclusion Test: I-V characteristic is asymmetric
%--------------------------------------------------------------------------
% First simple test whether data may be Langmuir characteristic
% (average of left side must be larger than average of right side)

% if mean(y(1:round(length(x)/2))) <= mean(y(round(length(x)/2):end))
%   disp(['I-V characteristic is not asymmetric enough, i.e. most ' ...
%       'probable not evaluable.'])
%   iudata = [];
%   return
% end


% 3. Calculate the ion saturation current
%--------------------------------------------------------------------------
[PolyIsat, ~, i_isat, Vf] = subPlasmaIsatPoly(data, p_IsRg);


% 4. Exclusion Test: Ion saturation current was found
%    (store Iisat_max for n_i calculation from Bohm criterion)
%--------------------------------------------------------------------------
if isempty(PolyIsat)
  disp('Could not detect I_i,sat.')
  iudata = [];
  return
else
  % Determine maximum Ion saturation current
  Iisat_max = mean(polyval(PolyIsat, x(i_isat) ));
end
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXX --->>> Check ion saturation current detection:
% if p_chk
%   figeps(12,8,1,0,100); clf;
%   axes('position', [0.20 0.16 0.78 0.72])
%   title('I_{sat} voltage limits');
%   hold on
%     % Plot Rawdata
%     plot(x,y*1e3,'ko-')
%     ind_isat = i_isat(1):i_isat(2);
%     plot(x(ind_isat),y(ind_isat)*1e3,'bo-')
%     xIs = min(x):max(x);
%     yIs = polyval(PolyIsat, xIs);
%     % Plot Ion saturation current fit
%     plot(xIs,yIs*1e3,'r-')
%   hold off
%   set(gca, 'xlim', [data.voltage(1) data.voltage(end)], ...
%   'ylim', 1e3*[min(data.current) max(data.current)])
%   mkplotnice('probe voltage (V)', 'probe current (mA)', 12, '-25', '-50');
% end
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
% TEST Feb 14 2013
%  Smooth logarithm of Ie over 'p_Vavg' Volts (empirically found to work)
L_sm  = ceil(p_Vavg/dx);
%
dump.ind = i_isat(2):length(x);
dump.diff = diff_discrete(x(dump.ind), smooth( ye(dump.ind) , L_sm) );
%
[~, imax] = min( dump.diff );
% 
if imax==length(yelog(i_isat(2):end)) || imax==1
 disp('The maximum of smoothed DI is at the left or right margin! No Ev!')
  iudata = [];
  return
end



% 12. Fit Characteristic (fit parameters:  n_e, T_e, V_p)
%--------------------------------------------------------------------------
% 12.1 Define fit range (ind_fit)
%--------------------------------------------------------------------------
if p_Vint(1)==-Inf  % Left:  1
  i_fit0 = 1;
end
if p_Vint(2)==Inf   % Right: end
  i_fit1 = numel(x);
end
if p_Vint(1)>-Inf && p_Vint(1)<Inf  % Left: manual
  i_fit0  = findind(x,p_Vint(1));
end
if p_Vint(2)>-Inf && p_Vint(2)<Inf  % Right: manual
  i_fit1 = findind(x,p_Vint(2));
end
% Define fit and weight indices:
xindfit = i_fit0:i_fit1;

% Define xfit and yefit as the to be fitted range
xfit  = x(xindfit);
yefit = ye(xindfit);

% 12.2 Determine Fit Start Parameters
%--------------------------------------------------------------------------
% Calculate Derivative again used for fitting:
  Dye   = diff_discrete(x, ye);
  Dyefit= diff_discrete(xfit, yefit);
  Dyesm = filtmooth(Dyefit, par_ysm);
% Start parameter Te_V at Plasma Potential: Te_V = -I_e(V_p) / I_e'(V_p)
if isempty(p_fix.Te)
  [~, iVp] = min(Dyesm);
  Te_start = yefit(iVp) / Dyesm(iVp);
else
  Te_start = p_fix.Te;
end
% Start parameter at plasma potential: cne_start = ne_start/ne_pre_fac;
if isempty(p_fix.ne)
  ne_start = yefit(iVp)/sqrt(Te_start) /(-e0*A_p*sqrt(e0/(2*pi*m_e)));
  cne_start = ne_start / ne_prefac;
else
  ne_start = p_fix.ne;
  cne_start = ne_start / ne_prefac;
end

% Switch between fit procedure (A) and (B)
fit_start = [cne_start  Te_start];

%--------------------------------------------------------------------------
% 4.2.2 Determine Fit Boundaries
%--------------------------------------------------------------------------
% Set fit limits [n_e, Te_V, Vp]
if isempty(p_fix.nelim)
  cn1 = 1e-6*cne_start;
  cn0 = 1e6*cne_start;
else
  cn0 = p_fix.nelim(1) /ne_prefac;
  cn1 = p_fix.nelim(2) /ne_prefac;
end
  
if isempty(p_fix.Telim)
  T0 = 1e-2*Te_start;
  T1 = 1e+2*Te_start;
else
  T0 = p_fix.Telim(1);
  T1 = p_fix.Telim(2);
end

%==========================================================================
% 4.3 Fit Procedure
%--------------------------------------------------------------------------
% Definition of Fit Function [Demidov1999]:
% T_eV: electron temperature in Volt (not eV!)
% R_p: probe radius
% L_p: probe length
% r_ce: electron Larmor radius: r_ce = m_e*v_the/(e0*B)
% v_the = sqrt(2*e0*T_eV/m_e): thermal velocity of electrons
% e0: elementary charge
% A_p: effective probe surface
% V: probe potential
% Phi: plasma potential
%
% (Dir) Perpendicular to B-Field:
% f(eV) = -3*(m_e^2) *R_p *log(pi*L_p/(4*R_p)) / ...
%          (8*pi*(e0^3)*r_ce*V) * (dj_e/d_V)
% Maxwell distribution:
% f(eV) = n_e * (m_e/(2*pi*e0*T_eV))^(3/2) * exp(-(Phi-V)/T_eV)
%
% Rearrange for fit function f(x):
% Fit function with 3 fit parameters: a:f(n_e), b:Te_V, c:Phi
% f(x) = c(1)*sqrt(c(2))*exp(-(c(3)-x)/c(2))
%
%   (Dir) Parallel to B-Field:
%   f(eV) = -3*(m_e^2) *L_p / (32*(e0^3)*V) * (dj_e/d_V)


% Set the fit options for LSQCURVEFIT
options = optimset('Display','off','TolFun',1e-8);
% In case upper boundaries are smaller than lower boundaries or numbers are
% not real, do not fit, and output NaN.
if isreal([cn0 cn1 T0 T1])
  if cn1>=cn0 && T1>=T0
    % Find the best fit depending on the fit range and thus on V_plasma
    % Half-Length of Range variation (full range = 2*N2, see below):
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    % Fit Loop: Length of fit range is scanned
    %======================================================================
    % Deviations between rawdata and fits to dIe/dV and Ie are calculated
    fitscanrange = 2/3;  % # Input: From right to left (>0 && <1)
    jvec = 1:p_scanstepsize:round(fitscanrange *length(xindfit));
    ctr_jvec = 1:length(jvec);
    % Initialize Statistic Variables
    Rsqr1 = zeros(length(jvec),1); Rstd1 = Rsqr1;
    Rsqr2 = Rsqr1; Rstd2 = Rstd1; Vploop = Rsqr1;    
    
          %%%%%%%%%%%%%%%%
          if p_chk
          figure(1); clf;
          subplot(2,1,1);
          hold on
          plot(xfit,yefit,'bo');
          hp = plot(0,0,'r');
          hold off
          mkplotnice('V (V)', 'I (A)', 12, '-20', '-30');
          subplot(2,1,2);
          hold on
          hp2 = plot(1,1);
          hp3 = plot(1,1);
          hp4 = plot(1,1);
          title('Rsqr (r:R(DI_e), b:R(I_e), g: Sum)')
          hold off
          end
          %%%%%%%%%%%%%%%%

    ctr_j = 0;  % Define extra counter because steps of jvec can be >1
    ffit{length(jvec)} = [];
    for j = jvec
      ctr_j = ctr_j + 1;
      % Indices for the loop
      iloopfit = xindfit(1:end-j);
      % Define plasma potential (at end of fit region)
        vp = x(iloopfit(end));
        Vploop(ctr_j) = vp;
      % Define Fit Limits and Fitfunction
      ub = [cn1 T1];  lb = [cn0 T0];
      ffun = @(c,x) c(1)*(vp-x)./c(2).*exp(-(vp-x)./c(2));
      %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      % Fit whole electron current up to the plasma potential
      xx = x(iloopfit); yy = Dye(iloopfit);
      ffit{ctr_j} = lsqcurvefit(ffun, fit_start, xx, yy, lb, ub, options);
      %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      % Chi Square Deviation Rsqr1 and Standard Deviation of dIe/dV
      Dyefit = ffun(ffit{ctr_j},xx);
          sef = sum( (yy-Dyefit).^2 );
          sdm = sum( (yy-mean(yy)).^2 );
        Rsqr2(ctr_j) = 1 - sef/sdm;       % The larger the better
        Rstd2(ctr_j) = std(yy-Dyefit);    % The smaller the better
      % Chi Square Deviation Rsqr2 and Standard Deviation of Ie
      % Calculate Fit to Electron Current from fit to Derivative:
      [~, Iefit] = int_discrete(xx, ffun(ffit{ctr_j},xx) );
      yy2 = ye(iloopfit);
          sef = sum( (yy2-Iefit).^2 );
          sdm = sum( (yy2-mean(yy2)).^2 );
        Rsqr1(ctr_j) = 1 - sef/sdm;       % The larger the better
        Rstd1(ctr_j) = std(yy2-Iefit); % The smaller the better
        
          %%%%%%%%%%%%%%%%
          if p_chk
          subplot(2,1,1);
          hold on
          set(hp, 'visible', 'off')
          hp = plot(xx,Iefit,'r','LineWidth',2); drawnow
          hold off

          subplot(2,1,2); 
          set(hp2,'visible', 'off')
          set(hp3,'visible', 'off')
          set(hp4,'visible', 'off')
          hold on
          hp2 = plot(jvec,Rsqr1','b');
          hp3 = plot(jvec,Rsqr2','r');
          hp4 = plot(jvec,Rsqr1'+Rsqr2','g');
          hold off
          drawnow
          end
          %%%%%%%%%%%%%%%%
        
    end
    hold off


    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    % Take Parameters of Optimal Fit
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    % Combined Parameter: Chi of dIe/dV and of Ie
    % Find fit with the lowest deviation (best plasma potential)
    switch p_fitgoal
      case 1               % Best Fit for I_e
        choice_R = Rsqr1;
      case 2               % Best Fit for dI_e/dV
        choice_R = Rsqr2;
      case 3               % Best Fit for I_e && dI_e/dV
        choice_R = Rsqr1.*Rsqr2;
    end
    [ymax, iRmax] = max(choice_R);
      k = choice_R > 0.95*ymax;
        hev = ctr_jvec(k);
        if iRmax==1 || iRmax==numel(xindfit)
          disp(['PlasmaFitIUDemidov Error: Optimal Plasma Potential ' ...
            'outside of scanned range.'])
          iudata = [];
          return
        else
          kL = hev>=iRmax; kL = hev(kL);
          kS = hev<=iRmax; kS = hev(kS);
          if ~isempty(kL) && ~isempty(kS)
            kS = kS(1);
            kL = kL(end);
          else
     disp('No fit: Did not find plasma potential, and hence no opt. fit.')
     disp('Try:    Maybe increase the fit scan range!')
            iudata.ni  = NaN;
            iudata.ne  = NaN; iudata.neerr = NaN;
            iudata.Te  = NaN; iudata.Teerr = NaN;
            iudata.Vp  = NaN; iudata.Vperr = NaN;
            iudata.Vf  = NaN;
            iudata.Vf2 = NaN;
            iudata.Fli = NaN;
            iudata.RsqrIe  = NaN;
            iudata.RsqrDIe = NaN;
            iudata.EDF.y=NaN;
            return
          end
        end

    iudata.Vperr = [Vploop(kL) Vploop(kS)]; % Wdth Rrat1-peak below 95% max
    iudata.RsqrIe  = Rsqr1(iRmax);
    iudata.RsqrDIe = Rsqr2(iRmax);
    % ----------------------------- Calculate the Optimal Fit of DIe (=Dye)
      xifitres = xindfit(1:end-jvec(iRmax));
      xx = x(xifitres);
      yy = Dye(xifitres);
    % Newly define Plasma Potential
      vp = x(xifitres(end));
    % Newly define Fit function
      ffun = @(c,x) c(1)*(vp-x)./c(2).*exp(-(vp-x)./c(2));
    % Calculate the Optimal Fit of DIe (=Dye)
      Dyefit = ffun(ffit{iRmax},xx);
    % Calculate the Optimal Fit For Ie
      [~, Iefit] = int_discrete(xx, ffun(ffit{iRmax},xx) );
    % Extract Confidence Intervals
      a0 = [ffit{iRmax}(1) ffit{iRmax}(2)];
      options = statset('FunValCheck', 'Off');
      [a,r,~,COVB] = nlinfit(xx,yy,ffun,a0,options);
    c = nlparci(a,r,'covar',COVB);
    
du.i=0;
du.i=du.i+1;iudata.fitfun{du.i}='natconst';
du.i=du.i+1;iudata.fitfun{du.i}='vp = Vp.n';
du.i=du.i+1;iudata.fitfun{du.i}='ne_prefac = -(3*sqrt(pi)/4)*(r_p*log(pi*L_p/(4*r_p))/A_p )*(data.B/e0);';
du.i=du.i+1;iudata.fitfun{du.i}='ffun=@(c,x) c(1)*(vp-x)./c(2).*exp(-(vp-x)./c(2))';
du.i=du.i+1;iudata.fitfun{du.i}='x=-150:150;';
du.i=du.i+1;iudata.fitfun{du.i}='DIe: y=ffun([ne.n/ne_prefac Te.n],x)';
du.i=du.i+1;iudata.fitfun{du.i}='Ie: [~,y]=int_discrete(x,DIe);';
du.i=du.i+1;iudata.fitfun{du.i}='plot(x,Ie,''b'',x,DIe,''r'');';
      
      % Extract Electron Energy Distribution Function
      % (for the fit function, for the smoothed data and the raw data)
      % [Demidov1999 Eq. (9) used with the energy dependend Electron Larmor
      % radius: R_Le = sqrt(2*m_e/(eB^2) * (V_plasma - V_probe))]
      %
      % Electron Energy Distribution Function (raw data)
      % (*** Check the results!)
      iudata.EEDF.Demidov.raw.x = xx;
      iudata.EEDF.Demidov.raw.y = -3/(8*pi) * ...
            (m_e^2/(e0^3)) * ...
            r_p*log(pi*L_p/(4*r_p)) *data.B ./ ((vp-xx).^(3/2)) .* yy;
      % Electron Velocity Distribution Function (Fit)
      % (*** Check the results!)
      iudata.EEDF.Demidov.fit.x = xx;
      iudata.EEDF.Demidov.fit.y = -3/(8*sqrt(2)*pi) * ...
            (m_e^(3/2)/(e0^(5/2))) * ...
            r_p*log(pi*L_p/(4*r_p)) *data.B ./ ((vp-xx).^(3/2)) .* Dyefit;

    % Druyvesteyn Method
    DDye = diff_discrete(x, Dye);
    DDye = filtmooth(DDye, 10);
    DDye = DDye(xifitres);
      iudata.EVDF.Dru.x = xx;
    % EVDF: Formula (3) from [Druyvesteyn1930 Zeitschrift f. Physik A]
      iudata.EVDF.Dru.y = -4*m_e/(A_p*(e0^2)) .* (vp-xx) .* DDye;
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


    % Calculate Output parameters
    iudata.Te = ffit{iRmax}(2);
    iudata.Vp = vp;
    iudata.ni = Iisat_max/(0.61*e0*A_p*sqrt(e0*iudata.Te/m_i));
    iudata.ne = ffit{iRmax}(1) *ne_prefac;
    iudata.neerr = [c(1,2) c(1,1)] *ne_prefac;
    iudata.Teerr = [c(2,1) c(2,2)];
    iudata.Vf= Vf.V;
    % Parameters for plot
    fit.Te = iudata.Te;
    fit.ne = ffit{iRmax}(1);
    % Calculate the Ion Flux: Gamma_i = n_i*<v_th,i>
    % OLD iudata.Fli = iudata.ni*sqrt(2*e0*iudata.Te/m_i);
    % 20130502: Calculate the Ion Flux: Gamma_i = n_i*v_i,sound
    % 20130806: Daisuke: Calculate Ion Flux: Gamma_i = I_sat / (e0*A_p)
    iudata.Fli = 0.61*iudata.ni*sqrt(e0*iudata.Te/m_i);
    
    %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    % XXX  --->>> Check of Fitting
    % if p_chk
    %   figeps(33,25,1,0,0); clf
    %   ne_step = ffit(1) *ne_prefac;
    %   [~, Iefit] = int_discrete(xx, ffun(ffit,xx));
    %   clf;
    %   subplot(5,3,1:3); hold on
    %     plot(Dye, 'r-o');
    %     plot(Dyefit, 'b', 'LineWidth', 2); d=max(Dyefit)-min(Dyefit);
    %     set(gca, 'ylim', [min(Dyefit)-0.80*d max(Dyefit)+0.80*d])
    %   subplot(5,3,4:6); hold on
    %     plot( yefit, 'r-o');
    %     plot( Iefit, 'b', 'LineWidth', 2); d=max(Iefit)-min(Iefit);
    %     set(gca, 'ylim', [min(Iefit)-0.40*d max(Iefit)+0.40*d])
    %   ix = (0:numel(Rsqr1)-1)+xindfit(end-2*N2+1);
    %   subplot(5,3, 7); plot(ix, Rsqr1, 'k-');
    %   subplot(5,3,10); plot(ix, Rstd1, 'b-');
    %   subplot(5,3,13); plot(ix, Rrat1, 'm-');
    %   subplot(5,3, 8); plot(ix, Rsqr2, 'k-');
    %   subplot(5,3,11); plot(ix, Rstd2, 'b-');
    %   subplot(5,3,14);plot(ix,  Rrat2, 'm-');
    %   subplot(5,3, 9); plot(ix, RR,   'k-');
    % puttextonplot(gca,[0 1],5,-15,['ne=' sprintf('%.5e',ne_step)],0,8,'k');
    % input('  >>> Press Enter to Continue <<<')
    % end
    % XXX  <<<--- Check of Fitting
    %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF CONTINUES BELOW AT "NO-OUTPUT" %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXX  <<<--- Check fit of electron current between Vf and Vp
if p_chk
  % Plot fits and raw data
  PlasmaFitIUDemidov_plot
  
   %---------------------------------- Check if manual fitting is activated
   iudata = PlasmaFitIUDemidov_manual(x,xind,y,ye,Dye, ...
                                      ne_prefac,Iisat_max,m_i,iudata);
   return
end
% XXX  <<<--- Check ion saturation current
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manual fitting OR NO OUTPUT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else   % Found at least one lower boundary larger than a upper boundary.
   disp('No automatic fit has been found.')
   disp(' -> Lower and upper fit boundaries are not logic.')

   %---------------------------------- Check if manual fitting is activated
   iudata = PlasmaFitIUDemidov_manual(x,xind,y,ye,Dye, ...
                                      ne_prefac,Iisat_max,m_i,iudata);
   return
  end
  
else  % Found at least one imaginary number, i.e. a fit makes no sense.
  disp('No automatic fit has been found.')
  disp(' -> Boundary values of ne, Te and Vp are not real numbers.')
  PlasmaFitIUDemidov_nooutput

   %---------------------------------- Check if manual fitting is activated
   iudata = PlasmaFitIUDemidov_manual(x,xind,y,ye,Dye, ...
                                      ne_prefac,Iisat_max,m_i,iudata);
   return

end

end