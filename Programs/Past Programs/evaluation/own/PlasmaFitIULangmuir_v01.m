function iudata = PlasmaFitIULangmuir_v01(data, p_chk, p_fix, p_Vint, ...
  p_Vavg, p_IsRg, p_scanstepsize, p_fitmethod)
%==========================================================================
% function iudata = PlasmaFitIULangmuir_v01(data, p_chk, p_fix, p_Vint, ...
%   p_Vavg, p_IsRg, p_scanstepsize, p_fitmethod)
%
% Last Change: 2013-02-15 17:08 C. Brandt, San Diego
% This is alpha version: v11.1
% This version tries to be more stable and flexible.
% For Details look into the file: 2012-04-19_IUFit.pdf
%--------------------------------------------------------------------------
% PlasmaFitIUDemidov_v11 calculates the plasma parameters n_e (electron 
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
% p_fitmethod: 'fitIeGetNeNi': calculate ne and ni from Iisat and fit to Ie
%       'fitIGetNe'   : calculate ne=ni from fit to total I
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


% 1. Cut Characteristic (both sides 5%)
%--------------------------------------------------------------------------
% Extract the characteristic and the probe parameters
% Cut margin: 2% both sides
L5 = round(0.02*length(data.voltage));
ind = L5:length(data.voltage)-L5;
data.voltage = data.voltage(ind);
data.current = data.current(ind);

x = data.voltage;
xind = 1:numel(x);
y = data.current;
  ll = numel(y);
  dx = (x(end)-x(1)) / numel(x);


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




%--------------------------------------------------------------------------
% 4.2.1 Determine Fit Start Parameters
%--------------------------------------------------------------------------
% Start parameter Te_V at Plasma Potential: Te_V = -I_e(phi) / I_e'(phi)
% Calculate Derivative again used for fitting:
  Dye   = diff_discrete(x, ye);
  Dyefit= diff_discrete(xfit, yefit);
  Dyesm = filtmooth(Dyefit, par_ysm);

% Start parameter Te_V at Plasma Potential Phi: Te_V = -I_e(phi)/I_e'(phi)
if isempty(p_fix.Te)
  [~, iVp] = min(Dyesm);
  Te_start = yefit(iVp) / Dyesm(iVp);
else
  Te_start = p_fix.Te;
end

% Start parameter for n_e (use electron saturation current with Te_start)
% Define pre factor for fit parameter 'a' of ne (ne = a / ne_prefac)
ne_prefac = (-e0 *A_p *sqrt(e0/(2*pi*m_e)));
if isempty(p_fix.ne)
  ne_start = ye(iVp) / (ne_prefac*sqrt(Te_start));
else
  ne_start = p_fix.ne;
end


% Define fit start parameters
fit_start = [ne_start Te_start];


%--------------------------------------------------------------------------
% 4.2.2 Determine Fit Boundaries
%--------------------------------------------------------------------------
% Set fit limits [n_e, Te_V, Vp]
if isempty(p_fix.nelim)
  n0 = 1e-6*ne_start;
  n1 = 1e+6*ne_start;
else
  n0 = p_fix.nelim(1);
  n1 = p_fix.nelim(2);
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
% Between the floating potential and the plasma potential an exponential 
% increase of the electron current is assumed:
% e0: elementary charge
% A_p: effective probe surface
% V: probe potential
% phi: plasma potential
% Te_V: electron temperature in Volt (not in eV!)
% 
% The current to the probe is the sum of electron and ion current, whereas
% the ion current can be approximately seen as the saturation current.
% I(V) = I_e(V) + I_isat
%
% I_e(V)  = I_esat * exp(-(phi-V)/Te_V)
% I_esat  = -e0 *n_e *A_p *<v+>              % <v+> avg velocity to surface
%   -->   = -e0 *n_e *A_p *sqrt(e0*Te_V/(2*pi*m_e))
%
% I_isat  = 0.61 *e0 *n_e *A_p *sqrt( e0 *Te_V /m_i)
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% % Fit Method: 'fitIeGetNeNi' n_e from fit to I_e and n_i from I_isat
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% I_e(V) = -e0*n_e*A_p*sqrt(e0*Te_V/(2*pi*m_e))*exp(-(Phi-V)/Te_V)
%
% Rearrange for fit function f(x):
% I_e(V) = -e0*A_p*sqrt(e0/(2*pi*m_e))*n_e*sqrt(Te_V)*exp(-(Phi-V)/Te_V)
%
% Fit function with 3 fit parameters: f(x) = a*sqrt(b)*exp(-(c-x)/b)
%
%      ne(a) = a / (-e0 *A_p *sqrt(e0/(2*pi*m_e)))
%           ne_prefac = (-e0 *A_p *sqrt(e0/(2*pi*m_e)))
%    Te_V(b) = b
%        Phi = c
%
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% TODO *** (B) n_e=n_i of ion saturation current
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% I_e(V) = -e0*n_e*A_p*sqrt(e0*Te_V/(2*pi*m_e))*exp(-(Phi-V)/Te_V)
% From the ion saturation current with n_e=n_i
% and n_i = I_i,sat/(0.61*e0*A_p*sqrt(e0*Te_V/m_i))
% I_e(V) = -(1/0.61)*I_i,sat*sqrt(m_i/(2*pi*m_e))*exp(-(Phi-V)/Te_V)
% Rearrange for fit function f(x):
% Fit function with 2 fit parameters: a:Phi, b:Te_V
%     f(x) = v1*exp(-(a-x)/b)

%==========================================================================
% Fit Method: 'fitIeGetNeNi'
%--------------------------------------------------------------------------

%==========================================================================
% 2013-02-15 Fit Method: 'fitIGetNe'
%--------------------------------------------------------------------------
% Measured Current: I = Ie(V) + Ii(V)
% Compared to the large electron current between V_f and V_p the ion 
% current can be simplified by just using the relatively small ion
% saturation current.
% This means: I = Ie(V) + Iisat      with
%   Ie(V) = Iesat*exp(-e*(Phi-V)/Te_eV)
%   Iesat = -e*ne*A*sqrt(Te_eV/(2*pi*m_e));
%   Iisat = 0.61*e*n_e*A*sqrt(Te_eV/m_i)


% Set the fit options for LSQCURVEFIT
options = optimset('Display','off','TolFun',1e-8);
% In case upper boundaries are smaller than lower boundaries or numbers are
% not real, do not fit, and output NaN.
if isreal([n0 n1 T0 T1])
  if n1>=n0 && T1>=T0
    % Find the best fit depending on the fit range and thus on V_plasma
    % Half-Length of Range variation (full range = 2*N2, see below):
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    % Fit Loop: Length of fit range is scanned
    %======================================================================
    % Deviations between rawdata and fits to dIe/dV and Ie are calculated
    fitscanrange = 2/3;  % # Input: From right to left (>0 && <1)
    jvec = 1:p_scanstepsize:round(fitscanrange *length(xindfit));
    % Initialize Statistic Variables
    Rsqr1 = zeros(length(jvec),1); Rstd1 = Rsqr1;
    Vploop = Rsqr1;    
    
          %%%%%%%%%%%%%%%%
          clf;
          figeps(10,18,1);
          subplot(2,1,1);
          hold on
          plot(xfit,yefit,'bo');
          hp = plot(0,0,'r');
          hold off
          mkplotnice('V (V)', 'I (A)', 12, '-20', '-30');
          subplot(2,1,2);
          hold on
          hp2 = plot(1,1);
          title('Rsqr (b:R(I_e))')
          hold off
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
      % Switch between fit procedure (A) and (B): Fit Limits & Fit Fct.
      switch p_fitmethod
        case 'fitIeGetNeNi'
          lb = [n1*ne_prefac T0];   % lower boundary n1 since Ie is neg.
          ub = [n0*ne_prefac T1];   % upper boundary n0 since Ie is neg.
          ffun = @(c,x) c(1)*sqrt(c(2)).*exp(-(vp-x)./c(2));
          % Fit whole electron current up to the plasma potential
          xx = x(iloopfit); yy = ye(iloopfit);
        case 'fitIGetNe'
          ub = [n0 T1];
          lb = [n1 T0];
          ffun = @(c,x) c(1)*sqrt(c(2)).*exp(-(vp-x)./c(2));
          % Fit whole electron current up to the plasma potential
          xx = x(iloopfit); yy = y(iloopfit);
      end
      
      % XXX information: ysm is the total current
      
      %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      ffit{ctr_j} = lsqcurvefit(ffun, fit_start, xx, yy, lb, ub, options);
      %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      
      % Chi Square Deviation Rsqr2 and Standard Deviation of Ie
      % Calculate Fit to Electron Current from fit to Derivative:
      Ifit = ffun(ffit{ctr_j},xx);
      switch p_fitmethod
       case 'fitIeGetNeNi'
        yy2 = ye(iloopfit);
       case 'fitIGetNe'
        yy2 = y(iloopfit);
      end
      sef = sum( (yy2-Ifit).^2 );
      sdm = sum( (yy2-mean(yy2)).^2 );
      Rsqr1(ctr_j) = 1 - sef/sdm;     % The larger the better
      Rstd1(ctr_j) = std(yy2-Ifit);   % The smaller the better
        
          %%%%%%%%%%%%%%%%
          subplot(2,1,1);
          hold on
          set(hp, 'visible', 'off')
          hp = plot(xx,Ifit,'r','LineWidth',2); drawnow
          hold off

          subplot(2,1,2); 
          set(hp2,'visible', 'off')
          hold on
          %hp2 = plot(jvec,Rsqr1' ./ ,'b');
          %Rmeasure = diff_discrete(jvec',Rsqr1) ./ Rsqr1;
          hp2 = plot(jvec, Rsqr1,'b');
          hold off
          drawnow
          %%%%%%%%%%%%%%%%
        
    end
    hold off

   
%     %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%     % Take Parameters of Optimal Fit
%     %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%     % Combined Parameter: Chi of dIe/dV and of Ie
%     % Find fit with the lowest deviation (best plasma potential)
%     choice_R = Rmeasure;
%     [ymax, iRmax] = min(choice_R);
%       k = choice_R > 0.95*ymax;
%         hev = ctr_jvec(k);
%         if iRmax==1 || iRmax==numel(xindfit)
%           disp(['PlasmaFitIULangmuir Error: Optimal Plasma Potential ' ...
%             'outside of scanned range.'])
%           iudata = [];
%           return
%         else
%           kL = hev>iRmax; kL = hev(kL);
%           kS = hev<iRmax; kS = hev(kS);
%           if ~isempty(kL) && ~isempty(kS)
%             kS = kS(1);
%             kL = kL(end);
%           else
%      disp('No fit: Did not find plasma potential, and hence no opt. fit.')
%      disp('Try:    Maybe increase the fit scan range!')
%             iudata.ni  = NaN;
%             iudata.ne  = NaN; iudata.neerr = NaN;
%             iudata.Te  = NaN; iudata.Teerr = NaN;
%             iudata.Vp  = NaN; iudata.Vperr = NaN;
%             iudata.Vf  = NaN;
%             iudata.Vf2 = NaN;
%             iudata.Fli = NaN;
%             iudata.RsqrIe  = NaN;
%             iudata.RsqrDIe = NaN;
%             iudata.EDF.y=NaN;
%             return
%           end
%         end


    % Definition of the Optimal Fit: First Fit which reaches Rsqr1 > fitlev
    fitlev = 0.97;
    Rmin = min(Rsqr1); Rmax = max(Rsqr1); Rdiff = Rmax - Rmin;
    Rind = Rsqr1 > Rmin + fitlev * Rdiff;
    ih = 1:length(Rsqr1);
    dump.ind = ih(Rind);
    iRmax = dump.ind(1);
    iudata.Vperr = NaN; % Wdth Rrat1-peak below 95% max
    
    %iudata.Vperr = [Vploop(kL) Vploop(kS)]; % Wdth Rrat1-peak below 95% max
    iudata.RsqrIe  = Rsqr1(iRmax);
    iudata.RsqrDIe = NaN;
    % Calculate the Optimal Fit
      xifitres = xindfit(1:end-jvec(iRmax));
      xx = x(xifitres);
    % Define Plasma Potential  
      vp = x(xifitres(end));
      switch p_fitmethod
        case 'fitIeGetNeNi'
          yy = ye(xifitres);
          ffun = @(c,x) c(1)*sqrt(c(2)).*exp(-(vp-x)./c(2));
        case 'fitIGetNe'
          yy = y(xifitres);
          ffun = @(c,x) c(1)*sqrt(c(2)).*exp(-(vp-x)./c(2));
      end
    % Calculate Optimal Fit
      Ifit = ffun(ffit{iRmax},xx);
    % Extract Confidence Intervals
      a0 = [ffit{iRmax}(1) ffit{iRmax}(2)];
      options = statset('FunValCheck', 'Off');
      [a,r,~,COVB] = nlinfit(xx,yy,ffun,a0,options);
    c = nlparci(a,r,'covar',COVB);

    % Calculate Output parameters (depending on fit method p_fitmethod)
    switch p_fitmethod
     case 'fitIeGetNeNi'
      iudata.Te = ffit{iRmax}(2);
      iudata.Vp = vp;
      iudata.Vperr = [vp vp]; % Actually errors are not calculated
      iudata.ni = Iisat_max/(0.61*e0*A_p*sqrt(e0*iudata.Te/m_i));
      iudata.ne = ffit{iRmax}(1)/(-e0*A_p*sqrt(e0/(2*pi*m_e))); % checked
      iudata.neerr = [c(1,2) c(1,1)] / ne_prefac;
      iudata.Teerr = [c(2,1) c(2,2)];
     case 'fitIGetNe'
      iudata.Te = ffit{iRmax}(2);
      iudata.Vp = vp;
      iudata.Vperr = [vp vp]; % Actually errors are not calculated
      iudata.ni = Iisat_max/(0.61*e0*A_p*sqrt(e0*iudata.Te/m_i));
      iudata.ne = iudata.ni;
      iudata.neerr = [c(1,2) c(1,1)] / ne_prefac;
      iudata.Teerr = [c(2,1) c(2,2)];
    end
    % Store the floating potential
      iudata.Vf = Vf.V;
    % Calculate the Ion Flux: Gamma_i = n_i*<v_th,i>
      iudata.Fli = iudata.ni*sqrt(2*e0*iudata.Te/m_i);
    % Extract Electron Velocity Distribution Function (Druyvesteyn Method)
      DDye = diff_discrete(x, Dye);
      DDye = filtmooth(DDye, 10);
      DDye = DDye(xifitres);
      iudata.EVDF.Dru.x = xx;
    % EVDF: Formula (3) from [Druyvesteyn1930 Zeitschrift f. Physik A]
      iudata.EVDF.Dru.y = -4*m_e/(A_p*(e0^2)) .* (vp-xx) .* DDye;
    % Store Fit method: 
      iudata.fitmethod = p_fitmethod;
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    
    
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

% %-----------------------------------------------------------------
% % 4.4 Calculate the floating potential from fit an linear I_i,sat
% %-----------------------------------------------------------------
% % With smoothing, the potential, where the current vanishes is looked for.
% % Before adding the linear ion saturation current fct. check for problems:
% % Problem Case 1: Iisat has negative values -> make them to NaN
%   Iisat = polyval(PolyIsat, xx);
%   ind = Iisat<0;
%   Iisat(ind) = NaN;
% % Sum fit of electron current and linear fit of ion saturation current:
% yfitei = Iefit + Iisat;
%   il = find(yfitei> 0);
%   ir = find(yfitei<=0);
% % Test whether points around zero crossing exist (may be NaN)
% if numel(ir)<2
%   disp('PlasmaFitIUDemidov Error in finding Vf from the fit:');
%   disp('not enough points at negative current.');
%   iudata = [];
%   return
% end
% if ~isempty(il) && ~isempty(ir)
%   % 2nd order polynomial at zero crossing
%   ind = il(end-2):ir(2);                  % indices at zero crossing
%   pf = polyfit(x(ind), yfitei(ind), 2);   % polynomial
%   % Find y=zero crossings of parabula:
%   hp = pf(2)/pf(1); hq = pf(3)/pf(1);
%   x01 = -hp/2 + sqrt((hp^2)/4 - hq);
%   x02 = -hp/2 - sqrt((hp^2)/4 - hq);
%   % Find closest to I_e+I_i fit function
%   if abs(x01-x(il(end))) < abs(x02-x(il(end)))
%     % x01 is the zero crossing
%     iudata.Vf2 = x01;
%   else
%     % x02 is the zero crossing
%     iudata.Vf2 = x02;
%   end
%   
% else
%   iudata.Vf2 = NaN;
% end


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXX  <<<--- Check fit of electron current between Vf and Vp
if p_chk

figeps(16,9,5,32,100); clf;
hold on
switch p_fitmethod
  case 'fitIeGetNeNi'
    yshow = ye;
  case 'fitIGetNe'
    yshow = y;
end
  plot(xind    ,1e3*yshow   ,'o-','Color', rgb('NavajoWhite'))
  plot(xifitres,1e3*Ifit,'r-','LineWidth', 2)
hold off
d = max(yshow)-min(yshow);
set(gca, 'ylim', 1e3*[min(yshow)-0.02*d max(yshow)+0.02*d])
set(gca, 'xlim', [1 ll])
mkplotnice('probe voltage (V)', 'fit I_e (mA)', 12, '-25', '-40');
puttextonplot(gca, [0 0], 10, 90, ...
                ['R=' sprintf('%4.2f', iudata.RsqrIe)], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 70, ...
                ['T_e=' sprintf('%4.2f', iudata.Te) 'eV'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 50, ...
                ['n_e=' sprintf('%4.2g', iudata.ne) 'm^{-3}'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 30, ...
                ['V_p=' sprintf('%4.2f', iudata.Vp) 'V'], 0, 12, 'k');
puttextonplot(gca, [0 0], 10, 10, ...
                ['V_f=' sprintf('%4.2f', iudata.Vf) 'V'], 0, 12, 'k');

% input('  >>> Press Enter to Continue <<<')
end
% XXX  <<<--- Check ion saturation current
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


%%%%%%%%%%%%%%%
% "NO-OUTPUT" %
%%%%%%%%%%%%%%%
  else   % Found at least one lower boundary larger than a upper boundary.
   disp('No fit: Lower and upper boundaries are not logic.')
    iudata.ni  = NaN;
    iudata.ne  = NaN; iudata.neerr = NaN;
    iudata.Te  = NaN; iudata.Teerr = NaN;
    iudata.Vp  = NaN; iudata.Vperr = NaN;
    iudata.Vf  = NaN;
    iudata.Vf2 = NaN;
    iudata.Fli = NaN;
    iudata.Rsqr= NaN;
    iudata.Rsqr2=NaN;
    return
  end
else  % Found at least one imaginary number, i.e. a fit makes no sense.
  disp('No fit: Boundary values of ne, Te and Vp are not real numbers.')
    iudata.ni  = NaN;
    iudata.ne  = NaN; iudata.neerr = NaN;
    iudata.Te  = NaN; iudata.Teerr = NaN;
    iudata.Vp  = NaN; iudata.Vperr = NaN;
    iudata.Vf  = NaN;
    iudata.Vf2 = NaN;
    iudata.Fli = NaN;
    iudata.Rsqr= NaN;
    iudata.Rsqr2=NaN;
  return
end

end