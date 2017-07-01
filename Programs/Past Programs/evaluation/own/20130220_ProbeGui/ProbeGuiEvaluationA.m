%==========================================================================
% Parameters: Probe, Probe Circuit, DAQ & Plasma
%--------------------------------------------------------------------------
% Plasma parameters
%--------------------------------------------------------------------------
% Gas Type
natconst;
str = get(handles.listbox1, 'String');
num = get(handles.listbox1, 'Value');
elementstring = cell2mat(str(num));
switch elementstring
 case 'Argon'
  plasma.m_ion = m_Ar./u0;
 case 'Deuterium'
  plasma.m_ion = m_D ./u0;
 case 'Helium'
  plasma.m_ion = m_He./u0;
 case 'Hydrogen'
  plasma.m_ion = m_H ./u0;
end

% Magnetic field on the axis
plasma.B = 1e-3 * str2double( get(handles.EditB,'String') );
%--------------------------------------------------------------------------  
% Probe Paramters (Langmuir probe)
%--------------------------------------------------------------------------
% probe shape: cylindrical rod
% Total A (m^2) is given, but we need r, l
  probe.r = 1e-3 * str2double( get(handles.EditProbeR,'String') );
  probe.l = 1e-3 * str2double( get(handles.EditProbeL,'String') );
% assume A is calculated as: A = 2*pi*r*l + pi*(r^2)
  probe.A = 2*pi *probe.r *probe.l + pi*(probe.r^2);
%==========================================================================



disp('Calculating plasma paramters ...')
%>>>=======================================================================
% Evaluation of single characteristics
%--------------------------------------------------------------------------
% This m-file is for radial IV-evaluation using 'PlasmaFitIULangmuir'

% # Input: Probe Theory

% find Pisces-A rawdata file
f = dir('Langmuir_iu*.dat');
lf= length(f);

% Preallocating radial plasma parameters
ne.n = zeros(lf,1);   ne.e = zeros(lf,2);
Te.n = zeros(lf,1);   Te.e = zeros(lf,2);
Vp.n = zeros(lf,1);   Vp.e = zeros(lf,2);
Vf.n = zeros(lf,1);   Vf.e = zeros(lf,2);
RIe  = zeros(lf,1);   RDIe = zeros(lf,1);
% EDF.y = zeros(numel(A(:,1)),length(f));
% EDF.y_raw = zeros(numel(A(:,1)),length(f));
% EDF.x = EDF.y;

%==========================================================================
% Calculate plasma parameters from I-V-characteristic
%--------------------------------------------------------------------------
A = load(handles.data.FileName);
% Assign input variables
  data.voltage = A(:,1);
  data.current = A(:,2);
  data.r = probe.r;
  data.l = probe.l;
  data.m_ion = plasma.m_ion;
  data.B = plasma.B;
% Show plots of the fits
  chk = 1;
% Set the start parameters and the fit boundaries
  fixpar.Te    = [];
  fixpar.Telim = [0.5 50];
  fixpar.ne    = [];
  fixpar.nelim = [];
  fixpar.Vp    = [];
% Use whole characteristic (up to V_p)
%  Vint = [-Inf, NaN];
  Vint = [-Inf, +Inf];
% # Input: Smooth over 'vol_avg' Volts
  vol_avg = 1;
% Ion saturation current range (length) [in Volt]
  p_IsRg = 2;
% Switch between fitting optimal I_e or dI_e/dV (only for Demidov)
  p_fitgoal = 1; % Best fit for (1) I_e (2) dI_e/dV (3) Sum of both
% Stepsize of fit range scan
  p_scanstepsize = str2double( get(handles.EditScanStepsize,'String') );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch handles.data.ProbeTheory
 case 'Langmuir'
  p_ni = 'fitIeGetNeNi';   % other possibility 'fitIGetNe'
  iudata = PlasmaFitIULangmuir_v01(data, chk, fixpar, Vint, vol_avg, ...
           p_IsRg, p_scanstepsize, p_ni);
         
 case 'Demidov'
  iudata = ProbeGui_PlasmaFitIUDemidov_v11(data, chk, fixpar, Vint, vol_avg, ...
           p_IsRg, p_fitgoal, p_scanstepsize, handles);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================

return

%==========================================================================
% Store variables to radial positions
%--------------------------------------------------------------------------
if isempty(iudata)
  ne.n(i,1) = NaN; ne.e(i,:) = NaN;
  ni.n(i,1) = NaN;
  Te.n(i,1) = NaN; Te.e(i,:) = NaN;
  Vp.n(i,1) = NaN; Vp.e(i,:) = NaN;
  Vf.n(i,1) = NaN; Vf.e(i,:) = NaN;
  Fli.n(i,1)= NaN;
  RIe(i,1)  = NaN;  RDIe(i,1) = NaN;
else
  ne.n(i,1) = iudata.ne; ne.e(i,:) = iudata.neerr;
  ni.n(i,1) = iudata.ni;
  Te.n(i,1) = iudata.Te; Te.e(i,:) = iudata.Teerr;
  Vp.n(i,1) = iudata.Vp; Vp.e(i,:) = iudata.Vperr;
  Vf.n(i,1) = iudata.Vf;
  Fli.n(i,1)= iudata.Fli;
  RIe(i,1)  = iudata.RsqrIe;
  RDIe(i,1) = iudata.RsqrDIe;
end
%==========================================================================


rvec = probe.rvec;

switch data.ProbeTheory
  
  case 1 % Langmuir
    % Save radial profiles
        str_p_fitgoal = 'RsqrFitLev';
  
  case 2 % Demidov
    % Save radial profiles
    switch p_fitgoal
      case 1
        str_p_fitgoal = 'RIe';
      case 2
        str_p_fitgoal = 'RDIe';
      case 3
        str_p_fitgoal = 'Sum';
    end
end

savefn = [data.ProbeTheory '_' str_p_fitgoal '.mat'];
save(savefn,'rvec','ne','ni','Vf','Vp','Te','Fli','RIe','RDIe')