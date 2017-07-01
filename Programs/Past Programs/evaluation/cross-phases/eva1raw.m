% DEFINITIONS for 20080215
% FILL IN ALL DEFINITION VARIABLES (INPUT#)
clear all

% INFO VISO Sony A6907:
% Wenn Ain 1V übersteigt muss VISO so gewählt werden, dass Aout 1V nicht
% übersteigt!
% Bsp: VISO=0.10 -> Ain=1V; A6907: 1V/div; Aout=0.10V
% Bsp: VISO=0.05 -> Ain=1V; A6907: 2V/div; Aout=0.05V
% Bsp: VISO=1.00 -> Ain=1V; A6907: 100mV/div; Aout=1V
% Bsp: VISO=1.00 -> Ain=1V; A6907: 100mV/div; Aout=1V
% VISO=1/(Anzeige*10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('meas_settings.mat')
% INPUT# DISTANCE of TRANSPORT Vfl-probes [m]
  Tdist = input('distance Vf12-probes [m]             (0.0092): ');
  if isempty(Tdist); Tdist = 0.0092; end
% INPUT# EMISSIVE PROBE - AMPLIFIERS
  Vep =   input('isolation amplifier - emissive probe    (0.5): ');
  if isempty(Vep); Vep = 0.5; end
  VISO=1/(10*Vep); Vep=VISO;
% INPUT# TRANSPORT PROBE - Vfl1
  Vt1 =   input('isolation amplifier - Gamma probe Vf1   (0.5): ');
  if isempty(Vt1); Vt1 = 0.5; end
  VISO=1/(10*Vt1); Vt1=VISO;
% INPUT# TRANSPORT PROBE - Vfl2
  Vt2 =   input('isolation amplifier - Gamma probe Vf2   (0.5): ');
  if isempty(Vt2); Vt2 = 0.5; end
  VISO=1/(10*Vt2); Vt2=VISO;
% INPUT# TRANSPORT PROBE - Isat AC
  Vti =   input('isolation amplifier - Gamma probe Isat  (0.5): ');
  if isempty(Vti); Vti = 0.5; end
  VISO=1/(10*Vti); Vti=VISO;
% INPUT# EXCITER VOLTAGE_OUT
  VEX =      input('exciter voltage out reducer            (1/10): ');
  if isempty(VEX); VEX = 1/10; end
% INPUT# PULSE ON AT WHICH MEASUREMENT? (default: 2 co )
index_pulse =input('measurement no with pulse ("-1": nopulse)    : ');
  if isempty(index_pulse); 
    index_pulse = -1;
  disp('No pulse will be detected, the whole time row will be evaluated!');
  end
save meas_settings.mat Tdist Vep Vt1 Vt2 Vti VEX index_pulse
else
  load('meas_settings.mat');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUT# DEFINITION OF CHANNELS BOARD 1
% ch1-8 Current Monitor Box - El Exciter
% INPUT# DEFINITION OF CHANNELS BOARD 2
  ch_lp1 = 1;
  ch_lp2 = 2;
  ch_tp1 = 3;
  ch_tpi = 4;
  ch_tp2 = 5;
  ch_ep  = 6;
  ch_int = 7;
  ch_vex = 8;

% LOAD RAW DATA
  disp('... loading and save probe data to eva_raw.mat');
  a=dir('B1*'); al = length(a);
  b=dir('B2*'); bl = length(b);

for i=1:al
  [A tt] = readmdf(a(i).name);
  [B tB] = readmdf(b(i).name);

% Sampling Frequency
  fsample = 1/(tt(2)-tt(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PULSE START & END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if i==index_pulse
  [t0, t1, i0, i1] = pulse_limits(tB, B(:,8));
  ipuls = i0:i1;
end

% SAVE DATA TO VARIABLES + Take Care of amplifier settings
%==========================================================================
% EMISSIVE PROBE:
  EPAC{i} = B(:,ch_ep)./Vep;

% TRANSPORT PROBE Vf1
  T1AC{i} = B(:,ch_tp1)./Vt1;
%    T1AC{i} = T1AC{i}-mean(T1AC{i});
%    T1AC{i} = T1AC{i}./std(T1AC{i});
% TRANSPORT PROBE Vfl2
  T2AC{i} = B(:,ch_tp2)./Vt2;
%    T2AC{i} = T2AC{i}-mean(T2AC{i});
%    T2AC{i} = T2AC{i}./std(T2AC{i});
% TRANSPORT PROBE - Isat AC
  TnAC{i} = B(:,ch_tpi)./Vti;
%    TnAC{i} = TnAC{i}-mean(TnAC{i});

% INTERFEROMETER
  INTF{i} = B(:,ch_int);

% EXCITER VOLTAGE Electrode #1
  EXVE{i} = B(:,ch_vex) ./VEX;

figeps(10,24,1)
  subplot(611); plot(tt(1:1e3)*1e3,T1AC{i}(1:1e3), 'r');
  subplot(612); plot(tt(1:1e3)*1e3,T2AC{i}(1:1e3), 'b');
  subplot(613); plot(tt(1:1e3)*1e3,TnAC{i}(1:1e3), 'k');
  subplot(614); plot(tt(1:1e3)*1e3,EPAC{i}(1:1e3), 'g');
  subplot(615); plot(tt(1:1e3)*1e3,INTF{i}(1:1e3), 'm');
%  subplot(616); plot(tt(1:1e3)*1e3,EXCM{i}(1:1e3));
  clf;
clear A B;

end
clear a b;

if isempty(ipuls)
  disp('NO PULSE DETECTED! TAKE WHOLE TIME TRACE!!!');
  ipuls = 1:lengtj(B);
end

save ev1raw.mat tt fsample EPAC Tdist T1AC T2AC TnAC INTF EXVE ipuls
disp('continue with <<eva2eva>>')
clear all