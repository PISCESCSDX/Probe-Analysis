% FILL IN ALL DEFINITION VARIABLES (INPUT#)

% INPUT# DISTANCE of TRANSPORT Vfl-probes [m]
  Tdist = 0.01065;
% INPUT# CURRENT MONITOR SENSITIVITY
  VCM = 0.1; % sensitivity [V/A]
% INPUT# EMISSIVE PROBE - AMPLIFIERS
  VAM502=20; VISO=0.10; Vep=VAM502*VISO;
% INPUT# TRANSPORT PROBE - Vfl1
  VAM502=20; VISO=0.05; Vt1=VAM502*VISO;
% INPUT# TRANSPORT PROBE - Vfl2
  VAM502=20; VISO=0.05; Vt2=VAM502*VISO;
% INPUT# TRANSPORT PROBE - Isat AC
  VAM502=10; VISO=0.05; Vtnac=VAM502*VISO;
% INPUT# TRANSPORT PROBE - Isat DC
  VAM502= 5; VISO=0.10; Vtndc=VAM502*VISO;
% INPUT# DEFINITION OF CHANNELS BOARD 1
% INPUT# DEFINITION OF CHANNELS BOARD 2
  ch_tnac= 1;
  ch_tndc= 2;  
  ch_tv1 = 3;
  ch_tv2 = 4;
  ch_cm  = 5;
  ch_int = 6;
  ch_ep  = 7;

% LOAD RAW DATA
  disp('... loading and save probe data to eva_raw.mat');
  b=dir('B2*'); bl = length(b);

for i=1:bl
  [B tt] = readmdf(b(i).name);

% SAVE DATA TO VARIABLES + Take Care of amplifier settings
%==========================================================================
% EMISSIVE PROBE:
  EPAC{i} = B(:,ch_ep)./Vep;

% TRANSPORT PROBE Vf1
  T1AC{i} = B(:,ch_tv1)./Vt1;
    T1AC{i} = T1AC{i}-mean(T1AC{i});
    T1AC{i} = T1AC{i}./std(T1AC{i});
% TRANSPORT PROBE Vfl2
  T2AC{i} = B(:,ch_tv2)./Vt2;
    T2AC{i} = T2AC{i}-mean(T2AC{i});
    T2AC{i} = T2AC{i}./std(T2AC{i});
% TRANSPORT PROBE - Isat AC
  TnAC{i} = B(:,ch_tnac)./Vtnac;
    TnAC{i} = TnAC{i}-mean(TnAC{i});
% TRANSPORT PROBE - Isat DC
  TnDC{i} = B(:,ch_tndc)./Vtndc;

% INTERFEROMETER
  INTF{i} = B(:,ch_int);

% EXCITER CURRENT MONITOR
  EXCM{i} = B(:,ch_cm) ./VCM;

% figeps(10,24,1)
%   subplot(611); plot(tt(1:1e3)*1e3,T1AC{i}(1:1e3), 'r');
%   subplot(612); plot(tt(1:1e3)*1e3,T2AC{i}(1:1e3), 'b');
%   subplot(613); plot(tt(1:1e3)*1e3,TnAC{i}(1:1e3), 'k');
%   subplot(614); plot(tt(1:1e3)*1e3,EPAC{i}(1:1e3), 'g');
%   subplot(615); plot(tt(1:1e3)*1e3,INTF{i}(1:1e3), 'm');
%   subplot(616); plot(tt(1:1e3)*1e3,EXCM{i}(1:1e3));
%   clf;
clear A B;

end
clear a b;

save ev1raw.mat tt EPAC Tdist T1AC T2AC TnAC TnDC INTF EXCM
disp('continue with <<eva2eva>>');