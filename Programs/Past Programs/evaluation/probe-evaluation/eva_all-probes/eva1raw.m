% LOAD RAW DATA
  disp('... loading and save probe data to eva_raw.mat');
  a=dir('B1*');
  b=dir('B2*');
  c=dir('cou*');
  [A1 tt] = readmdf(a(1).name);
  [A2 tt] = readmdf(a(2).name);
  [A3 tt] = readmdf(a(3).name);
  [B1 tt] = readmdf(b(1).name);
  [B2 tt] = readmdf(b(2).name);
  [B3 tt] = readmdf(b(3).name);
  [C1 tt] = readmdf(c(1).name);
  [C2 tt] = readmdf(c(2).name);
  [C3 tt] = readmdf(c(3).name);
  clear a b c;

% SAVE DATA TO VARIABLES + Take Care of amplifier settings
%==========================================================================
% EMISSIVE PROBE
  VAM502=10;
  EPAC = [B1(:,1) B2(:,1) B3(:,1)]./VAM502;
  EPDC = [B1(:,2) B2(:,2) B3(:,2)]./VAM502;
% TRANSPORT PROBE
  VAM502=10;
  T1AC = [C1(:,60) C2(:,60) C3(:,60)]./VAM502;
  T1DC = [C1(:,60) C2(:,60) C3(:,60)]./VAM502;
  T2AC = [B1(:,3) B2(:,3) B3(:,3)]   ./VAM502;
  T2DC = [B1(:,4) B2(:,4) B3(:,4)]   ./VAM502;
  TnAC = [C1(:,61) C2(:,61) C3(:,61)]./VAM502;
  TnDC = [C1(:,62) C2(:,62) C3(:,62)]./VAM502;
% BDOT PROBE
  BP1 = [B1(:,6) B2(:,6) B3(:,6)];
  BP2 = [B1(:,7) B2(:,7) B3(:,7)];
    BP2 = BP2./10; % 1:10 Spannungsteiler benutzt
  BP3 = [B1(:,8) B2(:,8) B3(:,8)];
% INTERFEROMETER
  INTF = [C1(:,64) C2(:,64) C3(:,64)];
% EXCITER CURRENTS
  VAM502 = 5;
  EXC(1,:,:) = A1;
  EXC(2,:,:) = A2;
  EXC(3,:,:) = A3;
  EXC = EXC./VAM502;
% EXCITER CURRENT MONITOR % 0.1V / A
  EXCM = [C1(:,63) C2(:,63) C3(:,63)].*10;
clear A1 A2 A3 B1 B2 B3 C1 C2 C3;

save eva_raw.mat tt EPAC EPDC T1AC T1DC T2AC T2DC TnAC TnDC BP1 BP2 BP3 INTF EXC EXCM
disp('continue with <<eva2eva>>');