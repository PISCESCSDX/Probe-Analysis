% LOAD RAW DATA
  disp('... loading and save probe data to eva_raw.mat');
  b=dir('B2*');
  [B1 tt] = readmdf(b(1).name);
  [B2 tt] = readmdf(b(2).name);
  [B3 tt] = readmdf(b(3).name);
  clear a b;

% SAVE DATA TO VARIABLES + Take Care of amplifier settings
%==========================================================================
% EMISSIVE PROBE
  VAM502=10; VISO=1/10; V=VAM502*VISO;
  ch=1; EPAC = [B1(:,ch) B2(:,ch) B3(:,ch)]./V;
  ch=2; EPDC = [B1(:,ch) B2(:,ch) B3(:,ch)]./V;
% INTERFEROMETER
  ch=4; INTF = [B1(:,ch) B2(:,ch) B3(:,ch)];
% EXCITER CURRENT MONITOR % 0.1V / A
  V=0.1;
  ch=3; EXCM = [B1(:,ch) B2(:,ch) B3(:,ch)]./V;
clear B1 B2 B3;

save eva_raw.mat tt EPAC EPDC INTF EXCM
disp('continue with <<eva2eva>>');