disp('... evaluate data of eva_raw.mat - save to eva_eva.mat');

load eva_raw.mat

% EVALUATION of PROBE DATA
%==========================================================================

% vector for the spectrum mean (psdmean)
  mpara = [0 0 0 1];

% EMISSIVE PROBE - check: ok
% write time traces to EP
  EP.AC = EPAC;
  mEP = mean(EPDC);
  sAC = size(EPAC);
    for i=1:sAC(2)
      M(1:sAC(1), i) = mEP(1,i);
      EP.DC(:,i) = M(:,i)+EPAC(:,i);      
    end
% calculate the spectrum
  for i=1:size(EPAC,2)
    % Sampling Frequency
    Fs=1/(tt(2)-tt(1));
    % Spectrum via Use of PWELCH
    [mPxx{i} EP.f]=psdmean(EPAC(:,i),Fs,16,0.5,16, mpara);
  end;
  EP.amp = cell2mat(mPxx);


% INTERFEROMETER - check: no
% Calculate the electron density-profile I.r  I.n
for i=1:3
  [I.ndl{i} I.r{i} I.n{i}] = interf160(mean(INTF(:,i)));
end;


% CURRENT MONITOR - check: ok
% Compare the measured current with the measured current of the
% shunt-system! Get amplitude of the currents.
EX.CM = EXCM;
for i=1:size(EXCM,2);
  [frq amp pha] = su_fft(tt, EXCM(:,i));
  % find the maximum peak and frequency
  [maxamp indi] = max(amp);
  EX.CMf(i)=frq(indi);
  EX.CMamp(i)=maxamp;
end;

save eva_eva.mat tt EX EP I
disp('continue with <<eva3plot>>');

clear all;