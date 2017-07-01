disp('... evaluate data of eva_raw.mat - save to eva_eva.mat');

load eva_raw.mat

% EVALUATION of PROBE DATA
%==========================================================================

% vector for the spectrum mean (psdmean)
  mpara = [0 0 0 1];

% EMISSIVE PROBE - check: ok
% write time traces to EP
  EP.AC = EPAC;
  EP.DC = EPDC;  
% calculate the spectrum
  for i=1:size(EPAC,2)
    % Sampling Frequency
    Fs=1/(tt(2)-tt(1));
    % Spectrum via Use of PWELCH
    [mPxx{i} EP.f]=psdmean(EPAC(:,i),Fs,16,0.5,16, mpara);
  end;
  EP.amp = cell2mat(mPxx);


% % BDOT PROBE - check: ok
% % calculate the spectrum
% BP1 BP2 BP3
  for i=1:size(BP1,2)
    % Sampling Frequency
    Fs=1/(tt(2)-tt(1));
    % Spectrum via Use of PWELCH
    [Pxx1{i} BP.f]=psdmean(BP1(:,i),Fs,16,0.5,16,mpara);
    [Pxx2{i} BP.f]=psdmean(BP2(:,i),Fs,16,0.5,16,mpara);
    [Pxx3{i} BP.f]=psdmean(BP3(:,i),Fs,16,0.5,16,mpara);
  end;
  BP.amp1 = cell2mat(Pxx1);
  BP.amp2 = cell2mat(Pxx2);
  BP.amp3 = cell2mat(Pxx3);


% INTERFEROMETER - check: no
% Calculate the electron density-profile I.r  I.n
for i=1:3
  [I.ndl{i} I.r{i} I.n{i}] = interf160(mean(INTF(:,i)));
end;


% TRANSPORT PROBE - check: no - gucken, ob wie die richtung des E-feldes
% ist
% CALCULATION of TRANSPORT GammaTilde
  % distance of floating potential probes
    dprobe = 0.011;
  % r-postion of tripel-rpobe
    rpos = 0.05;
  % E-field
    efield = (T2AC-T1AC)./dprobe;
  % calculate the density (from interferometer)
    % find the rpos in Interferometer-r-vector
    b=find( I.r{i} >= rpos );
    % get the appropriate densities
    in = cell2mat(I.n);
    dens = in(b(1),:);
  % calculate the timerow of the transport
    for i=1:size(efield,2)
      TP.GammaTilde(:,i)=efield(:,i).*dens(i);
    end;
% CALCULATION of TRANSPORT SPECTRUM
  for i=1:size( TP.GammaTilde ,2)
    % Sampling Frequency
    Fs=1/(tt(2)-tt(1));
    % Spectrum via Use of PWELCH
    [Pxx{i} freq]=psdmean( TP.GammaTilde(:,i) ,Fs,16,0.5,16,mpara);
  end;
  TP.f = freq;
  TP.amp = cell2mat(Pxx);


% EXCITER-CURRENTS - check: no
% calculate the CURRENTS
  Rshunt = 0.68; % resistance of the shunts /Ohm
  EX.I = EXC./Rshunt;
% calculate the amplitudes and phases of the Exciter-Currents
  for i=1:size(EXC,1)
    for j=1:size(EXC,3)
      [frq amp pha] = su_fft(tt, EXC(i,:,j)');
    % find the maximum peak and frequency
      [maxamp indi] = max(amp);
      EX.EXCamp(i,j)=maxamp;
      EX.EXCfrq(i,j)=frq(indi);
      EX.EXCpha(i,j)=pha(indi);
    end;
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

save eva_eva.mat tt EX EP I BP TP
disp('continue with <<eva3plot>>');

clear all;