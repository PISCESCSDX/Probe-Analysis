disp('... evaluate data of ev1raw.mat - save to ev2eva.mat');

load ev1raw.mat
al = length(T1AC);

% EVALUATION of PROBE DATA
%==========================================================================

% VECTOR FOR PSDMEAN:
  mpara = [5 4 2 1];

% EMISSIVE PROBE - check: ok
% write time traces to EP
  EP.AC = EPAC;
  sAC = size(EPAC);
% calculate the spectrum
  for i=1:size(EPAC,2)
    % Spectrum via Use of PWELCH
    [mPxx{i} EP.f]=psdmean(EPAC(:,i),fsample,16,0.5,16, mpara);
  end;
  EP.amp = cell2mat(mPxx);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSPORT PROBE: TRANSPORT, PHASE-SHIFTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('... evaluate transport and phase shifts (n,phi)');
for i=1:al
  filebase{i} = [num_prezero(i,al) '.eps'];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ATTENTION:
  % make sure, that both floating signals have the same standard deviation
  % otherwise: the electric field is completely senseless calculated
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ind = (4e4:10e4);
  T1 = T1AC{i}(ind);
  T2 = T2AC{i}(ind);
  Tn = TnAC{i}(ind);
  T1 = T1 - mean(T1);  
  T2 = T2 - mean(T2);
  Tn = Tn - mean(Tn);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CALCULATE the E-field
  % E = - grad \phi = -d\phi / dx
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  E = - (T2-T1) / Tdist;
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CALCULATE TRANSPORT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Gamma.sig{i} = E .* Tn;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CALCULATE TRANSPORT - SPECTRUM
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [Gamma.fft{i} Gamma.f] = psdmean(Gamma.sig{i},fsample,15,0.5,16, mpara);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PDF GAMMA
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [pdf_x{i} pdf_y{i}] = pdf(Gamma.sig{i}, 100, 0);
  S=int_discrete(pdf_x{i}, pdf_y{i});
  Gamma.pdfx{i} = pdf_x{i};
  Gamma.pdfy{i} = pdf_y{i}/S;
  %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % PLOT TRANSPORT
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   i, figeps(10,24,1)
%   hold on
%   subplot(511); plot(T1(1:1e3));
%   subplot(512); plot(T2(1:1e3));
%   subplot(513); plot(Tn(1:1e3));
%   subplot(514); plot(E(1:1e3));  
%   subplot(515); plot(Gamma{i}(1:1e3));
%   print('-depsc2', ['V1-V2-n-E-T_' num2str(i) '.eps']);
%   hold off
%   clf;
%   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PHASE BUSINESS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp_num(i,al);
  Vf1 = T1;
  Is  = Tn;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CALCULATE CROSS POWER SPECTRAL DENSITY
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  win = [];   % if empty: 8 windows default
  ovlp = [];  % if empty: 50% overlap
  nfft = [];  %2^11;  % 
  [Pxy,frq] = cpsd(Vf1, Is, win,  ovlp, nfft, fsample);
  i_f = find(frq<30e3);
  %
  Phase.f{i}  = frq(i_f);
  Phase.A{i}  = abs(Pxy(i_f));
  Phase.ph{i} = angle(Pxy(i_f)) / pi;
  %
  % MAX OF Ampl-SPECTRUM
  cutpercent=2;
  l10per = round(length(Phase.A{i}) *cutpercent/100);
  ind = 1+l10per: length(Phase.A{i}) -l10per;
  
  frqcut  = Phase.f{i}(ind);
  ampcut  = Phase.A{i}(ind);
  phascut = Phase.ph{i}(ind);
  [spMax i_sM] = max( ampcut );
  Phase.f_Max{i} = frqcut(i_sM);
  Phase.phMax{i} = phascut(i_sM) *180;
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % OLD PHASE EVALUATION VIA CORRELATION FUNCTION CFF
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %     [frq amp pha] = su_fft(tr, Is);
  % % Cross Correlation Function
  %     [cpp_v1i, tau] = ccf(Vf1,  Is, fsample);
  % % % CALCULATE THE MAIN FREQUENCY
  % %     [fr, Iam, Iph] = su_fft(tau', cpp_v1i);
  % %     plot(fr, Iam);
  % % PLOT Correlation functions
  %     dtau = ccfmax1(cpp_v1i, tau, 180);
  %       [frq amp pha] = su_fft(tr, Is);
  %       [mval mi] = max(amp);
  %       T = 1/frq(mi);
  %       ccfpha(i) = dtau/T*360,
  % %     d=5e2;
  % %     ind = (round(size(tau,2)/2)-d:round(size(tau,2)/2)+d);
  % %     fig(10,20,2)
  % %     subplot(311);
  % %       plot(tr(ind), Vf1(ind), 'b');
  % %       title('Vf1');
  % %     subplot(312);
  % %       plot(tr(ind), Is(ind), 'r');
  % %       title('Is');
  % %     subplot(313);
  % %       plot(tau(ind), cpp_v1i(ind), 'g');
  % %       title('cross correlation Vf1-Is');
  % %     clf
end


% % INTERFEROMETER - check: no
% % Calculate the electron density-profile I.r  I.n
% for i=1:3
%   [I.ndl{i} I.r{i} I.n{i}] = interf160(mean(INTF(:,i)));
% end;


% % CURRENT MONITOR - check: ok
% % Compare the measured current with the measured current of the
% % shunt-system! Get amplitude of the currents.
% EX.CM = EXCM;
% for i=1:size(EXCM,2);
%   [frq amp pha] = su_fft(tt, EXCM(:,i));
%   % find the maximum peak and frequency
%   [maxamp indi] = max(amp);
%   EX.CMf(i)=frq(indi);
%   EX.CMamp(i)=maxamp;
% end;

save ev2eva.mat filebase tt Gamma Phase
disp('continue with <<eva3plot>>');

clear all;