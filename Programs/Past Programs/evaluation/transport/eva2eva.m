disp('... evaluate data of ev1raw.mat - save to ev2eva.mat');

load ev1raw.mat
al = length(T1AC);

% EVALUATION of PROBE DATA
%==========================================================================

% vector for the spectrum mean (psdmean)
  mpara = [0 0 0 1];

% % EMISSIVE PROBE - check: ok
% % write time traces to EP
%   EP.AC = EPAC;
%   mEP = mean(EPDC);
%   sAC = size(EPAC);
%     for i=1:sAC(2)
%       M(1:sAC(1), i) = mEP(1,i);
%       EP.DC(:,i) = M(:,i)+EPAC(:,i);      
%     end
% % calculate the spectrum
%   for i=1:size(EPAC,2)
%     % Sampling Frequency
%     Fs=1/(tt(2)-tt(1));
%     % Spectrum via Use of PWELCH
%     [mPxx{i} EP.f]=psdmean(EPAC(:,i),Fs,16,0.5,16, mpara);
%   end;
%   EP.amp = cell2mat(mPxx);

% TRANSPORT PROBE
% CALCULATE the E-field
% E = - grad \phi = -d\phi / dx
for i=1:al
  % make sure, that both floating signals have the same standard deviation
  % otherwise: the electric field is completely senseless calculated
  T1AC{i} = T1AC{i} - mean(T1AC{i});
  T2AC{i} = T2AC{i} - mean(T2AC{i});  
  E{i} = - (T2AC{i}-T1AC{i}) / Tdist;
% CALCULATE TRANSPORT
  Gamma{i} = E{i} .* TnAC{i};
%   i,
%   figeps(10,24,1)
%   hold on
%   subplot(511); plot(T1AC{i}(1:1e3));
%   subplot(512); plot(T2AC{i}(1:1e3));
%   subplot(513); plot(TnAC{i}(1:1e3));
%   subplot(514); plot(E{i}(1:1e3));  
%   subplot(515); plot(Gamma{i}(1:1e3));
% %   print('-depsc2', ['V1-V2-n-E-T_' num2str(i) '.eps']);
%   hold off
%   clf;
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

save ev2eva.mat tt Gamma
disp('continue with <<eva3plot>>');

clear all;