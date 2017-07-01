disp('... evaluate data of ev1raw.mat - save to ev2eva.mat');

load ev1raw.mat
al = length(T1AC);

% EVALUATION of PROBE DATA
%==========================================================================

% VECTOR FOR PSDMEAN:
  mpara = [5 4 2 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMISSIVE PROBE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  EP.AC = EPAC;
  sAC = size(EPAC);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CALCULATE f-SPECTRUM
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i=1:size(EPAC,2)
    % Spectrum via Use of PWELCH
    [mPxx{i} EP.f] = psdmean(EP.AC{i}(ipuls),fsample,16,0.5,16, mpara);
  end;
  EP.amp = cell2mat(mPxx);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSPORT PROBE: TRANSPORT, PHASE-SHIFTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('... evaluate transport and phase shifts (phi,n)');
for i=1:al
  filebase{i} = [num_prezero(i,al) '.eps'];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ATTENTION:
  % make sure, that both floating signals have the same standard deviation
  % otherwise: the electric field is completely senseless calculated
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ind = ipuls;
  % # FILL IN INDEX
  taux = tt(ind);
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
  E = - (T1-T2) / Tdist;    % TAKE CARE CONVENTION!!! of GAMMA DIRECTION
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CALCULATE TRANSPORT (y entspricht r, x entspricht phi)
  % GAMMA DIRECTION:   GAMMA positive ==> IN (against CONVENTION)
  % LOOK FROM SOURCE TO PUMP    GAMMA-DIR: along y! ExB>0 => IN
  %                                                 ExB<0 => OUT
  % y^  B:(x)
  %  |              y>0: INSIDE
  %  |  V2   V1
  %  +----------->
  %     |Tdist|  x
  %                 y<0: OUTSIDE
  % GAMMA:
  Gamma.sig{i} = E .* Tn;
  % IN THIS COORDINATE SYSTEM: E IS CALCULATED: E = - (V1-V2)/Tdist
  % IF E is positive: ExB shows inwards (look at B-direction).
  % But DUE to the CONVENTION of GAMMAS DIRECTION then Gamma has to be
  % negative.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % APPLY CONVENTION !!!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Gamma.sig{i} = - Gamma.sig{i};
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
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
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % PLOT TRANSPORT
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  % CALCULATE CROSS POWER SPECTRAL DENSITY - CPSD
  % CROSS PHASES AND CROSS POWER SPECTRUM
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % TAKE THE NEXT LINES TO CALCULATE ONE
%   disp_num(i,al);
%   n=100; olap=0.5;
%   [Amean, phmean, f, phstd] = cpsdmean(T2, Tn, fsample, n, olap, 60e3);
%   Phase.f{i}     = f;
%   Phase.A{i}     = Amean;
%   Phase.ph{i}    = phmean;
%   Phase.phstd{i} = phstd;
  % OR TAKE THESE LINES TO CALCULATE THE MEAN OF BOTH PROBES
  disp_num(i,al);
  n=100; olap=0.5; % n=100 is ok
  [Amean1, phmean1, f, phstd1] = cpsdmean(T1, Tn, fsample, n, olap, 60e3);
  [Amean2, phmean2, f, phstd2] = cpsdmean(T2, Tn, fsample, n, olap, 60e3);
  Phase.f{i}     = f;
  Phase.A{i}     = 0.5*(Amean1+Amean2);
  Phase.ph{i}    = 0.5*(phmean1+phmean2);
  Phase.phstd{i} = 0.5*(phstd1+phstd2);
end

save ev2eva.mat filebase tt Gamma Phase
disp('continue with <<eva3plot>>');

clear all;