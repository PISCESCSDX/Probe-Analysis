function [arn] = arnoldtongue_manually(numlist, fdw, PulseFile, kfon)
%20080226-Di-13:47 Brandt
%function [le ri] = arnoldtongue_manually(numlist, fdw)
% IN: numlist: no. of cou0007.MDF files (e.g. = 7)
%     fdw: approx. self sustained oscillator frequency [Hz]
%          (determines delta_f in kf-spectrogram)
%     PulseFile: No (of numlist) (for pulse detection and current
%     determination)
%           -1: no PulseFile, no PC-Frisch-files
%OUT: Arnold tongue
% EX: A = arnoldtongue_manually([1:26], 6000, 24);

if nargin < 4; kfon=1; end
if nargin < 3; PulseFile = numlist( ceil(end/2) ); end
if nargin < 2; fdw = 5000; end

  cou = dir('cou*.MDF');
  a = dir('B1*.MDF');
  b = dir('B2*.MDF');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD EXCITER SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
excdat = load('exc.mat');



if PulseFile ~= -1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND PULSE ON/OFF TIMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[B2 B2t] = readmdf(b(numlist(PulseFile)).name);
[t0 t1]  = pulse_limits(B2t, B2(:,8), 0.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET EXCITER AMPLITUDE ( =CURRENT)
% Get Currents from Current Monitor signals (look at fdw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(numlist)
  [Ex Et]  = readmdf(a(numlist(i)).name);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % GET START AND END INDEX
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  i_a = find(Et>t0); i_a = i_a(1); i_b = find(Et<t1); i_b = i_b(end);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % MEAN over 8 Rogowski coils
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:8
    Acoil(j) = std(Ex(i_a:i_b, j));
    [A_out(j), ph_out] = calibration_cm(j, fdw, Acoil(j), 0);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CURRENT AMPLITUDE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Curr(i) = mean(A_out);
end
end % IF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE KF-spectrograms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if kfon == 1
for i=1:length(numlist)
  savename1 = ['kfdata_' num_prezero(numlist(i), max(numlist)) '.mat'];
  disp(' ');
  disp(['evaluation of: ' cou(numlist(i)).name]);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % LOAD cou*.MDF FILES
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [A At] = readmdf(cou(numlist(i)).name);
  % POINTS PER FFT-WINDOW: f_abtast/fdw
  fwinpts = 100*round( (1/fdw)/(At(2)-At(1)) );
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CALCULATE KF-spectrogram
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [kf.t kf.m kf.f kf.mmat kf.fmat] = kfspectrogram(At, A', fwinpts, 60e3);
  save(savename1, 'kf');
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE SYNCHRONIZATION LIMITS SEMI-MANUALLY
% (Input: manually determined values, and then automatically determine
% exact values)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(numlist)
  savename  = ['arn_man' num_prezero(numlist(i), max(numlist)) '.mat'];
  savename1 = ['kfdata_' num_prezero(numlist(i), max(numlist)) '.mat'];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % TEST FILE ALREADY EXIST
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~exist(savename, 'file')
    load(savename1); savename1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT f- and m-spectrogram (FOR DETERMINATION OF f_sy-LIMITS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (PLOT ALSO WITH BOUNDARIES WHERE EXCITER IS OFF)
    figeps(16,20,1)
    subplot(211); subplotfspec(kf.t, kf.f, 20*log10(kf.fmat)); freezeColors
    subplot(212); subplotkspec(kf.t, kf.m, kf.mmat');
    print_adv([1 1], '-r300', ...
                      ['mode-f_' num_prezero(i,length(numlist)) '.eps']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD EXCITER VARIABLES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Aex   = excdat.cl(numlist(i), 3);
    mex   = excdat.cl(numlist(i), 5);
    fex_0 = excdat.cl(numlist(i), 6);
    fex_1 = excdat.cl(numlist(i), 7);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIND SYNCHRONIZATION RANGE manually
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    arn.ta = input('time where synchronization begins [ms]: ');
      arn.ta=arn.ta/1e3;
    arn.tb = input('time where synchronization ends [ms]: ');
      arn.tb=arn.tb/1e3;
    inp_fa = input('frequency where synchronization begins [Hz]: ');
    inp_fb = input('frequency where synchronization ends [Hz]: ');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DETERMINE FREQUENCIES MORE PRECISELY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIND TIME INDICES AND exact frequency positions from spline spectrum
      i_ta = findind(kf.t, arn.ta);
      i_tb = findind(kf.t, arn.tb);
      [arn.fa Aexact] = findpkinfspec(kf.f, kf.fmat(:,i_ta), inp_fa);
      [arn.fb Aexact] = findpkinfspec(kf.f, kf.fmat(:,i_tb), inp_fb);
        disp(['fsy_begin = ' num2str(arn.fa) ' [Hz]' ...
              'fsy___end = ' num2str(arn.fb) ' [Hz]']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RE-PLOT f- and m-spectrogram (WITH SYNCH BOUNDARIES)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (PLOT ALSO WITH BOUNDARIES WHERE EXCITER IS OFF)
    figeps(16,20,1)
    clf;
    subplot(211); subplotfspec(kf.t, kf.f, 20*log10(kf.fmat)); freezeColors
        line([arn.ta*1e3 arn.ta*1e3], [0 30e3], 'LineWidth', 1);
        line([arn.tb*1e3 arn.tb*1e3], [0 30e3], 'LineWidth', 1);
    subplot(212); subplotkspec(kf.t, kf.m, kf.mmat');
    print_adv([1 1], '-r300', ...
                ['synch_mf_' num_prezero(numlist(i),numlist(end)) '.eps']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if PulseFile ~= -1
        arn.curr = Curr(i);
    end
    arn.Aex  = Aex;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE OUTPUT TO mat-file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(savename, 'arn', 'kf');
    input('<< PRESS ANY KEY >>'); disp(' ');
    clf
  end       % ----- if ~exist
end         % ----- for

end