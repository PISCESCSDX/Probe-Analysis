function [arn] = arnoldtongue_february(numlist, fex, kfon, flipon)
%20081106-Do-19:05 Brandt
%function [arn] = arnoldtongue_february(numlist, fex, kfon, flipon)
% IN: numlist: no. of cou0007.MDF files (e.g. = 7)
%     fex: [fstart fend]  in Hz
%     kfon: calculate kf-spectograms
%     flipon:  Flip Couronne Matrix if [A tt]=readmdf('couXXXX.MDF'); 
%              pcolor(A(1:1e3,:)'); shading flat
%              yields: ///!  Do not flip if it is \\\.
%OUT: Arnold tongue
% EX: 
% [arn] = arnoldtongue_february(1:1:20, 6.8e3, 15, Aex, fex, 1, 'flipon')

if nargin < 4; kfon=1; end
if nargin < 3; PulseFile = numlist( ceil(end/2) ); end

  cou = dir('cou*.MDF');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ### INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

fex_0 = fex(1);
fex_1 = fex(2);

% plasma pulse limit indices
  ind = 1:59e3;

load exc.mat
arn.Aex = cl(:,3);


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
  A  = A(ind, :);
  At = At(ind);
  % POINTS PER FFT-WINDOW: f_abtast/fdw
  % For Calculation of the spectrogram: Take 100 Wave-periods per
  % time-window!  =  100* period-length
  % OLD  fwinpts = 100*round( (1/fdw)/(At(2)-At(1)) );
  fwinpts = 1900;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CALCULATE KF-spectrogram
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  A=A';
  if strcmp(flipon, 'flipon'); A = flipud(A); end
  [kf.t kf.m kf.f kf.mmat kf.fmat] = kfspectrogram(At, A, fwinpts, 60e3);
  save(savename1, 'kf');
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE SYNCHRONIZATION LIMITS SEMI-MANUALLY
% (Input: manually determined values, and then automatically determine
% exact values)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=0;
status=0;
while status==0
  i=i+1;
  savename  = ['arn_man' num_prezero(numlist(i), max(numlist)) '.mat'];
  savename1 = ['kfdata_' num_prezero(numlist(i), max(numlist)) '.mat'];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % TEST FILE ALREADY EXIST
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~exist(savename, 'file')
    load(savename1);
    disp(['data: ' savename1])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT f- and m-spectrogram (FOR DETERMINATION OF f_sy-LIMITS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (PLOT ALSO WITH BOUNDARIES WHERE EXCITER IS OFF)
    figeps(16,20,1);
    subplot(211); subplotfspec(kf.t, kf.f, 20*log10(kf.fmat)); freezeColors
    subplot(212); subplotkspec(kf.t, kf.m, kf.mmat');
    % print_adv([1 1], '-r300', ...
    %                   ['mode-f_' num_prezero(i,length(numlist)) '.eps']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIND SYNCHRONIZATION RANGE manually
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    arn.ta = input('time where synchronization begins [ms]: ');
      arn.ta=arn.ta/1e3;
    arn.tb = input('time where synchronization ends [ms]: ');
      arn.tb=arn.tb/1e3;
%     inp_fa = input('frequency where synchronization begins [Hz]: ');
%     inp_fb = input('frequency where synchronization ends [Hz]: ');
    inp_fa = (fex_1 - fex_0)/kf.t(end)*arn.ta + fex_0;
    inp_fb = (fex_1 - fex_0)/kf.t(end)*arn.tb + fex_0;
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
    figeps(16,20,1);
    clf;
    subplot(211); subplotfspec(kf.t, kf.f, 20*log10(kf.fmat)); freezeColors
        line([arn.ta*1e3 arn.ta*1e3], [0 30e3], 'LineWidth', 1);
        line([arn.tb*1e3 arn.tb*1e3], [0 30e3], 'LineWidth', 1);
    subplot(212); subplotkspec(kf.t, kf.m, kf.mmat');
    % print_adv([1 1], '-r300', ...
    %             ['synch_mf_' num_prezero(numlist(i),numlist(end)) '.eps']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE OUTPUT TO mat-file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    in_again = input('<< Again?  y/n = 1/0 >>'); disp(' ');
    if in_again==1
      i=i-1;
    else
      save(savename, 'arn', 'kf');
    end
    clf
  end       % ----- if ~exist
  if i==length(numlist); status = 1; end
end         % ----- for

arn.Curr = 'empty, approximate from Aex';

end