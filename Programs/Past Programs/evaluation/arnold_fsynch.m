function [f_synch Mratio Aratio] = arnold_fsynch(tvec,mvec,fvec,mmat,fmat,fdw,fexvec,mex)
%20080225-Mo-12:48 Brandt
%function [fsynch] = arnold_fsynch(tvec,mvec,fvec,mmat,fmat,fdw,fexvec,mex)

% STEPS OF TVEC
  tl = length(tvec);
% MEX
  mex = abs(mex);

% # INPUT: Aratio LIMIT (5 means Aex must be 5 times higher than DW for
% Synch
  Ar_lim = 2.5;


% FREQUENCY RESOLUTION
  df = fvec(2)-fvec(1);
  disp([ 'frequency resolution: ' num2str(df) '[ Hz]'] )

% CHECK NOW for every tvec-step:%
% 1: A_mex > Am_max
%   IF TRUE ==> 2.
%   ELSE 4.
% 2. fex is maximal PEAK
%   IF TRUE ==> 3.
%   ELSE 4.
% 4. 

% DELTA f
  df = fvec(2)-fvec(1);
  resf = 10;
  sdf=round(df/resf);
  fv = fvec(1):sdf:fvec(end);

% CORRECT EXCITER FREQUENCY - GET ca.-OFFSETS
% ONLY SENSFUL IF MISMATCH IS BIG AT BEGIN
% if delta_f < df
  fd      = fdw-2*df:1:fdw+2*df;
  Adw     = spline(fvec,fmat(:,1), fd);
  [A i_A] = max(Adw);
  fdw0    = fd(i_A);
  % CORRECTION fex
  fe      = fexvec(1)-2*df:1:fexvec(1)+2*df;
  Aex     = spline(fvec,fmat(:,1), fe);
  [A i_A] = max(Aex);
  fex     = fe(i_A);
delta_fex = fex - fexvec(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGINN: ROUTINE DER SYNCHRONISATIONS-BESTIMMUNG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:tl
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VERSUCH NEU: VERHÄLTNIS EXCITER PEAK/ DW-PEAK
    % EXAKTE BESTIMMUNG DER PEAKHÖHEN
    % CORRECT EXCITER FREQUENCY - GET ca.-OFFSETS
    % ONLY SENSFUL IF MISMATCH IS BIG AT BEGIN
    % if delta_f < df
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % FIND MAX PEAK IN VICINITY OF DRIFT WAVE
      fdw = fdw0; status = 0;
      while status==0
        fd  = fdw-df/2:1:fdw+df/2; 
        Adw     = spline(fvec,fmat(:,i), fd);
        [AdwMax(i) i_A] = max(Adw);
          if i_A==1
            status = 0;
            fdw = fdw-df/2;          % shift new fdw to the left about df/2
          else
            if i_A==length(Adw)
              status = 0;
              fdw = fdw+df/2;       % shift new fdw to the right about df/2
            else
              status = 1;
              fdw = fd(i_A);
            end
          end
      end
    % CORRECTION from TOP: delta_fex
    fexvec(i) = delta_fex + fexvec(i);
    fe      = fexvec(i)-df/2:1:fexvec(i)+df/2;
    Aex     = spline(fvec,fmat(:,i), fe);
    [AexMax(i) i_A] = max(Aex);
    fex     = fe(i_A);
    % MAXIMUM OF SPECTRUM
    f_range = 500:1:30e3;
    [fspecMAX.A i_A] = max( spline(fvec,fmat(:,i), f_range) );
     fspecMAX.f = f_range(i_A);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PLOT m-MATRIX AND f-MATRIX
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figeps(10,10,1); clf;
  subplot(211); plot(mvec, mmat(:, i) )
  subplot(212); plot(fvec, fmat(:, i), 'r')
    xlm=get(gca, 'xlim'); ylm=get(gca, 'ylim');
  text( xlm(1)+0.8*(xlm(end)-xlm(1)), ylm(1)+0.8*(ylm(2)-ylm(1)), num2str(i));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SORT MODE MATRIX: mmat, and find mex amplitude
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [mA i_m] = sort( mmat(:, i) );
  % FIND INDEX of mvec with mvec=mex
      i_mex = find(mvec==mex);
  % FIND index in sorted index vector i_m of exciter-mode-nr mex
      i_e = find(i_m==i_mex);
  % m-EX/m-DW-RATIO: exciter mode amplitude / maximum amplitude
      Mratio(i) = mA(i_e)./mA(end);
      AmpEx(i) = mA(i_e);
      AmpMax(i)= mA(end);
  % AMPLITUDE RATIO
    Aratio(i) = AexMax(i)/AdwMax(i);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % IF EXCITER MODENR >= 85% of MAX MODENR ==> SYNCHRONIZATION
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if Mratio(i) >= 0.95
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IF RATIO Aex/Adw>=1 OR 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OLD:
    % [fA i_f] = max( fmat(:, i) );
    % % df-factor DELTA f factor
    % df_1 = 1.5;
    % if fexvec(i) < fvec(i_f)+  df_1*df && fexvec(i) > fvec(i_f)-df_1*df
    %
    % # INPUT Aratio > Ar_lim
    if ( (Aratio(i)> Ar_lim ) & (abs(fex-fspecMAX.f)<5) ) | ( (abs(fex-fdw)<5) & (abs(fex-fspecMAX.f)<5) )
      f_synch(i) = fex;
    else
      f_synch(i) = NaN;
    end
  %
  else
    f_synch(i) = NaN;
  end
%
%   [ (1:length(f_synch))' f_synch' Mratio' AexMax' AdwMax']
end

    [ (1:length(f_synch))' f_synch' Mratio' Aratio']
    figeps(10,10,2);  
      plot(Aratio.*Mratio);
    figeps(10,10,2);
      plot((Aratio.^2).*(Mratio.^20) - std((Aratio.^2).*(Mratio.^20)) )
%     input(' <<press any key>>');

end