function [tvec mvec fvec mmat fmat] = kfspectrogram(tt, A, fwinpts, fend)
%20080225-Mo-04:04 Brandt
%function [tvec mvec fvec mmat fmat] = kfspectrogram(tt, A, fwinpts, fend)
% Spatiotemporal Matrix to kf-Spectrogram.
% IN: tt: time vector
%     A: spatiotemporal matrix (64 rows, n columns)
%     fwinpts: samples per of fwindow
%OUT: tvec: time step vector
%     mvec: mode axis
%     fvec: frequency axis
%     mmat: time-mspec-matrix
%     fmat: time-fspec-matrix [dBV]
% EX:

if nargin<4; fend=30e3; end

% Window Overlap
    olap = 0.5;
% TIME INTERVAL
    dt = tt(2)-tt(1);
% CALCULATE STEPS
    dwin = round(olap*fwinpts);
    steps = floor((length(tt) - fwinpts)/dwin);
% CREATE t-vector
    tvec = (( 1:steps )-1).*(dwin-1)*dt;


% METHOD 1: CALCULATE KF-SPECTRA via kfspec
    disp('calculate kfspectra ...');
    for i=1:steps
      %------display
      disp_num(i,steps);
      % make index
      ind = (1:fwinpts)+(i-1)*dwin;
      % STOP
      [fvec mvec kfmat] = kfspec(A(:,ind), tt(ind), [0 8], fend, 95, 0.5);
      % STOP
      fmat(:,i) = mean(kfmat);
      mmat(:,i) = mean(kfmat');
    end

% % METHOD 2: CALCULATE KF-SPECTRA via SU-FFT
%     for i=1:steps
%         ind = (1:winpts)+(i-1)*dwin;
%         [freq, A(:,i), pha] = su_fft(tt(ind), sig(ind), 1);
% %       ind = (1:winpts)+(i-1)*dwin;
% %       % STOP
% %       [freq mvec kfmat] = kfspec(A(:,ind), tt(ind), [0 8], fend, 50, 0.5);
% %       % STOP
% %       mmat(:,i) = mean(kfmat);
% %       fmat(:,i) = mean(kfmat');
%     end

% % Amat in dBV
%     Amat = 20.*log10( Amat );

end