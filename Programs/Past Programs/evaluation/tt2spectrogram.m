function [tf frq Amat] = tt2spectrogram(tt, sig, winpts, fend)
%20080225-Mo-03:03 Brandt
%function [tf frq Amat] = tt2spectrogram(tt, sig, winpts, fend)
%TimeTraceToSpectrogram
% IN: tt: time vector
%     sig: signal
%     win_le: samples per of window
%OUT: tf: time step vector
%     frq: frequency axis
%     Amat: spectrogram amplitudes [dBV]
% EX:

if nargin<4; fend=30e3; end


% Window Overlap
    olap = 0.5;
% TIME INTERVAL
    dt = tt(2)-tt(1);

% CALCULATE STEPS
    dwin = round(olap*winpts);
    steps = floor( (length(tt) - winpts)/dwin );

% CREATE tf-vector
    tf = (( 1:steps )-1).*(dwin-1)*dt;

% CALCULATE SPECTRA
    for i=1:steps
        ind = (1:winpts)+(i-1)*dwin;
        [freq, A(:,i), pha] = su_fft(tt(ind), sig(ind), 1);
    end

    %-------CUT to fend
    i_f = find(freq<fend);
    frq = freq(i_f);
    Amat = A(i_f,:);

% Amat in dBV
    Amat = 20.*log10( Amat );

end