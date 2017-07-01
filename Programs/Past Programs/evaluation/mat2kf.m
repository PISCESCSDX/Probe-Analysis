function [freq, mvec, kfmat, kfaxis] = mat2kf(A, tt, flim, mlim, wl, defch);
%function [freq, mvec, kfmat, kfaxis] = mat2kf(A, tt, flim, mlim, wl,
%defch);
% m-files needed:   readmdf, kfspec
% IN: A: mxn matrix
%     tt      time row
%     flim    limits of freq-axis: [fmin fmax] /Hz
%     mlim    limits of m-axis: [mmin mmax]
%     wl      window length
%     defch   vector with defect channels
%OUT: plot of the kf-spectrum
% EX: [freq, mvec, kfmat, kfaxis] = mat2kf(A,tt,[0 25e3],[-8 8],80);

if nargin<6; defch = []; end;
if nargin<5; wl = 80; end;
if nargin<4; mlim=[0 8]; end;
if nargin<3; flim=[0 25e3]; end;
if nargin<2; error('Input of A and tt is missing!'); end;

[freq mvec kfmat] = kfspec(A, tt, mlim, flim(2), wl, 0.5); 

flim = find(freq>=flim(1) & freq<flim(2));
freq=freq/1e3; fend=freq(end);
fa=freq(1); fb=fend; ma=mlim(1); mb=mlim(2);
kfmax = max(max(kfmat));
kfaxis= [fa fb ma mb 0 kfmax];

end