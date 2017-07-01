function [h1,h2,h3] = plotfft2colors(frq, A1, A2, xlb, ylb, fftlim, fs, color1, color2)
%function [h1,h2,h3] = subplotfft2colors(frq, A1, A2, xlb, ylb, fftlim, fs, color1,
%color2)
% IN: frq: frq axis /kHz
%     A1: amplitude (decadic)
%     A2: amplitude which will be shown in the background
%     xlb:
%     ylb:
%     fftlim:
%     fs: FontSize
%     color1, color2
%OUT: duo-color areaplot of 2 fft-spectra
% EX:

if nargin<3; help plotfftduo; return; end;
if nargin<4; xlb = 'f [kHz]'; end;
if nargin<5; ylb = 'S [dB]'; end;
if nargin<6; fftlim = [min(A1) max(A1)]; end;
if nargin<7; fs=12; end;
if nargin<8; color1=[]; end;
if nargin<9; color2=[]; end;

figeps(12,6,1)
    axes('Position', [0.15 0.22 0.68 0.75]);
    [h1,h2,h3] = subplotfft2colors(frq, A1, A2, xlb, ylb, fftlim, fs, color1, color2);
    
end