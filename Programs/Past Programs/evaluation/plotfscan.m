function [h1, h2] = plotfscan(fdiff, fvec, W, xlb, ylb, clb, fs, fftlim)
%function [h1, h2] = plotfscan(fdiff, fvec, W, xlb, ylb, clb, fs, fftlim)
% pcolor-plot of a discrete f-scan
% IN: fdiff:  frequency difference vector
%     fvec:   frequency vector /Hz
%     W:      amplitude matrix(fdiff, fvec)
% output:   h1.a: handle axis pcolor, h1.p: handle pcolor-plot
%           h2.a: handle axis colorbar, h2.p: handle colorbar-plot
% EXAMPLE:  [h1, h2] = plotfscan(fdiff, fvec, W, xlb, ylb, clb)

if nargin<8; fftlim = [];  end
if nargin<7; fs = 12;  end
if nargin<6; clb = ''; end
if nargin<5; ylb = ''; end
if nargin<4; error('Input arguments are missing!'); end

figeps(12,11,1)
    axes('Position', [0.15 0.18 0.68 0.77]);
    [h1, h2] = subplotfscan(fdiff, fvec, W, xlb, ylb, clb, fs, fftlim);
    
end