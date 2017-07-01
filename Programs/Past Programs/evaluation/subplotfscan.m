function [h1, h2] = subplotfscan(fdiff, fvec, W, xlb, ylb, clb, fs, fftlim)
%function [ha, hp] = subplotfscan(fdiff, fvec, W, xlb, ylb, clb, fs, fftlim)
% pcolor-plot of a discrete f-scan
% IN: fdiff:  frequency difference vector
%     fvec:   frequency vector /Hz
%     W:      amplitude matrix(fdiff, fvec)
% output:   ha: handle axis, hp: handle pcolor-plot
% EXAMPLE:  [ha, hp] = subplotfscan(fdiff, fvec, W, xlb, ylb, clb)

if nargin < 8 | isempty(fftlim);
    fftlim(1) = min(min(W));
    fftlim(2) = max(max(W));
end;
if nargin < 7; fs = 12; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCOLOR PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1.p = pcolor(fdiff, fvec, W);
h1.a = gca;
shading flat
mkplotnice(xlb, ylb, fs, -20);
set(h1.a, 'clim', fftlim);

if ~strcmp(clb, '-1')
  h2.p = our_colorbar(clb, 14, 10, 0.015, 0.010);
else
  h2.p = 0;
end
h2.a = gca;


end