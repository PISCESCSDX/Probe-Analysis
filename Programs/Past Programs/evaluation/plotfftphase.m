function [h1, h2] = plotfftphase(frq, Axy, Phxy, PhStd, Aylim, phlabel)
%function [] = plotfftphase(frq, Axy, Phxy, PhStd)
% m-files needed:
% IN: frq: freq axis [kHz]
%     Axy:  CPSD amplitude
%     Phxy: CPSD phase
%     PhStd: optional phase standard deviation
%     Aylim:  [ylim1 ylim2]
%OUT: plot of the cpsd-spectrum and cpsdphase-spectrum
%     h1.a handle axis1, h1.p handle plot1 ...
% EX: plotfftphase(frq, Axy, Phxy, PhStd)

if nargin<6 || isempty(phlabel); phlabel='phase (n,\phi) [\pi]'; end;
if nargin<2; error('Input arguments are missing!'); end

fonts = 12;

figeps(10,9,1); clf
    axes('Position', [0.22 0.60 0.70 0.36]);
    [h1.a h1.p] = subplotfft(frq, Axy, ' ', 'S [dB]', [1 60], Aylim);
    
    axes('Position', [0.22 0.15 0.70 0.36]);
    [h2.a h2.p] = subplotphase(frq, Phxy, PhStd, '', phlabel, [1 60]);

end