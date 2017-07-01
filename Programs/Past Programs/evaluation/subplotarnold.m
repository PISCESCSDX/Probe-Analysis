function [hxl hyl] = subplotarnold(f1, A1, f2, A2, xlb, ylb, markercol)
%function [hxl hyl] = subplotarnold(f1, A1, f2, A2, xlb, ylb, markercol)
%
% IN: f1, f2: frequency vector 1, 2
%OUT: A1, A2: amplitude vector
% EX:  subplotzebra(tt, fvec, W, 0);

fs = 12;

if nargin<7; markercol = 'k'; end;
if nargin<5; error('Input arguments are missing!'); end;

MarkSize = 2;

hold on
  plot(f1, A1, 'ko', 'MarkerSize', MarkSize, ... 
    'MarkerEdgeColor', markercol, 'MarkerFaceColor', markercol);
  plot(f2, A2, 'ko', 'MarkerSize', MarkSize, ...
    'MarkerEdgeColor', markercol, 'MarkerFaceColor', markercol);
hold off

if ~strcmp(xlb, '-3')
  [hxl hyl] = mkplotnice(xlb,ylb, fs, -20);
end
   
end