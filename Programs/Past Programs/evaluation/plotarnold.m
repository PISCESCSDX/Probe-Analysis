function [] = plotarnold(f1, A1, f2, A2, xlb, ylb);
%function [] = plotarnold(f1, A1, f2, A2, xlb, ylb);
%
% IN: f1, f2: frequency vector 1, 2
%OUT: A1, A2: amplitude vector
% EX: 

if nargin<6; error('Input arguments are missing!'); end;

figeps(8,8,1)
    axes('Position', [0.22 0.20 0.7 0.75]);
    subplotarnold(f1, A1, f2, A2, xlb, ylb);

end