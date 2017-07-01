function [hp] = plotprof(n_x, n_y, exp10, sty);
%function [] = plotprof(n_x, n_y, exp10);
% m-files needed:   subplotprof
% IN: n_x: x-vector
%     n_y: y-profile-vector
%     exp10: 10^(exp10)
%     sty: plotstyle, i.e. 'r-'
%OUT: hp: handle of the plot
% EX: plotprof(n_x, n_y)

if nargin<4; sty='ko' ; end;
if nargin<3; exp10=0 ; end;
if nargin<2; error('Input arguments are missing!'); end;


  figeps(8,12,1);
%-- make pictures look the same on printer and on screen
  wysiwyg;

  axes('Position', [0.20 0.18 0.7 0.75]);
  hp = subplotprof(n_x, n_y, exp10, sty);

end