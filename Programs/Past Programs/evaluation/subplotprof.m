function [hp] = subplotprof(n_x, n_y, exp10, sty);
%function [hp] = subplotprof(n_x, n_y, exp10);
% IN: n_x: x-vector
%     n_y: y-profile-vector
%     exp10: 10^(exp10)
%     sty: plotstyle, i.e. 'r-'
% OUT: hp: handle of the plot
% EX:  subplotprof(n_x, n_y, 0, 'ko');

if nargin<4; sty='ko'; end;
if nargin<3; exp10=0; end;
if nargin<2; error('Input arguments are missing!'); end;

fonts = 14;

% PLOT
    hp = plot(n_x, n_y, sty);
    box on;
    set(gca, 'fontsize', fonts);
% LABELS
  xlabel('time [ms]', 'Fontsize', fonts);
  hyl = ylabel(ylbl, 'Fontsize', fonts);
% set ticklength
  set(gca, 'ticklength', [0.02 0.025]);

end