function kenprofiles(ne, pp, te, rsqr)
% function kenprofiles(ne, pp, te)
% Plots and the densiyt. plasma potential and temperature in one figure,
% and saves it as eps. Without radius vector.

if nargin < 4; 
  rsqr = -1; 
else
  a=find(rsqr>1 | rsqr<0.7);
  ne(a) = 0;
  pp(a) = 0;
  te(a) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ne, phi_p, Te
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fonts = 12;
figx = 20;
epsx = round(7.5*figx);
figeps(figx, 7, 1);
axpos = axes_positioning(gcf, 1/3*[1 1 1], 1, [1.0 2.0 1.5 0.5], 2.0);

% PLOT ne
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('position', axpos{1})
plot(ne, 'ko', 'Markersize', 2)
set(gca, 'xlim', [1 length(ne)]);
set(gca, 'Color', 'none')
mkplotnice('-1', 'n_e [m^{-3}]', fonts, -30);
puttextonplot(gca, 85, 90, '(a)', 0, fonts);

% PLOT phi_p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('position', axpos{2})
plot(pp, 'ko', 'Markersize', 2)
set(gca, 'xlim', [1 length(pp)]);
mkplotnice('-1', '\phi_p [V]', fonts, -30);
puttextonplot(gca, 85, 90, '(b)', 0, fonts);

% PLOT T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('position', axpos{3})
plot(te, 'ko', 'Markersize', 2)
set(gca, 'xlim', [1 length(te)]);
mkplotnice('-1', 'T_e [eV]', fonts, -30);
puttextonplot(gca, 85, 90, '(c)', 0, fonts);

% PRINT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print('-depsc2', '-r300', ['profiles_x' num2str(epsx) '.eps']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end