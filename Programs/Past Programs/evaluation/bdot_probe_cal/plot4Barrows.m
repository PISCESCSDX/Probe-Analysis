function [] = plot4Barrows(Bmat, fkHz, saveps)
%function [] = plot4Barrows(Bmat, fkHz, saveps)
%NEED:
% IN:
%OUT:
% EX:

if nargin < 3; saveps = 0; end
if nargin < 2; error('Input arguments are missing!'); end;

% fig(9,6,1);
%   axes('Position', [0.17 0.22 0.74 0.7]);
%   subplotzebra(tt, phi, A);

    lblvec= {'x [mm]', 'y [mm]'};
    phi= [0:1:3]'.*-0.5.*pi;
    ph_fac= exp(i.*phi);
    for i1= 1:4
        Bx4(:, :, i1)= real(Bmat.Bx .*ph_fac(i1));
        By4(:, :, i1)= real(Bmat.By .*ph_fac(i1));
    end
    tvec= [0:1:3].*0.25./(fkHz.*1e3);
    subplot4Barrows(Bmat.P1, Bmat.P2, Bx4, By4, lblvec, tvec)

   if saveps == 1
%      print_adv([0 1 1 1 1], '-r100', 'Ipar.eps', 95);
     print('-depsc2', '-r300', 'Barrows.eps');
   end

end