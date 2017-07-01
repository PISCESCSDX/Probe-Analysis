function [] = plot4Ipar(X, Y, fkHz, Bmat, Ipar, saveps)
%function [] = plot4Ipar(X, Y, fkHz, Bmat, Ipar, sav)
%NEED:
% IN:
%OUT:
% EX:

if nargin < 6; saveps = 0; end
if nargin<5; error('Input arguments are missing!'); end;

figeps(14,10,1);
    X= interp2(X, 2);
    Y= interp2(Y, 2);
    Ipar= interp2(Ipar, 2);

    lblvec= {'x [mm]', 'y [mm]', 'j [mA/cm^2]'};
    phi= [0:1:3]'.*-0.5.*pi;
    ph_fac= exp(i.*phi);
    for i1= 1:4
        Ipar_4(:, :, i1)= real(Ipar .*ph_fac(i1));
    end
    tvec= [0:1:3].*0.25;
    % t in s: tvec= [0:1:3].*0.25./(fkHz.*1e3);
    subplot4Ipar(X, Y, Ipar_4, lblvec, tvec);

  if saveps == 1
    print_adv([0 1 1 1 1], '-r100', 'Ipar.eps', 95);
  end

end