function [] = plot4Ipar(X, Y, fkHz, Bmat, Ipar)
%function [] = plot4Ipar(X, Y, fkHz, Bmat, Ipar)
%NEED:
% IN:
%OUT:
% EX:
    fs= 15;
  % INTERPOLATE DATA
    X= interp2(X, 2);
    Y= interp2(Y, 2);
    Ipar= interp2(Ipar, 2);
  % SET LABEL VEC
    lblvec= {'x [mm]', 'y [mm]', 'j [mA/cm^2]'};
    phi= [0:1:3]'.*-0.5.*pi;
    ph_fac= exp(i.*phi);
    for i1= 1:4
        Ipar_4(:, :, i1)= real(Ipar .*ph_fac(i1));
    end
    tvec= [0:1:3].*0.25;
    % t in s: tvec= [0:1:3].*0.25./(fkHz.*1e3);
    subplot4Ipar(X, Y, Ipar_4, lblvec, tvec);


function subplot4Ipar(P1, P2, QuaMat, lbl_vec, tvec)
    pos1= P1(: , 1);
    if pos1(end)== pos1(1), pos1= P1(1, :)'; end
    pos2= P2(1 , :)';
    if pos2(end)== pos2(1), pos2= P2(:, 1); end
    

    clf;
    i1=1;
        ah(i1)= axes('Position', apos(i1, :));
        pcolor(P1, P2, QuaMat(:, :, i1));
        axis equal
        shading flat
        caxis([-1.*maxamp, maxamp]');
%        text(xpos, ypos, ['\phi= ', num2str(tvec(i1).*360, '%3.0f'), ' °'], 'FontSize', fs);
        text(xpos, ypos, ['t= ', num2str(tvec(i1).*1, '%1.2f'), 'T'], 'FontSize', fs);
    colormap(pastell(128));
    xlabel(ah(4), lbl_vec{1}, 'FontSize', fs);
    xlabel(ah(3), lbl_vec{1}, 'FontSize', fs);
    ylabel(ah(3), lbl_vec{2}, 'FontSize', fs);
    ylabel(ah(1), lbl_vec{2}, 'FontSize', fs);
    ah(5)= our_colorbar(lbl_vec{3}, fs, 3);

    set(ah(3), 'FontSize', fs)
    set(ah(1), 'FontSize', fs, 'XTick', []);
    set(ah(4), 'FontSize', fs, 'YTick', []);
    set(ah(2), 'FontSize', fs, 'XTick', [], 'YTick', []);
    set(ah(5), 'position', [x0+2*dx+1.5*spalt, y0, 0.018, 2*dy+spalt])
    drawnow;
end