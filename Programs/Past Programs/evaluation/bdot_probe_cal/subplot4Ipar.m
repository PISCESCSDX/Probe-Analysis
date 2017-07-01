function subplot4Ipar(P1, P2, QuaMat, lbl_vec, tvec)

figeps(14,10,1);
    minval= 0;
    maxval= 0;
    minval= min(min(min(min(QuaMat))), minval);
    maxval= max(max(max(max(QuaMat))), maxval);
    maxamp= max(abs(maxval), abs(minval));
    
    pos1= P1(: , 1);
    if pos1(end)== pos1(1), pos1= P1(1, :)'; end
    pos2= P2(1 , :)';
    if pos2(end)== pos2(1), pos2= P2(:, 1); end
    
    %------------
    x0= 0.14; y0= 0.18; dx= 0.35; dy= 0.38; spalt= 0.02;
    apos= [x0, y0+dy+spalt, dx, dy; x0+dx+spalt, y0+dy+spalt, dx, dy; ...
           x0, y0, dx, dy; x0+dx+spalt, y0, dx, dy];
    xpos= pos1(1) + 0.05.*(max(pos1) - min(pos1));
    ypos= pos2(1) + 0.9.*(max(pos2) - min(pos2));
    hlxpos= pos1(1) - 0.5.*(max(pos1) - min(pos1));
    hlypos= pos2(1) + 2.1.*(max(pos2) - min(pos2));
    
    fs= 15;
    clf;
    for i1= 1:4
        ah(i1)= axes('Position', apos(i1, :));
        pcolor(P1, P2, QuaMat(:, :, i1));
        axis equal
        shading flat
        caxis([-1.*maxamp, maxamp]');
%        text(xpos, ypos, ['\phi= ', num2str(tvec(i1).*360, '%3.0f'), ' °'], 'FontSize', fs);
        text(xpos, ypos, ['t= ', num2str(tvec(i1).*1, '%1.2f'), 'T'], 'FontSize', fs);
    end
    colormap(pastell(128));
    xlabel(ah(4), lbl_vec{1}, 'FontSize', fs);
    xlabel(ah(3), lbl_vec{1}, 'FontSize', fs);
    ylabel(ah(3), lbl_vec{2}, 'FontSize', fs);
    ylabel(ah(1), lbl_vec{2}, 'FontSize', fs);
    ah(5)= our_colorbar(lbl_vec{3}, fs, 3);
    
    set(ah(3), 'FontSize', fs)
    set(ah(1), 'FontSize', fs, ...
               'XTick', []);
    set(ah(4), 'FontSize', fs, ...
               'YTick', []);
    set(ah(2), 'FontSize', fs, ...
               'XTick', [], 'YTick', []);
    set(ah(5), 'position', [x0+2*dx+1.5*spalt, y0, 0.018, 2*dy+spalt])
    
    drawnow;
end