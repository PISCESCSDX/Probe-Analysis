function subplot4Barrows(P1, P2, Bx4, By4, lbl_vec, tvec)

    
    pos1= P1(: , 1);
    if pos1(end)== pos1(1), pos1= P1(1, :)'; end
    pos2= P2(1 , :)';
    if pos2(end)== pos2(1), pos2= P2(:, 1); end
    x_st= pos1(1)   - 10;
    x_en= pos1(end) + 10;
    y_st= pos2(1)   - 10;
    y_en= pos2(end) + 10;
    
    %------------
    x0= 0.12; y0= 0.12; dx= 0.38; dy= 0.38; spalt= 0.02;
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
        quiver(P1, P2, Bx4(:, :, i1), By4(:, :, i1));
        text(xpos, ypos, ['t= ', num2str(1e3.*tvec(i1),'%04.2f'),...
             ' ms'], 'FontSize', fs);
        axis([x_st, x_en, y_st, y_en]);
    end
    xlabel(ah(4), lbl_vec{1}, 'FontSize', fs);
    xlabel(ah(3), lbl_vec{1}, 'FontSize', fs);
    ylabel(ah(3), lbl_vec{2}, 'FontSize', fs);
    ylabel(ah(1), lbl_vec{2}, 'FontSize', fs);
    
    set(ah(3), 'FontSize', fs)
    set(ah(1), 'FontSize', fs, ...
               'XTick', []);
    set(ah(4), 'FontSize', fs, ...
               'YTick', []);
    set(ah(2), 'FontSize', fs, ...
               'XTick', [], 'YTick', []);
    
    drawnow;
end