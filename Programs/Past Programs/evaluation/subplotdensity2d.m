function subplotdensity2d(xvec, yvec, mat, limc, xlb, ylb, fonts)
%

pcolor(xvec, yvec, mat);
shading interp;
colormap cbpastell(2);
caxis([-limc limc]);

mkplotnice(xlb, ylb, fonts, -30);
set(gca, 'ticklength', 0.02*[1 1]);
axis square;

end