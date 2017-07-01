function PlasmaFitIUDemidov_manual_v01


% ----------------------------- Calculate the Optimal Fit of DIe (=Dye)
  xifitres = xindfit(1:end-jvec(iRmax));
  xx = x(xifitres);
  yy = Dye(xifitres);
% Newly define Plasma Potential
  vp = x(xifitres(end));
% Newly define Fit function
  ffun = @(c,x) c(1)*(vp-x)./c(2).*exp(-(vp-x)./c(2));
% Calculate the Optimal Fit of DIe (=Dye)
  Dyefit = ffun(ffit{iRmax},xx);
% Calculate the Optimal Fit For Ie
  [~, Iefit] = int_discrete(xx, ffun(ffit{iRmax},xx) );
% Extract Confidence Intervals
  a0 = [ffit{iRmax}(1) ffit{iRmax}(2)];
  options = statset('FunValCheck', 'Off');
  [a,r,~,COVB] = nlinfit(xx,yy,ffun,a0,options);
c = nlparci(a,r,'covar',COVB);