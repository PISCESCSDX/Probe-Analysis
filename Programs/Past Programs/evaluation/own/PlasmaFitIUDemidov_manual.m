function iudata = PlasmaFitIUDemidov_manual(x,xind,y,ye,Dye, ...
  ne_prefac,Iisat_max,m_i,iudata) %#ok<INUSL>
% Aug-09-2013, Christian Brandt, San Diego


disp(' --- Starting Manual Fit Mode --- ')
status_quit = 0;
% Load natural constants
natconst

% Assign values to this m-file's fit parameters
if isnan(iudata.ne)
  disp('  Fit Info: iudata is NaN, you need to start from scratch.')
  iudata.Te = 1;
  iudata.ne = 1e16;
  iudata.Vp = 20;
  iudata.Vf =  0;
  iudata.RsqrIe = +Inf;
  iudata.RsqrDIe = +Inf;
  iudata.Fli = 0;
  iudata.ni = Iisat_max/(0.61*e0*iudata.probe.A*sqrt(e0*iudata.Te/m_i));
end

fit.Te = iudata.Te;
fit.ne = iudata.ne / ne_prefac;
fit.vp = iudata.Vp;

%==========================================================================
% While 
while status_quit==0
  
 par = input(['Change parameter (Enter: quit, 1:ne, 2:Vp, 3:Te, 4:Vf,' ...
              ' 5:load par, 9:no output): ']);
 % Break if Quitting
 if isempty(par)
   load('xchange.mat')
   % Save all data to file (can be used to fit different later)
   
   savefn = mkstring([num2str(xchange.shotno) '_fit'],'0',xchange.i, ...
            xchange.Ncharact,'.mat');
   save(savefn,'x','xind','y','ye','Dye','ne_prefac','Iisat_max', ...
               'm_i','iudata');
   break
 end
 
 % No Output
 if par==9
     PlasmaFitIUDemidov_nooutput
     break
 end

 % Load parameters of IV-charact. before
 if par==5
  % Load saved variables from open function: csdxprobeseep.m
  load('xchange.mat')
  savefn = mkstring([num2str(xchange.shotno) '_fit'],'0',xchange.i-1, ...
            xchange.Ncharact,'.mat');
  
  % Assign Values
  paLoad = load(savefn);
  fit.ne = paLoad.iudata.ne / ne_prefac;
  iudata.ne = paLoad.iudata.ne;
  fit.vp = paLoad.iudata.Vp;
  iudata.Vp = fit.vp;
  fit.Te = paLoad.iudata.Te;
  iudata.Te = fit.Te;
  iudata.Vf = paLoad.iudata.Vf;
  
 else
   
       % Enter new parameters value
       if par==4
         figure(4); plot(x,y,'-x'); title('manual V_f detection')
       end
       par2= input('  Enter new value: ');
       close(figure(4))

       % Put new parameter to store variable
       switch par
         case 1
           fit.ne = par2 / ne_prefac;
           iudata.ne = par2;
         case 2
           fit.vp = par2;
           iudata.Vp = fit.vp;
         case 3
           fit.Te = par2;
           iudata.Te = fit.Te;
         case 4
           iudata.Vf = par2;
       end
 end

%---------------------------------------------- Perform new fit calculation
  ivp = findind(x,fit.vp);
  xifitres = 1:ivp;
  xx =   x(xifitres);
  yy = Dye(xifitres);
% Newly define Fit function
  vp = fit.vp;
  ffun = @(c,x) c(1)*(vp-x)./c(2).*exp(-(vp-x)./c(2));
% Calculate the Optimal Fit of DIe (=Dye)
  Dyefit = ffun([fit.ne fit.Te],xx);
% Calculate the Optimal Fit For Ie
  [~, Iefit] = int_discrete(xx, Dyefit );
% Extract Confidence Intervals
  a0 = [fit.ne fit.Te];
  options = statset('FunValCheck', 'Off');
  [a,r,~,COVB] = nlinfit(xx,yy,ffun,a0,options);
c = nlparci(a,r,'covar',COVB);
% Recalculate fit regression coefficients
Dyefit = ffun([fit.ne fit.Te],xx);
    sef = sum( (yy-Dyefit).^2 );
    sdm = sum( (yy-mean(yy)).^2 );
  iudata.RsqrDIe = 1 - sef/sdm;       % The larger the better
% Recalculate fit regression coefficients
yy2 = ye(xifitres);
    sef = sum( (yy2-Iefit).^2 );
    sdm = sum( (yy2-mean(yy2)).^2 );
iudata.RsqrIe = 1 - sef/sdm;       % The larger the better
%---------------------------------------------- Calculate Plasma Parameters
iudata.ni = Iisat_max/(0.61*e0*iudata.probe.A*sqrt(e0*iudata.Te/m_i));
iudata.neerr = [c(1,2) c(1,1)] *ne_prefac;
iudata.Teerr = [c(2,1) c(2,2)];
iudata.Vperr = [NaN NaN];
iudata.Fli = 0.61*iudata.ni*sqrt(e0*iudata.Te/m_i);
%----------------------------------------------------------------- Plot Fit
PlasmaFitIUDemidov_plot
end
%==========================================================================

end