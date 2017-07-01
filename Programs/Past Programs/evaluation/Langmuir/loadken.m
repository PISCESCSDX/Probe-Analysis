function [fit,par,Rsqr,dii_Rsqr] = loadken(fname,plt)
% [FIT,PAR,R] = LOADKEN(FNAME) loads a complete dataset of a fittet probe 
% characteristic from "fname" and return the fit and the main paramters.
%
% [FIT,PAR, Rsqr] = LOADKEN(FNAME,PLT) returns the results and plots the data 
% for plt = 1.

if nargin == 1
  plt = 0;
end

if ~exist(strcat(fname,'.mod'),'file')
    tmp = load(strcat(fname,'.dat')); % old: '*.ken'
    fit = zeros(size(tmp,1),1).*NaN;
    par.pp = NaN;
    par.n1 = NaN;
    par.t1 = NaN;
    par.n2 = NaN;
    par.t2 = NaN;
    par.fp = NaN;
    Rsqr   = NaN;
    return
end

% Loading probe charakteristics and fit data
data    = load(strcat(fname,'.mod'));
uu      = data(:,1);
ii_fit  = data(:,2);
ii_filt = data(:,3);
ii      = data(:,4);

% Loading parameterfiles and assign values to te and ne
data   = load(strcat(fname,'.par'));
par.pp = data(2,1);
par.n1 = data(2,2);
par.t1 = data(2,3);
par.n2 = data(2,4);
par.t2 = data(2,5);
par.fp = data(2,6);

% Loading the derivatives filtered and modelled.
data   = load(strcat(fname,'.der'));
dii_filt = data(:,3);
dii_fit  = data(:,2);

% Loading the derivatives filtered and modelled.
data = load(strcat(fname,'.fit'));
uu_f      = data(:,1);
dii_fit2  = data(:,2);
dii_filt2 = data(:,3);

% Statistic values I-V-characteristic
ind = find(ii_fit);
if ~isempty(ind)
    sse = sum((ii(ind) - ii_fit(ind)).^2);
    sst = sum( (ii(ind) - mean(ii(ind))).^2 );
else
    sse = 1; sst = 1;
end
% quality of fit measure: 1-Chi-squared
Rsqr = 1-sse/sst;


% Statistic values dI/dV
ind = find(dii_fit);
if ~isempty(ind)
    dii_sse = sum( (dii_filt(ind) - dii_fit(ind)).^2 );
    dii_sst = sum((dii_filt(ind) - mean(dii_filt(ind))).^2);
else
    dii_sse = 1; dii_sst = 1;
end
% quality of fit measure: 1-Chi-squared
dii_Rsqr = 1 - dii_sse/dii_sst;

% Calculate a derivtive of a moving average charakteristic
dii = [ii(1);diff(smooth(ii,5))];

% Output
fit = ii_fit;

if plt==0
  return
end

% Create a legend with most importand parameters
if par.n2 ~= 0 && par.t2 ~=0
    str=cell(7,1);
    str(1) = {sprintf('fl-pot  = %1.2f_{ } V',par.fp)};
    str(2) = {sprintf('pl-pot = %0.1f V_{ }',par.pp)};
    str(3) = {sprintf('t_{e1}      = %0.2f eV',par.t1)};
    str(4) = {sprintf('n_{1}      = %1.2g m^{-3}',par.n1)};
    str(5) = {sprintf('t_{e2}      = %0.2f eV',par.t2)};
    str(6) = {sprintf('n_{2}      = %1.2g m^{-3}',par.n2)};
    str(7) = {sprintf('R-sqr = %1.3g (1=max)',1-sse/sst)};
else
    str=cell(5,1);
    str(1) = {sprintf('fl-pot  = %1.2f_{ } V',par.fp)};
    str(2) = {sprintf('pl-pot = %0.1f V_{ }',par.pp)};
    str(3) = {sprintf('t_{e1}      = %0.2f eV',par.t1)};
    str(4) = {sprintf('n_{1}      = %1.2g m^{-3}',par.n1)};
    str(5) = {sprintf('R-sqr = %1.3g (1=max)',1-sse/sst)};
end;

% Formatting yaxis and label
if max(abs(ii)) > 1E-6 && max(abs(ii)) < 1E-3
    yl = 'current (\muA)';
    cf = 1E6;
elseif max(abs(ii)) > 1E-3 && max(abs(ii)) < 1E-0
    yl = 'current (mA)';
    cf = 1E3;
else
    yl = 'current (A)';
    cf = 1;
end

% Plotting
fig = findobj('name','Probe Fit Result');
if isempty(fig)
    figure('numbertitle','off','name','Probe Fit Result')
    set(gcf,'menubar','none','toolbar','figure','papertype','a4',...
            'paperpositionmode','auto','units','pixels')
    set(gcf,'position',get(gcf,'position')+[1 -150 1 150])
else
    set(0,'currentfigure',fig)
    clf;
end
axes('position',[.16 0.38 .8 .55])
hold on
plot(uu,-ii.*cf,'.','markersize',1)
plot(uu,-ii_fit.*cf,'r')
plot(uu,-ii_filt.*cf)
set(gca,'fontsize',18)
ylabel(yl)
if 2*max(-ii) < max(-ii_fit) 
  mx=2*max(-ii);
else
  mx=1.2*max(-ii);
end
if max(ii) > .1*mx
  mn=-2*max(ii);
else
  mn=-.1*mx;
end
axis([round(min(uu)) round(max(uu)) mn.*cf mx.*cf])
ax = axis;
text((ax(2)-ax(1)).*.05+ax(1),ax(4)-(ax(4)-ax(3)).*.2,str','fontsize',14)
title(strcat('file://',pwd,'/',fname),'Interpreter','none','fontsize',14)


axes('position',[.16 0.10 .8 .22])
hold on
plot(uu,-dii*1E7,'.','markersize',1)
plot(uu,-dii_filt,'b')
plot(uu,-dii_fit,'r')
plot(par.pp-uu_f,-dii_filt2,'k')
plot(par.pp-uu_f,-dii_fit2,'m')
if 2*max(-dii_filt) < max(-dii_fit) 
  mx=2*max(-dii_filt);
else
  mx=1.2*max(-dii_filt);
end
if max(dii_filt) > .3*mx 
  mn=-2*max(dii_filt);
else
  mn=-.3*mx; 
end
axis([round(min(uu)) round(max(uu)) mn mx])
set(gca,'fontsize',18)
ylabel('dI/dU')
xlabel('voltage (V)')
puttextonplot(gca, 5, 80, ['dii-R-sqr= ' num2str(dii_Rsqr)]);

if isempty(Rsqr)
  Rsqr = -1;
end


end