function  [xf,yf,fitpar] = fit_gauss2(x,y,xn)
%function  [xf,yf] = fit_gaussbent(x,y,xn)
%  x: x-vector
%  y: y-vector
% xn: x-vector for fit

if nargin <3; xf = interp(x, 10); else xf = xn; end
if nargin==0; help fitgauss; return,  end


% USE NORMALIZED DATA FOR FIT
fitpar.norm = max(abs(y));
yn = y/fitpar.norm;

% Find position of maximum
[~,pos] = max(y);
x0 = x(pos);

fitpar.eq = 'a*exp(-abs((x-b)/c).^d)+e';
ft = fittype(fitpar.eq);
sp =[max(yn), x0,     min(x), 2,  0];
up =[1.1, x0+1.1, x0+5.1,  2.1,  +1];
low=[0.0, x0-1.1, 0     ,  2.0,  -1];
cf = fit(x, yn, ft, 'startpoint', sp, 'lower', low, 'upper', up);
a=cf.a;
b=cf.b;
c=cf.c;
d=cf.d;
e=cf.e;

fitpar.abcde = [a,b,c,d,e];
fitpar.FWHM = 2*c;

fprintf('FWHM: %f\n',          fitpar.FWHM);
fprintf('fitted center: %f\n', b);

x = xf;
yf = fitpar.norm*eval(fitpar.eq);

end