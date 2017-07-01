function [fitfct, ymax] = fct_fitgauss(x,y,limits_lo, limits_up,start_val)
%==========================================================================
%function [fitfct, ymax] = fct_fitgauss(x,y,limits_lo, limits_up,start_val)
%--------------------------------------------------------------------------
% Fit for 5 parameter function  a*exp(-abs((x-b)./c).^d)+e
% Thus: limits_up, limits_lo and start_val must contain 5 values.
%--------------------------------------------------------------------------
% EX: 
%==========================================================================

if nargin<5; start_val = [0.5  0  4  2  0]; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make column vectors if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(x,2) > size(x,1)
  x = x';
end

if size(y,2) > size(y,1)
  y = y';
end

ymax = max(y);
y1 = y/ymax;
fitgauss = fittype('a*exp(-abs((x-b)./c).^d)+e');
opts = fitoptions(fitgauss);
set(opts,'start',start_val);
opts.Upper = limits_up;
opts.Lower = limits_lo;
fitfct = fit(x, y1, fitgauss, opts);

end