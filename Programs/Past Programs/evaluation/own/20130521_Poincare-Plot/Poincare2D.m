function [px, py, info] = Poincare2D(s1, s2, s3, zeroplane)
%==========================================================================
%function [px, py] = Poincare2D(s1, s2, s3)
%--------------------------------------------------------------------------
% May-21-2013, C. Brandt, San Diego
% Poincare2D computes the Poincare plane from a 3D Poincare plot. The plane
% where the intersection points are taken must be provided by 'zeroplane'.
% zeroplane 1 for example means where s1 = 0.
%--------------------------------------------------------------------------
%INPUT:
% s1, s2, s3: signals 1 to 3 for the 3D plot
% zeroplane: number of the zero plane (1, 2 or 3)
%OUTPUT:
% px, py: intersection points on (x,y)-plane
%--------------------------------------------------------------------------
% EXAMPLE: 
%==========================================================================

ctr = 0;
% Define the zero plane
szero = eval(['s' num2str(zeroplane)]);
switch zeroplane
  case 1
    sp1 = s2;
    sp2 = s3;
    info = 's(t+0\tau)=0';
  case 2
    sp1 = s3;
    sp2 = s1;
    info = 's(t+1\tau)=0';
  case 3
    sp1 = s1;
    sp2 = s2;
    info = 's(t+2\tau)=0';
end

% Pre-allocation of output signals
% px = NaN(1,length(szero)-1); py = px;
for i = 2:length(szero)
  if sign(szero(i-1)) ~= sign(szero(i))
    ctr = ctr + 1;
    px(ctr) = ( sp1(i) + sp1(i-1) ) / 2;
    py(ctr) = ( sp2(i) + sp2(i-1) ) / 2;
  end
end

if ctr == 0
  px = []; py = [];
end

end