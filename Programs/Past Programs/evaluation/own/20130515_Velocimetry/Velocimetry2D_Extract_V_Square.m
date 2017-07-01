function [cpos,vvec] = Velocimetry2D_Extract_V_Square(povec,dtvec)
%==========================================================================
%function [cpos,vvec] = Velocimetry2D_Extract_V_Square(posvec,dtvec)
%--------------------------------------------------------------------------
% May-16-2013, C. Brandt, San Diego
% Velocimetry2D_Extract_V_Square calculates the velocity vector 'vvec' at
% the center position 'cpos' of four pixels at the positions 'povec'.
% Between two pixel pairs the time delay 'dtvec' is needed.
% --
% The assumptions for the used equations are that an object passes the 
% four pixels with a constant velocity and phase front of the object is
% straight (not bent), like a infinite wave front. This means the object
% is assumed to be much larger than the covered area of the pixels.
%--------------------------------------------------------------------------
%IN:
% povec: structure with (x,y)-positions p1, p2, p3, p4 of the chosen pixels
% dtvec: structure with the 2 time delays tab, tcd
%        (whereas each time delay means: tab=t_b-t_a, and so on)
%OUT:
% cpos: center position of velocity vector
% vvec: 2D velocity vector components (at position posvec)
%--------------------------------------------------------------------------
%EXAMPLE: (velocity vector assumed: v = 86.6ex -50ey)
% povec.p1=[-50   0];
% povec.p2=[+50   0];
% povec.p3=[  0 -50];
% povec.p4=[  0 +50];
% dtvec.tab = +0.866;
% dtvec.tcd = -0.500;
% [cpos,vvec] = Velocimetry2D_Extract_V_Square(povec,dtvec);
%EXAMPLE: (velocity vector assumed: v = 100ex)
% povec.p1=[-50   0];
% povec.p2=[+50   0];
% povec.p3=[  0 -50];
% povec.p4=[  0 +50];
% dtvec.tab =  1.0;
% dtvec.tcd =  0.0;
% [cpos,vvec] = Velocimetry2D_Extract_V_Square(povec,dtvec);
%==========================================================================

% Shorten the variables
tab = dtvec.tab;
tcd = dtvec.tcd;
ax  = povec.p1(1); ay = povec.p1(2);
bx  = povec.p2(1); by = povec.p2(2);
cx  = povec.p3(1); cy = povec.p3(2);
dx  = povec.p4(1); dy = povec.p4(2);

% OUTPUT: Calculate the center position of the velocity vector
cpos = (povec.p1 + povec.p2 + povec.p3 + povec.p4)/4;

% Calculate the velocity in y-direction
% Calculate Prefactor of Vx = F * Vy
F = (tab*(-cy+dy)-tcd*(-ay+by)) ./ (tcd*(-ax+bx)-tab*(-cx+dx));
% Calculate Vy
Vy = ((-ax+bx)*F + (-ay+by)) ./ (tab.*((F.^2) + 1));
% Calculate Vx
Vx = F.*Vy;

% OUTPUT: Velocity Vector
vvec = [Vx, Vy];

end