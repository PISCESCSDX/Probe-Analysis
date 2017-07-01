function [cpos,vvec] = Velocimetry2D_Extract_V_Triangle(povec,dtvec)
%==========================================================================
%function [cpos,vvec] = Velocimetry2D_Extract_V_Triangle(posvec,dtvec)
%--------------------------------------------------------------------------
% May-15-2013, C. Brandt, San Diego
% Velocimetry2D_Extract_V_Triangle calculates the velocity vector 'vvec' at
% the center position 'cpos' of three pixels at the positions 'povec'.
% Between each pixel pair the time delay 'dtvec' is needed.
% --
% The assumptions for the used equations are that an object passes the 
% three pixels with a constant velocity and phase front of the object is
% straight (not bent), like a infinite wave front. This means the object
% is assumed to be much larger than the covered area of the pixels.
%--------------------------------------------------------------------------
%IN:
% povec: structure with (x,y)-positions p1, p2, p3 of the chosen pixels
% dtvec: structure with the 3 time delays t21, t32, t13
%        (whereas each time delay means: tab=t_b-t_a, and so on)
%OUT:
% cpos: center position of velocity vector
% vvec: 2D velocity vector components (at position posvec)
%--------------------------------------------------------------------------
%EXAMPLE: (velocity vector assumed: v = 86.6ex -50ey)
% povec.p1=[0 0];
% povec.p2=[100 0];
% povec.p3=[50 100];
% dtvec.tab = +0.866;
% dtvec.tbc = -0.933;
% dtvec.tca = +0.067;
% [cpos,vvec] = Velocimetry2D_Extract_V_Triangle(povec,dtvec);
%==========================================================================

% Shorten the variables
tab = dtvec.tab;
tbc = dtvec.tbc;
tca = dtvec.tca;
% THINK ABOUT THE PIXEL DIRECTION!
% povec = [3 1] means x=1, y=3!
ax  = povec.p1(2); ay = povec.p1(1);
bx  = povec.p2(2); by = povec.p2(1);
cx  = povec.p3(2); cy = povec.p3(1);

% OUTPUT: Calculate the center position of the velocity vector
cpos = (povec.p1 + povec.p2 + povec.p3)/3;

% Calculate the velocity in y-direction
% Calculate Prefactor of Vx = F * Vy
F = (by*(tab+tbc)+cy*(tbc+tca)-ay*(tab+2*tbc+tca)) ./ ...
  ( -bx*(tab+tbc)-cx*(tbc+tca)+ax*(tab+2*tbc+tca) );
% Calculate Vy
Vy = ((ax-cx)*F + (ay-cy)) ./ (tca.*((F.^2) + 1));
% Calculate Vx
Vx = F.*Vy;

% OUTPUT: Velocity Vector
vvec = [Vx, Vy];

end