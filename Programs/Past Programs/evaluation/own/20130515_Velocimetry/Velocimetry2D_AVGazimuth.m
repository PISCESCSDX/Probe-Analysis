function [rvec, vazavg] = Velocimetry2D_AVGazimuth(filename, cp)
%==========================================================================
%function [rvec, vazavg] = Velocimetry2D_AVGazimuth(filename, cp)
%--------------------------------------------------------------------------
% Jun-04-2013, C. Brandt, San Diego
% This m-file is supposed to calculate the azimuthally averaged velocity.
% The center pixel needs to be defined.
%--------------------------------------------------------------------------
%INPUT
%  filename: string of filename containing the velocimetry data
%  cp: coordinates of center position [centerrow centercolumn]
%OUTPUT
%  rvec: radial pixel vector
%  vazavg: structure containing:
%    .vrad: radial velocity (negative: inward)
%    .vtheta (positive: anti-clockwise)
%    .freq frequency vector
%--------------------------------------------------------------------------
%EXAMPLE
% filename = 'velo_triangle_d1_t0001-5000_winrat0.20.mat';
% cp = [66 66];
% [rvec, vvec] = Velocimetry2D_AVGazimuth(filename, cp);
%==========================================================================

% Load velocimetry data
load(filename);

% Get size of vector array
% !!! As defined in function Velocimetry2D_Extract_V_Triangle
% vvec 1st dimension: y-position
% vvec 2nd dimension: x-position
% vvec 3rd dimension: frequency
% vvec 4th dimension: velocity vector x,y-components: 1:x, 2:y
% (checked the above on June-11-2013, Christian)
% ---
% According to function Velocimetry2D_Extract_V_Triangle
% Vx really means velocity in horizontal direction in a picture.
% Because: "Point(3,1) means x=1, y=3"  !!!
sizex = size(vvec,2); %#ok<NODEF>
sizey = size(vvec,1);
 Nfre = size(vvec,3);

% Determine lowest distance from center position to boundary of data
% (This is the maximum of the radial vector.)

% Distances from center position:
 distLeft = cp(2) - 1;
distRight = sizex - cp(2);
   distUp = sizey - cp(1);
 distDown = cp(1) - 1;

% Minimal distance:
mindist = matmin([distLeft distRight distUp distDown]);

% --- Create radial-azimuthal array (zero in the center) ---
   N = 64;
 phi = ((1:N+1)'/N)*2*pi;
rvec = (0:1:mindist)';
% Meshgrid
[Mrad,Mphi] = meshgrid(rvec,phi);

% --- Calculate interpolated points at azimuthal array ---
% Help array = azimuthal array (zero center) + cp
xi = Mrad.*cos(Mphi) + cp(1);  % xi means horizontal
yi = Mrad.*sin(Mphi) + cp(2);  % yi means vertical

% Old coordinates
x = (1:1:sizex)';  % x means horizontal
y = (1:1:sizey)';  % y means vertical
[Y,X] = meshgrid(y,x);

% Calculate interpolated radial-azimuthal array
vhor = zeros(size(xi,1),size(xi,2),Nfre);
vver = vhor;
for i=1:Nfre % *** Maybe error in this loop ??? X Y xi yi correct?
%   vhor(:,:,i) = interp2(X,Y,vvec(:,:,i,1),xi,yi);
%   vver(:,:,i) = interp2(X,Y,vvec(:,:,i,2),xi,yi);
  vhor(:,:,i) = interp2(Y,X,vvec(:,:,i,1),yi,xi); % because vvec(:,,,) vertic
  vver(:,:,i) = interp2(Y,X,vvec(:,:,i,2),yi,xi); % vvec(,1,,) horizontal
end

% % Test plots
% pcolor(xi,yi,vhor(:,:,72)); set(gca,'clim',500*[-1 1])
% pcolor(xi,yi,vver(:,:,72)); set(gca,'clim',500*[-1 1])

% --- Calculate azimuthally averaged velocity ---
for i=1:Nfre
  % Coordinate system used for vvec:
  %  vvec(:,-,-,-) is Y
  %  vvec(-,:,-,-) is X
  %  phi is increasing from x to y
  %  HOWEVER, velocity is defined as vvec(-,-,-,1) horizontal
  %  HOWEVER, velocity is defined as vvec(-,-,-,2) vertical
  % x
  % ^
  % |
  % |
  % |
  % |
  % |
  % |
  % +------------------>y
  vabs = sqrt( (vhor(:,:,i).^2) + (vver(:,:,i).^2) );
  alpha = Mphi + atan2( vver(:,:,i), vhor(:,:,i) );
  urad(:,:,i) = vabs .* sin(alpha);
  uphi(:,:,i) = vabs .* cos(alpha);
end

vazavg.vrad     = urad;
vazavg.vtheta = uphi;
vazavg.freq = pa1.freq;

end