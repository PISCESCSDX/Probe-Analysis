function [VF] = Velocimetry2D_TriangleMethod(filename, indtime, ... 
  d,crco,pix2r,calcphstd,chk)
%==========================================================================
% function [VF] = Velocimetry2D_TriangleMethod(filename, indtime, ... 
%   d,crco,pix2r)
%--------------------------------------------------------------------------
% May-17-2013, C. Brandt, San Diego
% Velocimetry2D_TriangleMethod calculates the velocity vector field of a
% movie. It uses the whole image size.
% The results are saved to the file 'velo_triangle_d*_t1*-t2*.mat'.
% Last change: Jun-13-2013
%--------------------------------------------------------------------------
%IN:
% filename: string of filename of movie file (cine)
% indtime: indices of time interval
% d: # of pixels between the bottom pixel pair E 1,3,5,7, ...
% crco: structure for cross correlation paramters:
%     .winrat: window length ratio to signal length [0..1]
%     .fend: cut spectrum at frequency fend
% pix2r: conversion factor from pixel to meter (meter = pixel * pix2r)
% calcphstd: calculate standard deviation of time lag; 1 yes, 0 no (faster)
% chk: structure
%    .checked 1: checkplot activated
%    .frequency: plot frequency
%OUT:
% VF: structure containing variables of the velocity field
%   cpos(ver-dim,hor-dim,center position of velocity vector (1:x, 2:y))
%   vvec(ver-dim,hor-dim,freq-vec,velocity vector (1:x, 2:y))
%  vcpsd(ver-dim,hor-dim,freq-vec,CPSD of 1:x 2:y v-component)
%--------------------------------------------------------------------------
%EXAMPLE:
% d = 1;
% crco.n = 5; crco.fend = 20e3;
% pix2r = 1.5e-3;
% indtime = 1:5000;
% VF = Velocimetry2D_TriangleMethod(filename, indtime, d, crco, pix2r);
%==========================================================================

if nargin<7; chk.checked = 0; end
  
% DISPLAY MESSAGE
disp('*** Calculating Velocity Field using Triangle Method (3 pixels) ***')

% Extract movie image informations
cineinfo = cineInfo(filename);
wdth = cineinfo.Width;
hght = cineinfo.Height;

% Define pixels for which the time series will be extracted
ctr = 0;
jvec = 1:wdth;
kvec = 1:hght;
pix = NaN(wdth*hght,2);  % Preallocation
for k=kvec
  for j=jvec
    ctr = ctr+1;
    pix(ctr,:) = [k j];  % !!! j:along width, k:along height
                         % --> pix index runs from left to right, upwards
  end
end

% Load pixel data
[tt,P] = pixel2tt(filename,pix,'checkplot-off',indtime);
fs = 1/(tt(2)-tt(1));

% Parameters for Cross-Correlation calculation
winrat = crco.winrat;
  fend = crco.fend; 
  olap = 0.5;
  kbin = 20;

% Extract the length of the frequency spectrum (for array preallocation)
[~,freq] = fftwindowparameter(size(P,1),winrat,olap,fs,fend);
lfreq = length(freq.cut);

% Pre-Allocate matrices for position and velocity
 cpos = zeros(length(kvec),length(jvec),2);
 vvec = zeros(length(kvec),length(jvec),lfreq,2);
vcpsd = zeros(length(kvec),length(jvec),lfreq);

% LOOP: Go through all pixel-triangles and calculate the time lags
% Determine the center-pixel matrix (positions where the vectors are 
% calculated)
lk = numel(kvec);
% Vertical step pixel distance (height of equal sided triangle)
dpixy = ceil(d*sqrt(3)/2);
jvec2 = 1:wdth-d; kvec2 = 1:hght-dpixy;
figeps(14,14,1, 0,100);
if chk.checked; figeps(14,14,2,40,100); end
% kvec2=54:1:74; jvec2=kvec2;
[grid.horizontal, grid.vertical] = meshgrid(jvec,kvec);
for k=kvec2
  for j=jvec2

    % Determine the triple pixels
    ind1 = (k-1)*lk+1 + (j-1);
    ind2 = ind1 + d;
    ind3 = ((k-1)+dpixy)*lk + j+floor(d/2);
    
    % Extract the 3 signals
    s1 = P(:,ind1);
    s2 = P(:,ind2);
    s3 = P(:,ind3);

  %-------------------------------------------Test Pixel Triangle Positions
% figure(1)
pic = zeros(128,128);
pic( pix(ind1,1), pix(ind1,2) ) = 0.5; % checked 20130612-10h24
pic( pix(ind2,1), pix(ind2,2) ) = 0.7;
pic( pix(ind3,1), pix(ind3,2) ) = 0.9;
pcolor(pic); shading flat
title(['vertical horizontal: ' num2str([k j])])
drawnow
  %------------------------------------------------------------------------
    
  %-----------------------------------------------Test Time Delay Direction
  % dt=1e-5; fs = 1/dt; n=5; olap = 0.5; fend = 20e3; kbin=20;
  % t=(0:1e4)*dt; f=3000;
  % s1 = sin(2*pi*f*t) + 0.01*randn(1,length(t));
  % s2 = sin(2*pi*f*t - pi/2) + 0.01*randn(1,length(t));
  % [A,ph,frq,phstd] = cpsdmean(s1,s2,fs,5,0.5,10e3);
  % [dt1,pa] = Velocimetry2D_Extract_TimeDelay(s1,s2,fs,n,olap,fend,kbin);
  % % s2 after s1: must give a positive time lag 
  % %              of dtau = pi/2 / (2*pi*f) = +8.333e-5
  % %              and indeed: dt1(154) = +8.34e-5
  %------------------------------------------------------------------------

  % Calculate time delays between pixel pairs AB(1,2), BC(2,3) and CD(3,1)
  % Triangle Definition:  C                     ^
  %                      A B  <> horizontal  &  v vertical
    [dt1, pa1] = Velocimetry2D_Extract_TimeDelay(s1,s2, ...
      fs,winrat,olap,fend,kbin,calcphstd);
    [dt2, pa2] = Velocimetry2D_Extract_TimeDelay(s2,s3, ...
      fs,winrat,olap,fend,kbin,calcphstd);
    [dt3, pa3] = Velocimetry2D_Extract_TimeDelay(s3,s1, ...
      fs,winrat,olap,fend,kbin,calcphstd);

    dtvec.tab = dt1;
    dtvec.tbc = dt2;
    dtvec.tca = dt3;
    
  % Coordinates
    povec.p1 = pix(ind1,:) *pix2r;
    povec.p2 = pix(ind2,:) *pix2r;
    povec.p3 = pix(ind3,:) *pix2r;

    disp(num2str([j k]))
    [cp,vv] = Velocimetry2D_Extract_V_Triangle(povec,dtvec); % checked
    i_vertical = k-kvec(1)+1;
    i_horizontal = j-jvec(1)+1;
    cpos(i_vertical,i_horizontal,:)   = cp;
    vvec(i_vertical,i_horizontal,:,:) = vv;
    % Save sum of CPSD-Spectrum of pixel triple
    vcpsd(i_vertical,i_horizontal,:,:) = (pa1.cpsd+pa2.cpsd+pa3.cpsd)/3;

    %%%%%%%%%%%%% CHECK PLOT QUIVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chk.checked
    ind = find(pa1.freq>chk.frequency);
    figure(2);
    vx = vvec(:,:,ind(1),1); vy = vvec(:,:,ind(1),2);
    scale = 2000;
    ind = abs(vx)>2*scale; vx(ind) = 0; vy(ind) = 0;
    ind = abs(vy)>2*scale; vx(ind) = 0; vy(ind) = 0;
    hp = quiver(grid.horizontal,grid.vertical, 5*vx/scale,5*vy/scale);
    axis square
    set(hp,'autoscale','off')
    set(gca,'xlim',[jvec2(1) jvec2(end)],'ylim',[kvec2(1) kvec2(end)])
    drawnow
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
end

t1 = mkstring('','0',indtime(  1),9999,'');
t2 = mkstring('','0',indtime(end),9999,'');
savefn = [filename(1:end-5) '_velo_triangle_d' num2str(d) '_t' t1 '-' t2 ...
  '_winrat' sprintf('%.2f',winrat) '.mat'];
save(savefn,'grid','cpos','vvec','vcpsd','pa1','pix','jvec','pix2r')

% Store variables in output variable
VF.vec   = vvec;
VF.pos1  = cpos;
VF.vcpsd = vcpsd;
VF.crco  = pa1;
VF.pix   = pix;
VF.pix2r = pix2r;

end