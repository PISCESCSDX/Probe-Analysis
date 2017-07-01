function [VF] = Velocimetry2D_SquareMethod(filename, d, crco, pix2r)
%==========================================================================
%function [VF] = Velocimetry2D_SquareMethod(filename, d, crco, pix2r)
%--------------------------------------------------------------------------
% May-17-2013, C. Brandt, San Diego
% Velocimetry2D_SquareMethod calculates the velocity vector field of a
% movie. It uses the whole image size.
% The results are saved to the file 'velo_square.mat'.
%--------------------------------------------------------------------------
%IN:
% filename: string of filename of movie file (cine)
% d: # of pixels between each pixel pair E 1,3,5,7, ...
% crco: structure for cross correlation paramters:
%     .n: number of windows
%     .fend: cut spectrum at frequency fend
% pix2r: conversion factor from pixel to meter (meter = pixel * pix2r)
%OUT:
% VF: structure containing variables of the velocity field
%--------------------------------------------------------------------------
%EXAMPLE:
% filename = '15172.cine';
% d = 1;
% crco.n = 5; crco.fend = 20e3;
% pix2r = 1.5e-3;
% [VF] = Velocimetry2D_SquareMethod(filename, d, crco, pix2r)
%==========================================================================


% DISPLAY MESSAGE
disp('*** Calculating Velocity Field using Square Method (4 pixels) ***')

% Extract movie image informations
cineinfo = cineInfo(filename);
wdth = cineinfo.Width;
hght = cineinfo.Height;

% Define pixels for which the time series need to be extracted
ctr = 0;
jvec = 1:hght;
kvec = 1:wdth;
pix = NaN(wdth*hght,2);   % Preallocation
for j=jvec
  for k=kvec
    ctr = ctr+1;
    pix(ctr,:) = [j k];
  end
end

% Load pixel data
chk = 'checkplot-off';
% *** [tt,P] = pixel2tt(filename,pix,chk);
load deleteme.mat
fs = 1/(tt(2)-tt(1));

% Parameters for Cross-Correlation calculation
   n = crco.n;
fend = crco.fend; 
olap = 0.5;
kbin = 20;

% Pre-Allocate matrices for position and velocity
cpos1 = zeros(length(jvec),length(kvec),2);
vvec1 = zeros(length(jvec),length(kvec),250,2);
cpos2 = zeros(length(jvec),length(kvec),2);
vvec2 = zeros(length(jvec),length(kvec),250,2);

% Go through all pixel-squares and calculate the time lags
dh= (d+1)/2;

for j=40:41
  disp(num2str([j length(jvec)-1]));
  for k=60:61

% ORIGINAL:
% for j=1+d:hght-d
%   disp(num2str([j length(jvec)-1]));
%   for k=1+d:wdth-d

%     ind1 = (j-1)*numel(jvec) + (k-1); % pixel 1: (j  ,k-1)
%     ind2 = (j-1)*numel(jvec) + (k+1); % pixel 2: (j  ,k+1)
%     ind3 = (j-2)*numel(jvec) +  k;    % pixel 3: (j-1,k  )
%     ind4 = (j  )*numel(jvec) +  k;    % pixel 4: (j+1,k  )
    ind1 = (j-1   )*numel(jvec) + (k-dh); % pixel 1: (j  ,k-dh)
    ind2 = (j-1   )*numel(jvec) + (k+dh); % pixel 2: (j  ,k+dh)
    ind3 = (j-1-dh)*numel(jvec) +  k;     % pixel 3: (j-dh,k  )
    ind4 = (j-1+dh)*numel(jvec) +  k;     % pixel 4: (j+dh,k  )

    % Extract the 4 signals
    s1 = P(:,ind1);
    s2 = P(:,ind2);
    s3 = P(:,ind3);
    s4 = P(:,ind4);
    
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
  
  % Calculate time delays between pixel pairs (1,2) and (3,4)
    [dt1, ~] = Velocimetry2D_Extract_TimeDelay(s1,s2,fs,n,olap,fend,kbin);
    [dt2,pa] = Velocimetry2D_Extract_TimeDelay(s3,s4,fs,n,olap,fend,kbin);

    dtvec.tab = dt1;
    dtvec.tcd = dt2;

% *** Check coordinates
    povec.p1 = pix(ind1,:) *pix2r;
    povec.p2 = pix(ind2,:) *pix2r;
    povec.p3 = pix(ind3,:) *pix2r;
    povec.p4 = pix(ind4,:) *pix2r;

    [cp,vv] = Velocimetry2D_Extract_V_Square(povec,dtvec);
    
    cpos1(j,k,:)   = cp;
    vvec1(j,k,:,:) = vv;

    % % Diagonal elements
    % ind1 = (j  )*numel(jvec) + (k-1); % 1: j+1 k-1
    % ind2 = (j-2)*numel(jvec) + (k+1); % 2: j-1 k+1
    % ind3 = (j  )*numel(jvec) + (k+1); % 3: j+1 k+1
    % ind4 = (j-2)*numel(jvec) + (k-1); % 4: j-1 k-1
    % s1 = P(:,ind1);
    % s2 = P(:,ind2);
    % s3 = P(:,ind3);
    % s4 = P(:,ind4);
    % 
    % [dt1, ~] = Velocimetry2D_Extract_TimeDelay(s1,s2,fs,n,olap,fend,k);
    % [dt2,pa] = Velocimetry2D_Extract_TimeDelay(s3,s4,fs,n,olap,fend,k);
    % 
    % dtvec.tab = dt1;
    % dtvec.tcd = dt2;
    % 
    % povec.p1 = pix(ind1,:);
    % povec.p2 = pix(ind2,:);
    % povec.p3 = pix(ind3,:);
    % povec.p4 = pix(ind4,:);

    [cp,vv] = Velocimetry2D_Extract_V_Square(povec,dtvec);
    
    cpos2(j,k,:)   = cp*pix2r;
    vvec2(j,k,:,:) = vv*pix2r;

  end
  save('velo_square_corr1.mat','vvec1','vvec2','cpos1','cpos2',...
    'pa','pix', 'jvec','pix2r')
end

% Store variables in output variable
VF.vec1 = vvec1; VF.pos1 = cpos1;
VF.vec2 = vvec2; VF.pos2 = cpos2;
VF.crco = pa;
VF.pix = pix;
VF.pix2r = pix2r;

end