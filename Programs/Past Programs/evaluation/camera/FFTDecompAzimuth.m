function FFTDecompAzimuth(inp)
%==========================================================================
%function FFTDecompAzimuth(inp)
%--------------------------------------------------------------------------
% Read images from video file, extract azimuthal signals for a given set of
% radii and calculate the inverse Fourier transformation to decompose the
% mode numbers.
%
% (C) 10.09.2013, C. Brandt, San Diego
% (C) 07.05.2013, C. Brandt, San Diego
% - use ffteigenmode.m for better FFT decomposition
%---------- Normalization Info --------------------------------------------
% (1) The average fluctuation amplitude is calculated
%     'lightfluc.amp' (in camera intensity units (CIU)).
%     here: averaged: AvgAmpLfluc
% 
% (2) From that the average energy in the 2D azimuthal fluctuations is 
%     calculated: here: AvgEnLflucArea
%
% (3) NEW: The pictures from the cine files are loaded and the average is
%     subtracted.
%
% (4) Shift is substracted -> pure fluctuation data in Camera Intensity
%     units (CIU).
%
% (5) FFT decomposition and reconstruction.
%
% (6) Normalization of the 2D reconstruction result to 'AvgAmpLfluc'.
%
% (7) Radial mode energies are calculated by A^2
%     Normalizeation by: ***
%------------------------------------------------------------INPUT & OUTPUT
% inp.startframe, inp.endframe: start and end frame
% inp.cp:   definition of the radial position r=0, [vertical horizontal]
% inp.pix2r: conversion factor from pixel to millimeter
% inp.rmin: minimum radius
% inp.rmax: maximum radius
% inp.dr:   step radius
% inp.m:    mode vector (modes to decompose)
%==========================================================================


% sample frequency
fs = inp.info.frameRate;

%==========================================================================
% # INPUT
%==========================================================================
ilim = inp.startframe:1:inp.endframe;
% number of azimuthal "pixel probes"
num = inp.resphi;
% center pixel
cp = inp.cp;
% Calculate from pixel to radial units (mm)
pix2r = inp.pix2r;
% Define radial vector
rmin = inp.rmin;
rmax = inp.rmax;
dr   = inp.dr;
rvec = rmin:dr:rmax;
% mode vector
mvec = inp.m;

% Initialize variables
ctr = 0;

%--------------------------------------------------------------------------
% Average amplitudes of measured light fluctuations (std2 of each image)
% (cam-intensity units: CIU)
AvgAmpLfluc = mean(inp.lightfluc.std2);
% Average Energy in whole azimuthal area (in CIU*pixel^2)
% Use same definition as for the calculation of the energies!
AvgEnLflucArea = (AvgAmpLfluc/sqrt(2) * (2*pi*rvec*pix2r) ).^2 ;
% Average maximal amplitudes
AvgMaxLfluc = mean( inp.lightfluc.amp );
%--------------------------------------------------------------------------

%==========================================================================
% Temporal for-loop ('ilim' images)
%==========================================================================
En_m(length(rvec)) = 0;
En_cam(length(rvec)) = 0;
invm{length(mvec)} = zeros(length(mvec),num);
En_mr{length(mvec)} = zeros(1,length(rvec));
for i=ilim
  ctr = ctr+1;
  disp(num2str(ctr))
    
  %------------------------------------------
  % READ MEASUREMENT IMAGE
  pic = double(cineRead(inp.moviefile,i));
  %------------------------------------------
  % Subtract average
  pic = pic - inp.movieavg;
  %------------------------------------------
  
  % Preallocation
  c_r = 0;
  X = zeros(num,numel(rvec));
  Y = zeros(num,numel(rvec));
  for r=rvec
    c_r = c_r+1;
    
    % Extract artificial probe array, transpose for FFT
    %--------------------------------------------------
    % (Info about direction in CAMEXTRACTCOU.m
    %  17.09.2013
    %  - extract. azimut. array in mathemat. pos. sense: counter-clockwise
    %  -> phase increase: in (phi,A)-plot wave propagates to left side,
    %                     in camera view clockwise
    %  -> phase decrease: in (phi,A)-plot wave propagates to right side,
    %                     in camera view counter-clockwise
    [phi, Avec, ~] = camextractcou(pic, r, num, cp);
    Avec = Avec';  phi = phi';
    
    % Calculate Fourier transformation (phi/2 to obtain integer m)
    [~, amp, pha] = ffteigenmode(phi/2, Avec, 10);
    phase.amp{c_r}(i,:) = amp;
    phase.pha{c_r}(i,:) = pha;
    
    % Calculate Inverse Fourier transformation for each mode in mvec
    for im = 1:length(mvec)
      invm{im}(c_r,:) = amp(im)*cos(pi*mvec(im)*phi+pha(im));
    end
    
    % Calculate radial light profiles
    Esum = 0;
    for im = 1:length(mvec)
      % En_mr: radial dependent mode energy (unit: (amp*length)^2 )
      En_mr{im}(c_r)  = ( (amp(im)/sqrt(2)) *2*pi*(r*pix2r) )^2;
      Esum = Esum + En_mr{im}(c_r);
    end
    
    % Energy reconstructed from mode = mvec
    En_m(c_r) = Esum;
    % Energy in camera measurement (all FFT-modes: m=1-32)
    En_cam(c_r) = sum( ((amp/sqrt(2)) * 2*pi*(r*pix2r)).^2 );

    % Return to xy coordinates and plot with pcolor
    % (sin and cos are the other way around, because of the direction
    %  definition of camextractcou.m going clockwise from 12h00)
    X(:,c_r) = pix2r* r*sin(phi'*pi);
    Y(:,c_r) = pix2r* r*cos(phi'*pi);
  end



%==========================================================================
% MODE PLOTS
%==========================================================================
% Cut image up down and left right range
pic = pic(inp.cut.ud, inp.cut.lr);

% Calculate position vectors (center=0,0)
nX = pix2r*( inp.cut.lr-cp(2) );
nY = pix2r*( inp.cut.ud-cp(1) );

% Store camera raw image
mode.cam.x = nX';
mode.cam.y = nY';
mode.cam.pic = pic/AvgMaxLfluc;
mode.fs  = fs;
mode.i_t = ctr;

% Superimposed plot: Sum of all mode numbers taken into account
te = invm{1}';
for g=2:length(mvec)
  te = te + invm{g}';
end


%==========================================================================
% Store mode decomposed plots
%--------------------------------------------------------------------------
% Add one column to make 1st and last element periodic
x_show = Y; x_show(size(Y,1)+1, :) = x_show(1,:);
y_show = X; y_show(size(X,1)+1, :) = y_show(1,:);
% For length of mvec and +1 for storing the superimposed plot at the end
for im=1:length(mvec)+1
  if im<length(mvec)+1
    M_show = invm{im}'/AvgMaxLfluc;
    M_show(size(M_show,1)+1,:) = M_show(1,:);
  else
    M_show = te/AvgMaxLfluc;
    M_show(size(M_show,1)+1,:) = M_show(1,:);
  end
  %------------------------------------------------------------------------
  % INFO for the variable 'mode':
  %   M_show is normalized to AvgMaxLfluc, i.e., to the average maximum
  %   light fluctuation amplitudes (done in FCT_A_MOVIESTATISTIC.m).
  %   This has nothing to do with the degree of fluctuation (%)!
  mode.m{im}.x   = x_show;
  mode.m{im}.y   = y_show;
  mode.m{im}.pic = M_show;
  mode.mvec = mvec;
  %------------------------------------------------------------------------
  
  
  %==========================================================================
  % RADIAL LIGHT INTENSITY PROFILES
  %==========================================================================
  radprof.r      = rvec*pix2r;
  radprof.En_m   = En_m./AvgEnLflucArea;
  radprof.En_cam = En_cam./AvgEnLflucArea;

  % figure
  % hold on
  % plot(radprof.r, radprof.En_m,   'b')
  % plot(radprof.r, radprof.En_cam, 'r')
  % hold off
  % return
  
  if im<length(mode.mvec)+1
    radprof.En_mr{im} = En_mr{im}./AvgEnLflucArea;
    tt.en{im}(i) = sum(radprof.En_mr{im});
  end

end

savebase = [inp.savebase '/mkfftdec_'];
fn = mkstring(savebase, '0', i, inp.info.NumFrames, '.mat');
save(fn, 'mode', 'radprof');

end % for-loop: i=ilim


% Save Energy, Amplitude Variables
%==========================================================================
tt.tvec = 1/fs * ((1:length(ilim))'-1); %#ok<STRNU>
save([inp.savebase '/modeenergy.mat'],'tt')


% Save Phase Data Variables
%==========================================================================
% Phase difference between mi and mj
c_r = 0;
for r=rvec
c_r = c_r+1;
ph.phd12{c_r}=(unwrap(phase.pha{c_r}(:,2))-unwrap(phase.pha{c_r}(:,3)))/(2*pi);
ph.phd13{c_r}=(unwrap(phase.pha{c_r}(:,2))-unwrap(phase.pha{c_r}(:,4)))/(2*pi);
ph.phd14{c_r}=(unwrap(phase.pha{c_r}(:,2))-unwrap(phase.pha{c_r}(:,5)))/(2*pi);
ph.phd15{c_r}=(unwrap(phase.pha{c_r}(:,2))-unwrap(phase.pha{c_r}(:,6)))/(2*pi);
ph.phd23{c_r}=(unwrap(phase.pha{c_r}(:,3))-unwrap(phase.pha{c_r}(:,4)))/(2*pi);
ph.phd24{c_r}=(unwrap(phase.pha{c_r}(:,3))-unwrap(phase.pha{c_r}(:,5)))/(2*pi);
ph.phd25{c_r}=(unwrap(phase.pha{c_r}(:,3))-unwrap(phase.pha{c_r}(:,6)))/(2*pi);
ph.phd34{c_r}=(unwrap(phase.pha{c_r}(:,4))-unwrap(phase.pha{c_r}(:,5)))/(2*pi);
ph.phd35{c_r}=(unwrap(phase.pha{c_r}(:,4))-unwrap(phase.pha{c_r}(:,6)))/(2*pi);
ph.phd45{c_r}=(unwrap(phase.pha{c_r}(:,5))-unwrap(phase.pha{c_r}(:,6)))/(2*pi);
end
ph.tvec = (ilim-ilim(1))/fs*1e3;
ph.rvec = rvec;
save([inp.savebase '/phasedata.mat'], 'phase', 'ph');

% % save variables
% mode.info = 'cam: camera pictures, m: mode data (mvec)';
% radprof.info = 'r in cm, En_m: inverse fft (mvec), En_mr: radial mode energy, En_cam: camera pic energy';
% save('mkfftdec.mat','mode', 'radprof','ilim');
end