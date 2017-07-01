% Oct-15-2013, C. Brandt, San Diego 
% m-file for replacing mkcalc.m to pre-process video data
%
% Manual:
% 1. Check all ### INPUTs in this file and Run it
% 2. Output is a phase data matrix and a folder with fft decomposition file
%    for each frame

a = dir('*.cine');

for i_files=1:length(a)
% ### INPUT: filenames of movie and statistic file
mov.path  = '/home/csdx/20131011_CSDX/camera';
filebase = a(i_files).name(1:end-5);
disp(['*** FFT decomposing movie ' a(i_files).name])
mov.fn    = [filebase '.cine'];

if ~isdir(filebase);
  mkdir(filebase);
end
fftdec.savebase = filebase;

%%%%%%%%%%%%%%%%%%%%%%% Use this for calculating the pixel to mm conversion
% rvec = [48 40 32 24];
% load 19479_statistics
% pcolor( movstat.avg ); shading flat
% inp = input('Enter current position of longest probe (pix): ');
% rcam(1) = inp;
% load 19480_statistics
% pcolor( movstat.avg ); shading flat
% inp = input('Enter current position of longest probe (pix): ');
% rcam(2) = inp;
% load 19481_statistics
% pcolor( movstat.avg ); shading flat
% inp = input('Enter current position of longest probe (pix): ');
% rcam(3) = inp;
% load 19482_statistics
% pcolor( movstat.avg ); shading flat
% inp = input('Enter current position of longest probe (pix): ');
% rcam(4) = inp;
% plot(rvec, rcam, 'o-')
% xlabel('real position (mm)')
% ylabel('camera position (pixel)')
% p = polyfit(rvec,rcam,1);
% pix2r = 1/p(1);
% disp('The pixel to radial conversion is: ')
% disp(num2str(pix2r)) % 2.1622 mm/pixel
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% ### Input: If needed activate the statistic calculation
% fct_a_moviestatistic(filebase);

mov.pathfn = [mov.path '/' mov.fn];
mov.avg.fn = [filebase '_statistics.mat'];
mov.avg.pathfn = [mov.path '/' mov.avg.fn];

% ### INPUT: pixel to millimeter conversion (one pixel is ... mm)
fftdec.pix2r = 2.1622;

% Load statistic file
load(mov.avg.pathfn);


%==========================================================================
% Find the optimal center
%--------------------------------------------------------------------------
% # INPUT: pixel range to find center
rd = 10;
% # INPUT: Resolution of azimuthal angle
resphi = 128;
% load averaged image picture
matavg = movstat.avg;
% Find start value: maximum of mat
[v1, i1] = max(matavg);
[ ~, colmax] = max(v1);
rowmax = i1(colmax);
vrow = -rd:1:rd;  vcol = vrow;
int = zeros(numel(vrow),numel(vcol));
for irow = 1:length(vrow)
  for icol = 1:length(vcol)
    cptest = [rowmax+vrow(irow) colmax+vcol(icol)];
    [rvec, mavg] = Matrix2dAzimuthAvg(matavg, cptest, resphi);
    y = std(mavg,1)';
    % Integration of standard deviation is a measure of deviation:
    int(irow,icol) = int_discrete(rvec,y);
  end
end
% Find indices of local minimum of integration matrix 'int'
[v1, i1] = min(int);
[ ~, colmin] = min(v1);
rowmin = i1(colmin);
% Store new center position to variable
cp = [rowmax+vrow(rowmin) colmax+vcol(colmin)];
% Show Testplot
matavg = matavg./matmax(matavg);
matavg(cp(1),cp(2)) = 2;
pcolor(matavg); shading flat; axis square; set(gca,'clim',[0 1.2])
% *** inp = input('Is the center ok? (yes: Enter, no: 0)');
inp = [];
if ~isempty(inp); return; end
%==========================================================================


%==========================================================================
% Find the cut range of the rawdata frames
%--------------------------------------------------------------------------
% # Input: filename of movie
% Save info.mat with video information
fftdec.info = cineInfo(mov.pathfn);
fftdec.moviefile = mov.pathfn;
fftdec.statfile  = mov.avg.pathfn;
% Read first image to determine image size for pre-allocation
pic = double(cineRead(mov.pathfn,1));
disp(['Size of pic: ' num2str(size(pic))])
% Get dimensions of matrix
sizever = size(pic,1);
sizehor = size(pic,2);
% --- Determine lowest distance from center position to boundary of data --
% (This is the maximum of the radial vector.)
% Distances from center position:
 distLeft = cp(2) - 1;
distRight = sizehor - cp(2);
   distUp = sizever - cp(1);
 distDown = cp(1) - 1;
% Minimal distance:
mindist = matmin([distLeft distRight distUp distDown]);
% # INPUT: left-right range and up-down range
cut.ud = cp(1)-mindist:cp(1)+mindist;
cut.lr = cp(2)-mindist:cp(2)+mindist;
% plot image and check cut range
pcolor(pic(cut.ud,cut.lr)); shading flat; axis square
% *** inp = input('Is the cutting range ok?  (y: enter, n: 0)  ');
inp=[];
if ~isempty(inp); return; end
close all

pixpos.cutrange.lr = [cut.lr(1) cut.lr(end)];
pixpos.cutrange.ud = [cut.ud(1) cut.ud(end)];
pixpos.probecenter.lr = cp(2);
pixpos.probecenter.ud = cp(1);
save(mov.avg.fn, '-append','pixpos')
%==========================================================================


%==========================================================================
% FFT decomposition of single light fluctuation images
%--------------------------------------------------------------------------
% Do not use improcess8, instead just use camerastatistics and mksum!
% Define FFT data range
fftdec.startframe = 1;
fftdec.endframe   = fftdec.info.NumFrames;
% Define center pixel
fftdec.cp = cp;
% Define cut range
fftdec.cut.ud = cut.ud;
fftdec.cut.lr = cut.lr;
% Number of azimuthal pixel probes
fftdec.resphi = resphi;
% ### Input: Details of radial decomposition
fftdec.dr = 1;
fftdec.rmin = 1;
fftdec.rmax = mindist;
% ### Input: mode vector
fftdec.m = 0:5;

% Transfer light fluctuation quantities to variable 'fftdec'
fftdec.lightfluc.std2 = movstat.lightfluc.std2;
fftdec.lightfluc.amp = movstat.lightfluc.amp;
fftdec.movieavg = movstat.avg;

% Do the FFT Decomposition
FFTDecompAzimuth(fftdec);
%==========================================================================
end
return


%==========================================================================
% Play Cine Video
%--------------------------------------------------------------------------
% Read the first frame
fnum = 1;
cdata = cineRead(fn,fnum);

figeps(15,8,1);
for i=1:1000
  cdata = double(cineRead(fn, i));
  clf;
  pcolor(cdata); axis equal
  set(gca, 'clim', [1 3000])
  shading flat
  pause(0.1)
end
%==========================================================================