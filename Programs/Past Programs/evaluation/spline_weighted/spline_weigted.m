% TEST VECTOR
x = 1:10;
y = [1 3 6 4 12 17 21 31 29 40];

plot(x,y, 'o')

% NOW SPLINE IN THIS WAY THAT THE 2 BAD POINTS ARE WEIGHTED LESS THAT REST
% BAD POINTS: x,y=4,4;  x,y=8,31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEIGHTED SPLINE INTERPOLATION (RED)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # INPUT: WEIGHT SPECIAL DATA POINTS
  % WEIGHT OF GOOD POINTS = 100
  w=ones(1,length(x))*100;
  % REDUCE WEIGHT OF BAD POINTS TO 1
  w(4) = 1;
  w(8) = 1;
% SPLINE INTERPOLATION WITH SPAPS (GOOD FOR NOISE REDUCTION)
sp = spaps(x, y, 0200, w, 3);
% CREATE PLOT SPLINE DATA
SR = fnplt(sp);

hold on; plot(SR(1,:), SR(2,:), 'r-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOT WEIGHTED SPLINE INTERPOLATION (BLUE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp = spaps(x, y, 1);
% PLOT SPLINE
SR = fnplt(sp);
hold on; plot(SR(1,:), SR(2,:), 'b-')
