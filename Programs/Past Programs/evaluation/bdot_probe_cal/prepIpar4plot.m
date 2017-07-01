function [Ipa] = prepIpar4plot(X,Y,Ipar,phi_fix,Xcorr,Ycorr)
%function [Ipa] = prepIpar4plot(X,Y,Ipar,phi_fix,Xcorr,Ycorr)
% 1. Load evaluation data of bdot probe with
%    load B_Ipar.mat
% 2. phi123vec
% 3. A = prepIpar4plot(X,Y,Ipar,phi123vec)
% 4. h = pcolor(A.X, A.Y, A.mat);
%
% IN: X,Y: position vectors
%     Ipar: Ipar matrix
%     phi_fix: phase of the plot (0..1)
%     Xcorr,Ycorr: correction of the x- and y-axis (CENTERING)
%OUT:
% EX:

if nargin < 5; Xcorr=0; end
if nargin < 6; Ycorr=0; end

      % INTERPOLATE DATA
        X= interp2(X, 2);
        Y= interp2(Y, 2);
        % small correction of X, Y
          X=X+Xcorr;
          Y=Y+Ycorr;
        Ipar= interp2(Ipar, 2);
      % MAKE X,Y,I_par QUADRATIC
        maxx = max(max(X)); minx = min(min(X));
        maxy = max(max(Y)); miny = min(min(Y));
        limmax = min(abs([ maxx maxy minx miny ]));
        % NEW LIMITS:
        xylim = [-limmax limmax];
        % APPLY xylim to X,Y,I_par
        indx = find( X(1,:)>xylim(1) & X(1,:)<xylim(2) );
        indy = find( Y(:,1)>xylim(1) & Y(:,1)<xylim(2) );
        X = X(indy, indx);
        Y = Y(indy, indx);
        Ipar = Ipar(indy, indx);
      % PREPARE Ipar ACCORDING to phi123vec SHIFTED PLOTS
        tvec = phi_fix;
        phi = phi_fix*2*pi;
          ph_fac= exp(i.*phi);
        Ipar_fix(:, :)= real(Ipar .*ph_fac);
% CREATE OUTPUT
  Ipa.X = X;
  Ipa.Y = Y;  
  Ipa.mat = Ipar_fix;
  Ipa.clim = clim_sym(Ipar_fix);
% HINT FOR PLOT
	disp('continue with pcolor(Ipa.X, Ipa.Y, Ipa.mat)');
end