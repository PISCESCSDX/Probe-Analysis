function [PolyIsat, voltlim, ivoltlim, Vf] = ...
  subPlasmaIsatPoly(data,IsRg)
%==========================================================================
%function [PolyIsat, voltlim, ivoltlim, iVf] = subPlasmaIsatPoly(data)
%--------------------------------------------------------------------------
% SUBPLASMAISATPOLY calculates the linear approximation of the ion
% saturation current in any Langmuir current-voltage-characteristic.
% If no voltage limits are given, the limits are chosen automatically.
%--------------------------------------------------------------------------
% IN: data.voltage
%     data.current
%     data.voltlim: [minimum maximum]
%     IsRg: length of ion saturation current range (in Volt)
%OUT: PolyIsat: polynomial of first order of the ion saturation current
%     voltlim: lower and upper voltage limit for the ion saturation current
%     ivoltlim: indices of voltlim
%     iVf: index of the found floating potential
%--------------------------------------------------------------------------
% EX: [PolyIsat, voltlim, ivoltlim, iVf] = subPlasmaIsatPoly(data);
%            y = polyval(PolyIsat, x);
% 11.04.2012-10:21 C. Brandt, San Diego, new I_i,sat algorithm
% 05.04.2012-10:21 C. Brandt, San Diego, added iVf
% 04.04.2012-17:22 C. Brandt, San Diego
%==========================================================================

if nargin < 2
  IsRg = input('ion saturation current range (in Volt): ');
end

x = data.voltage;
dx = (x(end)-x(1)) / numel(x);
y = data.current;

% Define Left hand side Limit:
indcutleft = 1; % changed: 20130212

% Decide which range is the ion saturation current
indcutright = indcutleft + round( IsRg/dx );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot for LaTeX Documentation about the probe m-files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figeps(12,8,1,0,100); clf
% axes('position', [0.15 0.15 0.80 0.80])
% hold on
% plot(x,      y*1e3, 'o-', 'Color', 0.8*[1 1 1])
% plot(x,abs( filtmooth(y, floor(30/dx) ) )*1e3, 'r-')
% hold off
% set(gca, 'xlim', [x(1) x(numel(x))], 'ytick', -2:2)
% ylim = get(gca, 'ylim');
% line(x(indcutright)*[1 1], ylim, 'Color', 'k')
% mkplotnice('voltage (V)', 'I (mA)', 12, '-25', '-35');
% dy = -20;
% puttextonplot(gca, [0 1], 5, -15+0*dy, 'raw data', 0, 10, 0.8*[1 1 1]);
% puttextonplot(gca, [0 1], 5, -15+1*dy, 'y_{filt30}', 0, 10, 'r');
% print('-depsc2', 'subPlasmaIsatPoly.eps')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate linear regression of ion saturation part
% p = polyfit(x(indcutleft:indcutright), yfilt10(indcutleft:indcutright), 1);
p = polyfit(x(indcutleft:indcutright), y(indcutleft:indcutright), 1);
% Calculate linear regression of the whole characteristic
q = polyfit(x, y, 1);

% Exclude increasing ion saturation currents beyond a certain limit
% p at p(1)/q(1) < -0.05 seems to be too large(!)
if p(1)/q(1) > -2%-0.05
  PolyIsat = p;
  voltlim = [x(indcutleft) x(indcutright)];
  ivoltlim = [indcutleft indcutright];
else
  disp('  subPlasmaIsatPoly_v2: No Output: I_i,sat too strongly increasing!')
  PolyIsat = [];
  voltlim  = [];
  ivoltlim = [];
  Vf       = [];
  return
end

% Find Floating Potential
[~, Vf.i] = min( smoothBrochard(abs(y),4) );
Vf.V = x(Vf.i);

% Exclusion Tests: I-V characteristic needs to be asymmetric enough
%--------------------------------------------------------------------------
% Langmuir characteristic: Some Simple Symmetry Tests
% Average of left side must be larger than average of right side
  ye = y - polyval(p,x);

i_le = 1:round(length(x)/2);
i_ri = round(length(x)/2)+1:length(x);

ye_avg_le = mean(ye(i_le));
ye_avg_ri = mean(ye(i_ri));

% ye_std_le = std(ye(i_le));
% ye_std_ri = std(ye(i_ri));

ye_std_le = mean(sqrt(ye(i_le).^2));
ye_std_ri = mean(sqrt(ye(i_ri).^2));

if ye_avg_le <= ye_avg_ri
  disp('  subPlasmaIsatPoly_v2: No Output: Found I_e lower on left half.')
  PolyIsat = [];
  voltlim  = [];
  ivoltlim = [];
  Vf       = [];
elseif ye_avg_ri >= 0
  disp('  subPlasmaIsatPoly_v2: No Output: Positive I_e on right half.')
  PolyIsat = [];
  voltlim  = [];
  ivoltlim = [];
  Vf       = [];
elseif ye_std_ri/ye_std_le < 0.5          % found empiric -> may be wrong
  disp('  subPlasmaIsatPoly_v2: No Output: Std(I_e) too low on right half.')
  PolyIsat = [];
  voltlim  = [];
  ivoltlim = [];
  Vf       = [];
elseif abs(ye_avg_ri)/abs(ye_avg_le) < 1  % found empiric -> may be wrong
  disp('  subPlasmaIsatPoly_v2: No Output: Avg(I_e) too low on right half.')
  PolyIsat = [];
  voltlim  = [];
  ivoltlim = [];
  Vf       = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot for LaTeX Documentation about the probe m-files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figeps(12,8,1,0,100); clf
% axes('position', [0.15 0.15 0.80 0.80])
% hold on
% plot(x,1e3*y,'-o','Color', 0.8*[1 1 1])
% plot(x(indcutleft:indcutright),1e3*y(indcutleft:indcutright),'b-o')
% plot(x,1e3*polyval(p,x),'r')
% hold off
% set(gca, 'xlim', [x(1) x(numel(x))], 'ytick', -2:2)
% mkplotnice('voltage (V)', 'I (mA)', 12, '-25', '-35');
% print('-depsc2', 'subPlasmaIsatPoly_fit.eps')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%XXXXXXXXXXXXXXXXXXXXXXXXXX --- CHECK PLOTS ---  XXXXXXXXXXXXXXXXXXXXXXXXXX
% disp(['avg le ri: ' sprintf('%.1e',ye_avg_le) ...
%   ' ' sprintf('%.1e',ye_avg_ri)])
% disp(['avg ri/le: ' sprintf('%.1e',ye_avg_ri/ye_avg_le) ])
% 
% disp(['std le ri: ' sprintf('%.1e',ye_std_le) ...
%   ' ' sprintf('%.1e',ye_std_ri)])
% 
% disp(['std ri/le: ' sprintf('%.1e',ye_std_ri/ye_std_le) ])
% 
% % --->>> Check Plot
% figure(1); clf
% subplot(2,1,1)
% hold on
% plot(y,'ko')
% plot(indcutleft:indcutright,y(indcutleft:indcutright),'bo')
% plot(polyval(p,x),'r')
% hold off
% 
% if ~isempty(Vf)
%   puttextonplot(gca, [1 1], -100, -15, ['V_f=' num2str(Vf.V)],0,12,'k');
%   puttextonplot(gca, [1 1], -100, -30, ['i_f=' num2str(Vf.i)],0,12,'k');
% else
%   PolyIsat = [];
% end
%  
% subplot(2,1,2)
% plot(ye)
% %---<<<
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

end