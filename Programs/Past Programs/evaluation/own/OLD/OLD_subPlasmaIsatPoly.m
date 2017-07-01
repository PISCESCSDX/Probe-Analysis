function [PolyIsat, voltlim, ivoltlim, Vf] = subPlasmaIsatPoly(data)
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


x = data.voltage; dx = (x(end)-x(1)) / numel(x); indx = 1:numel(x);
y = data.current;

% Define Left hand side Limit:
indcutleft = 1+round(0.02*length(y));

% Find Right Side of the Ion Saturation Region:
%   Smooth strongly -> Take absolute value -> Minimum is Floating Pot.
%   -> Maximum of part up to Vf is end of Ion saturation region
% Smooth over 29 Volts (empiric value)
  yfilt30 = abs( filtmooth(y, floor(30/dx) ) );
 
fi = 1; status = 0;
while status==0
 hef = filtmooth(y, floor(fi/dx) );
 [x0, dir] = zeropoints(x, hef);  
 % Only take the zero crossings that go from + to - (direction = -1)
 ind = dir==-1;

 if ~isempty(ind)
   if sum(ind) == 1 
     % Floating potential:
     x0 = x0(ind);
     Vf.V = x0(1);

     a = x<Vf.V; a = indx(a);
     Vf.i = a(end);

     status = 1;

   else % 
     fi = fi + 1;
     if fi > 30       % Filter range to large. No Output!
       PolyIsat = [];
       voltlim = [];
       ivoltlim = [];
       Vf = [];
       return
     end
   end
 else
   PolyIsat = [];
   voltlim = [];
   ivoltlim = [];
   Vf = [];
   return
 end
end
   
% Decide which range is the ion saturation current
% Idea: Compare one strongly with one weakly filtered characteristic
    [~, imin] = min(yfilt30);
  ind = 1:imin;
  yfilt30 = yfilt30(ind);
  % Filter characteristic (over 10 Volts)
  yfilt10 = filtspline(x, y, floor(10/dx), 0.5);
  % Compare yfilt30 and yfilt10
  [~, imindiff] = min( yfilt30-yfilt10(ind) ); % added 2013-01-30
  ind = ind(1):imindiff;
    yfilt30 = yfilt30(ind);
  ind2 = yfilt30 > yfilt10(ind);
indcutright = max(ind(ind2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot for LaTeX Documentation about the probe m-files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figeps(12,8,1,0,100); clf
% axes('position', [0.15 0.15 0.80 0.80])
% hold on
% plot(x,      y*1e3, 'o-', 'Color', 0.8*[1 1 1])
% plot(x,abs(yfilt10*1e3), 'b-')
% plot(x,abs( filtmooth(y, floor(30/dx) ) )*1e3, 'r-')
% hold off
% set(gca, 'xlim', [x(1) x(numel(x))], 'ytick', -2:2)
% ylim = get(gca, 'ylim');
% line(x(indcutright)*[1 1], ylim, 'Color', 'k')
% mkplotnice('voltage (V)', 'I (mA)', 12, '-25', '-35');
% dy = -20;
% puttextonplot(gca, [0 1], 5, -15+0*dy, 'raw data', 0, 10, 0.8*[1 1 1]);
% puttextonplot(gca, [0 1], 5, -15+1*dy, 'y_{filt10}', 0, 10, 'b');
% puttextonplot(gca, [0 1], 5, -15+2*dy, 'y_{filt30}', 0, 10, 'r');
% print('-depsc2', 'subPlasmaIsatPoly.eps')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If right side limit is not found correctly:
if indcutright <= indcutleft
  disp('subPlasmaIsatPoly_v2: No Output: I_i,sat limits not detected!')
  PolyIsat = [];
  voltlim  = [];
  ivoltlim = [];
  return
end


% Calculate linear regression of ion saturation part
% p = polyfit(x(indcutleft:indcutright), yfilt10(indcutleft:indcutright), 1);
p = polyfit(x(indcutleft:indcutright), y(indcutleft:indcutright), 1);
% Calculate linear regression of the whole characteristic
q = polyfit(x, y, 1);

% Exclude increasing ion saturation currents beyond a certain limit
% p at p(1)/q(1) < -0.05 seems to be too large(!)
if p(1)/q(1) > -0.05
  PolyIsat = p;
  voltlim = [x(indcutleft) x(indcutright)];
  ivoltlim = [indcutleft indcutright];
else
  disp('subPlasmaIsatPoly_v2: No Output: I_i,sat too strongly increasing!')
  PolyIsat = [];
  voltlim  = [];
  ivoltlim = [];
  Vf       = [];
  return
end



% Exclusion Tests: I-V characteristic needs to be asymmetric enough
%--------------------------------------------------------------------------
% Langmuir characteristic: Some Simple Symmetry Tests
% Average of left side must be larger than average of right side
  ye = y - polyval(p,x);

i_le = 1:round(length(x)/2);
i_ri = round(length(x)/2)+1:length(x);

ye_avg_le = mean(ye(i_le));
ye_avg_ri = mean(ye(i_ri));

ye_std_le = std(ye(i_le));
ye_std_ri = std(ye(i_ri));


if ye_avg_le <= ye_avg_ri
  disp('subPlasmaIsatPoly_v2: No Output: Found I_e lower on left half.')
  PolyIsat = [];
  voltlim  = [];
  ivoltlim = [];
  Vf       = [];
elseif ye_avg_ri >= 0
  disp('subPlasmaIsatPoly_v2: No Output: Positive I_e on right half.')
  PolyIsat = [];
  voltlim  = [];
  ivoltlim = [];
  Vf       = [];
elseif ye_std_ri/ye_std_le < 0.5             % found empiric
  disp('subPlasmaIsatPoly_v2: No Output: Std(I_e) too low on right half.')
  PolyIsat = [];
  voltlim  = [];
  ivoltlim = [];
  Vf       = [];
elseif abs(ye_avg_ri)/abs(ye_avg_le) < 5     % found empiric
  disp('subPlasmaIsatPoly_v2: No Output: Avg(I_e) too low on right half.')
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


%XXXXXXXXXXXXXXXXXXXXXXXXXX --- CHECKS ---  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% disp(['avg le ri: ' sprintf('%.1e',ye_avg_le) ...
%   ' ' sprintf('%.1e',ye_avg_ri)])
% disp(['avg ri/le: ' sprintf('%.1e',ye_avg_ri/ye_avg_le) ])
% 
% disp(['std le ri: ' sprintf('%.1e',ye_std_le) ...
%   ' ' sprintf('%.1e',ye_std_ri)])
% 
% disp(['std ri/le: ' sprintf('%.1e',ye_std_ri/ye_std_le) ])
% 
% %--->>> Check Plot
% close all; figure
% subplot(2,1,1)
% hold on
% plot(y,'ko')
% plot(indcutleft:indcutright,y(indcutleft:indcutright),'bo')
% plot(polyval(p,x),'r')
% hold off
% 
% subplot(2,1,2)
% plot(ye)
% %---<<<
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

end