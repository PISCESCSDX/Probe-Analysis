function [t0, t1, i0, i1] = pulse_limits(tt, sig, std_threshold)
%function [t0 t1 i0 i1] = pulse_limits(tt, sig, std_threshold)
% Find the begin and end time of a pulse in a timetrace.
% ONLY 1 PULSE PER TIME TRACE! PULSES IN BETWEEN WON'T BE DETECTED.
% IN: tt: time vector
%     sig: time row
%     std_threshold: pulse registered if "sig" bigger than 
%                    std_threahold percentage of std
%OUT: t0, t1: time pulse start
%     i0, i1: indices
% EX:

if nargin < 2; error('More input!'); end
if nargin < 3; std_threshold = 0.1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PULSE DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  i_sig = find(abs(sig) > std_threshold);

  
if isempty(i_sig)
  t0 = tt(1);
  t1 = tt(end);
  i0 = 1;
  i1 = length(tt);
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PULSE START & END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  t0 = tt(i_sig(1));
  t1 = tt(i_sig(end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INDICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  i0 = i_sig(1);
  i1 = i_sig(end);
end
 

end