% DEMONSTRATION HOW THE PHASE of a signal can be determined

% CREATE artificial time signal
tt = (0:30e3)/10e6;
tl = length(tt);

% FREQUENCY
f = 3e3;

sig = sin(2*pi*f.*tt) + 0.1*randn(tl,1)';

% \citep{SunJ2008}: description of Hilbert transform
% The phase is determined by calculation of the analytical signal sa(t) of 
% a time series s(t). sa(t)=s(t)+i*h(t)
% h(t) is the Hilbert transform.
% The Matlab Hilbert function yields directly sa(t)!
% That means matlab: hi = hilbert(s(t)) = sa(t), is a complex function.
% The phase is simply determined by angle(hi)!

hi = hilbert(sig);
phi  = angle(hi);
phi0 = unwrap(phi)/pi+0.5;

hold on
plot(tt, sig, 'r')
plot(tt, phi0)
hold off