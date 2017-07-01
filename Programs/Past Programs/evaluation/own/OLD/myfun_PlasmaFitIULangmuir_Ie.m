function y = myfun_PlasmaFitIULangmuir_Ie(c,x)
%==========================================================================
% This is a function to be used in the fit lsqcurvefit
% I_e(V) = -e0*n_e*A_p*sqrt(e0*Te_V/(2*pi*m_e))*exp(-(Phi-V)/Te_V)
%
% Rearrange for fit function f(x):
%
%  I_e(V) = -e0*A_p*sqrt(e0/(2*pi*m_e))*n_e*sqrt(Te_V)*exp(-(Phi-V)/Te_V)
%
% Fit function with 3 fit parameters: a:n_e, b:Te_V, c:Phi
% f(x) = c(1)*sqrt(c(2))*exp(-(c(3)-x)/c(2))
% c(1): -e0*A_p*sqrt(e0/(2*pi*m_e))*n_e
% c(2): electron temperature
% c(3): plasma potential
%==========================================================================

y = c(1)*sqrt(c(2))*exp(-(c(3)-x)/c(2));

end