function autokenfit(fstr, pl, pr, par, B, tr, nvec, lovec, hovec, ...
  fvec, epsabs, epsrel, meth)
%==========================================================================
%function autokenfit(fstr, pl, pr, par, B, tr, nvec, lovec, hovec, ...
%  fvec, epsabs, epsrel, meth)
%--------------------------------------------------------------------------
% AUTOKENFIT fits all probe rawdata files with chfit.
%--------------------------------------------------------------------------
% IN: fstr: filename string of all probe data files, e.g. 'meas*.dat'
%     pl: probe length (m)
%     pr: probe radius (m)
%    par: [0, 1] probe orientation: perpendicular (0), parallel (1)
%      B: B-field (T)
%     tr: [string] truncate, ex: '-40:15' (V)
%     meth: 1: take best fit for IV-characteristic
%           2: take best fit for IV-characteristic 
%           3: take best fit for sum of IV- and dI/dV
%OUT: evaluation mat-files for all fits
%--------------------------------------------------------------------------
% EX: autokenfit(fstr,pl,pr,par,B,tr,nvec,lovec,hovec,...
%        fvec,epsabs,epsrel,meth)
%==========================================================================

if nargin<10
  meth = 1;
end

a = dir(fstr);
la = length(a);

fitn{la}  = [];
fitf{la}  = [];
fitlo{la} = [];
fitho{la} = [];
fitmat{la}= [];

% FOR all files
for i=1:la
  ctr=1;

rmat = zeros(length(fvec)*length(nvec)*length(hovec), 5);
% FOR all filterstrengths
for  f = fvec
% FOR number of moving avg points
for  n = nvec
% FOR voltage offset
for ho = hovec
% FOR voltage offset
for lo = lovec
  
  % Calculate fit: Do not save any files
  [iufit, fitpar] = ...
  chfit_kinetic_mag(a(i).name,pl,pr,par,B,tr,ho,n,f,epsabs,epsrel,0);

  rmat(ctr,1) = fitpar.nn;
  rmat(ctr,2) = fitpar.ff;
  rmat(ctr,3) = fitpar.lo;
  rmat(ctr,4) = fitpar.ho;
  if isempty(iufit)
    rmat(ctr,5) = -Inf;
    rmat(ctr,6) = -Inf;
  else
    rmat(ctr,5) = iufit.R;
    rmat(ctr,6) = iufit.RDer;
  end


  ctr=ctr+1;
end % FOR lo
end % FOR ho
end % FOR n
end % FOR f


% CHECK SUM ERROR (of IV-characteristic and dI/dV-curve)
%-------------------------------------------------------
[~, iR   ] = max(rmat(:,4));
[~, iRDer] = max(rmat(:,5));


% CHECK SUM ERROR (of IV-characteristic and dI/dV-curve)
%-------------------------------------------------------
SumRRDer = rmat(:,5) + rmat(:,6);
[~, iSum] = max(SumRRDer);

fitn{i}  = [rmat(iR,1) rmat(iRDer,1) rmat(iSum,1)];
fitf{i}  = [rmat(iR,2) rmat(iRDer,2) rmat(iSum,2)];
fitlo{i} = [rmat(iR,3) rmat(iRDer,3) rmat(iSum,3)];
fitho{i} = [rmat(iR,4) rmat(iRDer,4) rmat(iSum,4)];

% EVALUATE THE BEST FIT AGAIN
%----------------------------
chfit_kinetic_mag( a(i).name, pl, pr, par, B, tr, fitlo{i}(meth), ...
  fitho{i}(meth), fitn{i}(meth), fitf{i}(meth), epsabs, epsrel, 1);

fitmat{i} = rmat;

end

info = 'fitmat(:,i)=n, f, lo, ho, Rsqr, RDer'; %#ok<NASGU>
save('allmat.mat', 'fitmat', 'info');

end