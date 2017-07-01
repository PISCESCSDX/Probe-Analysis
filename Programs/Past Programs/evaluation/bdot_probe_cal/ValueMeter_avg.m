function [amp_out, pha_out, fkHz_out, amp_std_out, pha_std_out, tstspec]= ...
            ValueMeter_avg(time, sig, zeropad, fkHz, tlen)
%function [amp_out, pha_out, fkHz_out]= ...
%            ValueMeter_avg(time, sig, zeropad, fkHz, tlen)
%
%   output: amp_out         amplitude (NOT A_eff!) in units of sig
%           pha_out         Phase in deg
%           fkHz_out        exact frequency found in signal 1
%           amp_std_out     standard deviations of amplitudes
%           pha_std_out     standard deviations of phases
%           tstspec         array of averaged spectra (for visualization) 
%                           containing fields .fvec [kHz] .amp .phi
%
% cross phases are calculated relative to sig(:, 1); exactly the same
% frequency as found in sig1 is used
%       sig1 sig2 sig 3 ...
%   t
%   i   sig
%   m
%   e
%

   if nargin == 0
       help ValueMeter_avg;
       return
   end

   time= time - time(1);
   avs= floor(time(end) / tlen);
   sigs= size(sig, 2);
   
   ind_num= sum(time < tlen) - 1;
   for i1= 1:avs
       i_st= (i1 - 1) .* ind_num + 1;
       i_en= i1 .* ind_num;
       idx(i1, :)= [i_st:i_en];
   end

   %-----measure first signal----------------------------------
   for i1= 1:avs
       [amp(i1, 1), pha(i1, 1), fre(i1), spe(i1, 1)]= ...
                 ValueMeter(time(idx(i1, :)), ...
                            sig(idx(i1, :), 1), ...
                            zeropad, ...
                            fkHz .* 1000, ...
                            false);
   end
   
   %-----measure remaining signals relatively to first one-----
   fkHz_out= mean(fre) ./ 1000;
   for i2= 2:sigs
       for i1= 1:avs
           [amp(i1, i2), pha(i1, i2), ftrash, spe(i1, i2)]= ...
                    ValueMeter(time(idx(i1, :)), ...
                               sig(idx(i1, :), i2), ...
                               zeropad, ...
                               fkHz_out .* 1000, ...
                               true);
       end
   end
   
   %------before averaging phases add/sub 2pi as needed--------
   pha= pha - pha(:, 1) * ones(1, size(sig, 2));
   for i2= 1:sigs
       [y_hist, x_phi]= hist(pha(:, i2), 20);
       [mval, max_ind]= max(y_hist);
       phi_est= x_phi(max_ind);
       add= (pha(:, i2) - phi_est) < -1.*pi;
       sub= (pha(:, i2) - phi_est) > pi;
       pha(:, i2)= pha(:, i2) + 2.*pi.*add - 2.*pi.*sub;
   end
   
   %------test spectrum: averaged amplitudes and phases from 1st time
   %------series
   for i2= 1:sigs
       tstspec(i2).amp= zeros(size(spe(1,i2).ampvec));
       for i1= 1:avs
           tstspec(i2).amp= tstspec(i2).amp + spe(i1, i2).ampvec; 
       end
       tstspec(i2).amp=  tstspec(i2).amp ./ avs;
   end
   tstspec(1).fvec= spe(1).fvec' ./ 1000;   %Hz -> kHz
   for i2= 2:sigs
       tstspec(i2).amp=  tstspec(1).amp .* tstspec(i2).amp;
       tstspec(i2).phi=  spe(1, 1).phivec - spe(1, i2).phivec;
       tstspec(i2).phi=  mod(unwrap(tstspec(i2).phi), 2.*pi) - pi;
       tstspec(i2).fvec= spe(i2).fvec' ./ 1000;   %Hz -> kHz
   end
   
   
   amp_out= mean(amp, 1);
   pha_out= mean(pha, 1);
   amp_std_out= std(amp, 1);
   pha_std_out= std(pha, 1);
end

