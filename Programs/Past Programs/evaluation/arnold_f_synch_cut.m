function [ile fle iri fri] = arnold_f_synch_cut(f_synch, fdw)
%20080225-Mo-18:12 Brandt
%function [fsy_start fsy_end] = arnold_f_synch_cut(f_synch, fdw)
% Zusammenhängendes Synchronisationsgebiet rauscutten.

    [i_lo] = find(f_synch <  fdw );
    [i_up] = find(f_synch >= fdw );
    i_lo = i_lo(end);
    i_up = i_up(1);
    % limit to the left
    if i_lo > 1
        status=0;
        while status == 0
          i_lo = i_lo -1;
          if i_lo==1 | isnan(f_synch(i_lo))
            status = 1;
          end
        end
    end
    ile = i_lo+1; fle = f_synch(ile); 
    %
    %
    % limit to the right
    if i_up < length(f_synch)
        status=0;
        while status == 0
          i_up = i_up +1;
          if i_up==length(f_synch) | isnan(f_synch(i_up))
            status = 1;
          end
        end
    end
    iri = i_up-1; fri = f_synch(iri);
    
end