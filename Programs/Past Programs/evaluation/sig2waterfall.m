function [tvec, fkHz, s]= sig2waterfall(tsig, sig, df, fkHz_en)
%function [tvec, fkHz, s]= sig2waterfall(tsig, sig, df, fkHz_en)
%
% calculates subsequent fft after windowing original transient
% each fft is calculated with an own window (not rectangular)
%
% input:  tsig, sig          time [s] and signal, standing vectors
%         df                 frequency resolution [Hz] determines window
%                            size
%         fkHz_en            end of frequency range in kHz (is cutted out)
% output: tvec, fkHz         vectors for x- and y-axes
%         s


    tsig= tsig - tsig(1);
    dt_all= tsig(end)-tsig(1);
    w_num= dt_all/ (1/df);
    dt_wdw= dt_all ./ w_num;
    wdw_len= sum(tsig < dt_wdw) - 1;
    tsig_wdw= tsig(1:wdw_len);
    
    %--------------------------------------------calculate spectra
    for i1= 1:w_num
        i_st= (i1 - 1).* wdw_len + 1;
        i_en=  i1     .* wdw_len;
        sig_wdw= sig(i_st:i_en);
        [fvec, spec]= fft_spectrum(tsig_wdw, sig_wdw, 0, 2);
        s(:, i1)= spec;
        tvec(i1)= tsig( round((i_st + i_en)/2) );
    end
    
    %--------------------------------------------cut spectrum
    f_ind= (fvec < 1e3.* fkHz_en);
    fkHz= 1e-3.* fvec(f_ind);
    s= s(f_ind, :);

end