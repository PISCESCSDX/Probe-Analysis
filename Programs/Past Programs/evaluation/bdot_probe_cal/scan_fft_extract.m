function scan_fft_extract(base)
%function scan_fft_extract(base)
%
% Input:    basedir     top level directory of measurement; one subdir
%                       must be 'TR_Frisch' _OR_ number of directory
%                       given to scandir
%
    if nargin == 0
        help scan_fft_extract;
        return
    end

    %------------if a number is given as input change to directory
    struct= whos('base');
    if strcmp(struct.class, 'double')
        base= scandir(base);
    end

    if base(end)~='/', base= [base, '/'];end
    fstruct= dir([base, '*.def']);
    if size(fstruct, 1) ~= 1
        error('def-file not found or more than one.');
    end
    deffile= [base, fstruct.name];
    fstruct1= dir([base, 'xypos*.dat']);
    fstruct2= dir([base, 'rzpos*.dat']);
    if (size(fstruct1, 1) == 0 )&...
       (size(fstruct2, 1) == 0 )
        error('no rzpos*.dat and no xypos*.dat found');
    end
    
    hdl= conf_parser(deffile);

    fkHz   = hdl.fkHz;
    t_st   = hdl.t_st;
    t_en   = hdl.t_en;
    header = hdl.header;
    f_str  = hdl.f_str;
    ch_list= hdl.ch_list;
    anzch  = length(ch_list);

    frisch_dir= [base, 'TR_Frisch/'];
    pre_dir   = [base, 'process/'];
    diag_dir  = [base, 'process_diag/'];
    outname=  [pre_dir, 'c_amp_Bprobes1.dat'];
    outname2= [pre_dir, 'c_amp_Bprobes2.dat'];
    flist= dir([frisch_dir,'*.MDF']);

    make_folder(pre_dir);
    make_folder(diag_dir);

    %-----------load mdf-pos-list
    [file_list.name_mat, file_list.file_time, file_list.idx_vec]= ...
            mdf_load_filebundle_list(frisch_dir);
    anz= length(file_list.idx_vec);
    
    %-----------remove outfile if exists
    fs= dir(outname);
    if size(fs, 1)~= 0 
        delete(outname);
    end
        
    fprintf(1, '                ');
    for i1= 1:anz
        %--------------------------------------------display info
        fprintf(1, ['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b', ...
                    'file ', num2str(i1, '%05d'), ...
                        '/', num2str(anz, '%05d')]);
        
        %--------------------------------------read and calculate
        mdfnum= str2num(flist(i1).name(3:8));
        [mat, tvec, t_stamp]= mdf_load_64ch(frisch_dir, mdfnum, file_list);
        %[mat, tvec]= readmdf([frisch_dir, flist(i1).name]);    
        ind_st= sum(tvec <= t_st);
        ind_en= sum(tvec <= t_en);
        
        mat2= mat(ind_st:ind_en, ch_list);
           
        tvec2= tvec(ind_st:ind_en);
        [A1, phi1, f, s_A1, s_phi1, tstspec]= ...
            ValueMeter_avg(tvec2, mat2, 5, fkHz, 0.05);
        %time_spec_indicator(tvec, mat(:, ch_list(1)), ...
        %                    tvec2, mat2(:, 1), ...
        %                    tstspec, f);
        mat3= (mat2 - (ones(size(mat2, 1), 1) * mean(mat2)) )...
                   ./ (ones(size(mat2, 1), 1) * std(mat2));
        % maybe this has to be changed:
        % sorted not found on Akademgorodok + direct index "2:4"
        time_spec_indicator(tvec2(1:300), mat3(1:300, :), ...
                            tvec2(1:300), mat3(1:300, 1), ...
                            tstspec, f);
        outname_diag= [diag_dir, 'diag_mdf', num2str(mdfnum, '%06d'), ...
                       '.eps'];
        print('-depsc2', outname_diag);
        %---------------------------------------------save result
        outline= [i1, t_stamp, f, ...
                    A1(1), ...
                    reshape([A1(2:end); ...
                             phi1(2:end)], 1, 2.*anzch -2) ...
                    ];
        
        save_linefile(outname, header, outline, f_str);
        
    end
    fprintf(1, '\n');
    
    merge_process_poslist(outname, outname2);

end