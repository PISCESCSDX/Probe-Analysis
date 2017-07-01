function [meas_mat, timevec, t_stamp]= mdf_load_64ch(basedir, idx_nr, ...
                                                     file_list)
% function [meas_mat, timevec, t_stamp]= mdf_load_64ch(basedir, idx_nr)
% 
% load a measurement from transient recorder up to 64 channel
% choose measurement by basedir and the number of brd 1; find
% connected files within one measurement by preprocess/filebundle.dat
% set all ch in preprocess/mask_ch to NaN
%
% input       basedir         string containing dir to act on
%             idx_nr          running number of brd 1
% output      meas_mat        measured voltages [smpl x channel]
%             timevec         time in s [smpl x 1]
%             t_stamp         measurement time as matlab number
%

    if nargin < 3
        [name_mat, file_time, idx_vec]= ...
            mdf_load_filebundle_list(basedir);
    else
        name_mat = file_list.name_mat;
        file_time= file_list.file_time;
        idx_vec  = file_list.idx_vec;
    end

    if basedir(end) ~= '/', basedir= [basedir, '/']; end
    
    meas_cellmat= {};
    ind= idx_vec == idx_nr;
    smpl_num= [];
    
    for i= 1:size(name_mat, 2);
        name= name_mat{ind, i};
        if  ~strcmp(name, '[]');
            [meas_cellmat{i}, timevec]= readmdf([basedir, name]);
            smpl_num= size(meas_cellmat{i}, 1);
        else
            meas_cellmat{i}= [];
        end
    end
    
    meas_mat= [];
    for i= 1:size(name_mat, 2)
        if size(meas_cellmat{i}) ~= [0,0]
            meas_mat= [meas_mat, meas_cellmat{i}];
        else
            meas_mat= [meas_mat, NaN(smpl_num, 8)];
        end
    end

    %mask channel
    if size(dir([basedir, 'preprocess/mask_ch.dat']), 1) == 1
        [ln, ln_string]= find_ascii_token(...
                         [basedir, 'preprocess/mask_ch.dat'], ...
                         'mask ch: ' );
        mask_ch= sscanf(ln_string, '%d');
        meas_mat(:, mask_ch)= NaN;
    end
    
    %---------time stamp
    t_stamp= file_time(ind);
    
end