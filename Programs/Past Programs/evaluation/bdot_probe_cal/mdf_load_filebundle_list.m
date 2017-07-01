function [name_mat, timevec, idx_vec]= mdf_load_filebundle_list(basedir)
%function [name_mat, timevec, idx_vec]= mdf_load_filebundle_list(basedir)
%
%loads matrix of *.MDF files from preprocess/filebundle.dat, which
%was created by mdf_filebundle_list.m
%
%input          basedir         string containing dir to act on
%output         name_mat        cell matrix of filenames,
%                                rows = indep. measurements
%                                cols = files within one measurement
%               timevec         vector of times of measurement (brd 1)
%                               in matlab time format (days from 1.1.0000)
%               idx_vec         vector of numbers of board 1

    if basedir(end) ~= '/', basedir= [basedir, '/']; end
    
    % checking for / creating filebundle.dat
    f= dir([basedir, 'preprocess/filebundle.dat']);
    if size(f, 1)== 0
        fprintf(1, ['creating preprocess/filebundle.dat', ...
                   ' with mdf_filebundle_list.m\n']);
        mdf_filebundle_list(basedir);
    end


    fid= fopen([basedir, 'preprocess/filebundle.dat'], 'rt');
    tline= fgetl(fid);
    while tline(1)== '%'
        tline= fgetl(fid);
    end
    
    timevec_str= [];
    name_mat= {};
    idx_vec= [];
    while 1
        c_line= line2cellstr(tline);
        name_mat= [name_mat; c_line(5:end)];
        timevec_str= [timevec_str; c_line{2}, ' ', ...
                                   c_line{3}, ' ', ...
                                   c_line{4}];
        idx_vec= [idx_vec; sscanf(c_line{1}, '%d')];
        tline= fgetl(fid);
        if ~ischar(tline),   break,   end
    end
    
    %-----------------------------------check for empty columns
    for i1= 1:size(name_mat, 1)
        for i2= 1:size(name_mat, 2)
            name_exist(i1, i2)= strcmp(name_mat{i1, i2}, '[]');
        end
    end
    col_ind= (0 ==sum(name_exist, 1));
    name_mat= name_mat(:, col_ind);
    %----------------------------------and kill them-----------
    
    timevec= datenum(timevec_str, 'dd.mmm yyyy HH:MM:SS');
    
    fclose(fid);
    
end