%This pogram saves data from the lp_radialscan script in a one .txt file
%(comma delimited) format. The data will be saved in two columns in a file
%name V(section #).txt. The first column will be the Vin (bias voltage)
%data and the second column will be the Vout (output voltage or current)
%data. This file will have headings.
%
%Output files designed for Saiakt's IDL probe analysis program



[a b] = size(Vwrite);

for i = 1:1

    D(:,1) = Vwrite{i};
    D(:,2) = Iwrite{i};

    titlestr = sprintf('V%.f.txt', i);

    fid = fopen(titlestr, 'w');
    fprintf(fid, 'Vin   Vout\n');
    fclose(fid);
    
    dlmwrite(titlestr, D, '-append')
    %comma delimiter is default, '\t' for tab delimited

end

