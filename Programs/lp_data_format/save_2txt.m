%This pogram saves data from the lp_radialscan script in a two .txt file
%(comma delimited) format. The bias voltage data will be in a file named
%Vin(section #).txt and the output voltage data will be in a file named
%Vout(section #).txt
%
%Output files designed for Igor (CSDX program) probe anaylsis on PISCES lab
%computer. (Shota's code) and MATLAB analyze_lp function (Mike's code).

mkdir(datafilename);

[a b] = size(Vwrite);

for i = 1:b

filestrIN = sprintf('V%.fin.txt', i);
filestrOUT = sprintf('V%.fout.txt', i);

dlmwrite(filestrIN, Vwrite{i});
dlmwrite(filestrOUT, Iwrite{i});
movefile(filestrIN,datafilename);
movefile(filestrOUT,datafilename);

end

