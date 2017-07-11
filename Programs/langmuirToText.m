function langmuirToText(filebase,zeroFileBase)
%Part Two of the formatAndAverageLangmuirData program. Turns Langmuir
%data into two text documents that can be fed into the "analyze_lp" program

%Eliminates noise. File name is subject to change.

[V1,I1]= formatAndAverageLangmuirData(filebase);
[V0,I0]= formatAndAverageLangmuirData(zeroFileBase);
V=V1-V0;
I=I1-I0;

for i=1:length(V(1,:))
    dlmwrite([filebase '_Vin_' int2str(i) '.txt'],V(:,i),'delimiter','\t');
    dlmwrite([filebase '_Iout_' int2str(i) '.txt'],I(:,i),'delimiter','\t');
end

end

