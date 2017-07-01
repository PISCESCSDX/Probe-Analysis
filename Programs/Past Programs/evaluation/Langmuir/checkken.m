function checkken
% Checks all probe chrakteristics in a directory.
% manually in konsole with:
%chfit -fit=kinetic_mag -l=0.004 -r=0.0001 -b=0.0316 -n=2 -t=-200:100 meas_000000.dat 
%chfit -fit=kinetic_mag -l=0.004 -r=0.0001 -b=0.0316 -t=-190:190 -n=5 -f=1 -ho=2 meas_000002.dat

fls = dir;
if isstruct(fls)
    for i = 1:length(fls)
        tmp(i,:) = cellstr(fls(i).name);
    end;
    fls = tmp;
end;

for i = 1:length(fls)
    if ~isempty(strfind(fls{i},'.mod'));
        ans = 'y';
        while lower(ans) == 'y'
            loadken(fls{i}(1:end-4),1);
            ans = input('Plot again Y/N','s');
        end;
        
    end;
end;
