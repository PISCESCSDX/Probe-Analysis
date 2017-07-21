%% lp_radialscan summary 2017_07_17 - Brandon Schwendeman
%
%This program reads langmuir probe data in a .h5 file that was recorded
%while conducting a radial scan in CSDX. The probe was moved through the
%plasma while continously preforming voltage bias sweeps (ideally from high
%to low) and the data needs to already cut into sections for analysis. Each
%section contains multiple bias sweeps to be averaged over that will
%represent a radial location in the plasma. This will smooth out the final IV
%curves when analyzing. The data in the .h5 file should be in the form of a matrix, 
%where each column represents a section of the data (each column should be 
%a continous set of bias sweeps). The bias voltage should be in the "V"
%section of the .h5 file and the output voltage (across a resister, or
%current in amps) should be in the "I" section. Finally, the experimental
%voltage offset should already have been used to adjust the data in the .h5
%file.
%
%****This program can use either output voltage across a resister or the
%current. The program and function will not attempt to alter the data. Make
%sure the error file (noise file) and the data are in the same units and if
%using an output voltage, the data will need to be adjusted during the
%analysis.*****
%
%This program will analyze each individual section, locate each bias sweep,
%cut out the transitions between bias sweeps, and then average all the
%sweeps in a section to obtain Vin and Vout data for each section (radial
%location). It will follow the same procedure for both the data and the
%error/noise file, and then subtract the noise from the data.
%
%This program will prompt the user for:
%   Data file name (.h5 format)- .h5 file name and file path (if necessary)
%       that contains the probe data
%   Error file name (.h5 format) - .h5 file name and file path (if
%       necessary) that contains the electronic noise data (taken with no
%       plasma and magnetic field active)
%   Bias sweeps per radial section - the approximate number of bias sweeps
%       contained in each section, error on the high side to include all the
%       data
%   Index interval size for min/max search - this is the size of the index
%       interval used to search for the minimum and maximums of a bias sweep.
%       A higher value will make the program run faster, but might cause it to
%       miss data at the beginning and end of the section or miss transition
%       sections. Adjust based on the number of data points taken.
%   Voltage step for averaging data - This is the bias voltage step that
%       the output voltage (or current) will be averaged over. A larger value
%       will make the program run faster and smooth out more noise/flucuations
%       in the data. However, too large of a value will lower the resolution of
%       the IV curve and make it more difficult to identify the necessary
%       trends for probe analysis. Adjust as necessary.
%
%The final data will be in the CELL ARRAY variables:
%   Vwrite (input bias voltage data) &
%   Iwrite (output voltage or current, whichever units it was in
%   originally)
%Each column of the data will represent a radial section that multiple bias
%sweeps were averaged to obtain (the number of columns will equal the number
%of sections the data was cut into). The index values of Vwrite and Iwrite
%correspond to one another.


%% File info and program parameters input

%defualt responses to user prompts
datafilename = '(path, if necesarry)/filename.h5'; %'(path, if necesarry)/filename.h5'
errorfilename = 'B0_24914_port1.h5';  %'(path, if necesarry)/filename.h5'
%a file path is necessary if the file is not in the current MATLAB folder
sweepsperintervalSTR = '10';
minmaxindexstepSTR = '2000';
averagestepSTR = '0.1';   %Volts

defaultAns = {datafilename, errorfilename, sweepsperintervalSTR, minmaxindexstepSTR, averagestepSTR};

%set prompt responses to default (these can be changed by the user when
%prompted)
promptAns = defaultAns;

%set up prompt dialog box
    %prompt dialog
    prompt = {'Data file name (.h5 format):',...
          'Error file name (.h5 format):',...
          'Bias sweeps per radial section:',...
          'Index interval size for min/max search:',...
          'Voltage step for averaging data:'};
    %prompt title  
    name = 'File info/program parameters';
    %set number of lines for response boxes to one
    numlines = 1;
    %prompt user for file information and program parameters (defualt
    %answers will be displayed and can be changed
    promptAns = inputdlg(prompt, name, numlines, promptAns);

%convert data from input prompts to usable state for program
datafilename = char(promptAns(1));
errorfilename = char(promptAns(2));
sweepsperinterval = str2double(promptAns(3));
minmaxindexstep = str2double(promptAns(4));
averagestep = str2double(promptAns(5));

%% Import data and run the lp_radialscan_format function

%import data from .h5 file into two matrices. Vmatrix is Vin voltage
%biasing data and Imatrix is Vout data (voltage read across the resister)
%that is used to calculate the current flowing through the probe

Vmatrix = h5read(datafilename,'/V');
Imatrix = h5read(datafilename,'/I');

ErrV = h5read(errorfilename, '/V');
ErrI = h5read(errorfilename, '/I');

%run format function to cut out transition intervals between bias sweeps and
%average the data
[Vin, Vout] = lp_radialscan_format(Vmatrix, Imatrix, sweepsperinterval, minmaxindexstep, averagestep);
[Errin, Errout] = lp_radialscan_format(ErrV, ErrI, sweepsperinterval, minmaxindexstep, averagestep);

%% Calculation of electronic error and data adjustment

%calculate number of sections in error matrices for the loop
[a b] = size(Errout);

for i = 1:b

%create loop variable to manipulate, for loops adjusts each column (so each
%cut section) one at a time
Errinloop = Errin(:,i);   
Erroutloop = Errout(:,i);

%locate the last nonzero value in the data and delete the trailing zeros
i1 = find(Erroutloop, 1, 'last');
ErroutloopFIXED = double(Erroutloop(1:i1));

ErrinloopFIXED = double(Errinloop(1:i1));

%create a linear fit for the error data
%Errfit = fit(ErrinloopFIXED, ErroutloopFIXED, 'poly1');
[m,b]=LinearRegression(ErrinloopFIXED, ErroutloopFIXED);

%remove trailing zeros from Vout data and subtract the electronic noise
%data to create adjusted data matrix
i2 = find(Vout(:,i),1,'last');

%VoutAdjusted{i} = Vout(1:i2,i) - Errfit(Vin(1:i2,i));
VoutAdjusted{i}=zeros(i2,1);
for j=1:i2
VoutAdjusted{i}(j,1) = Vout(j,i) - (m*(Vin(j,i))+b);
end
end

%% Final data output for saving

%VoutAdjusted matrix already created in prevoious section
Iwrite = VoutAdjusted;

%delete trailing zeros from Vin data
[a b] = size(Vin);
for i = 1:b

A = Vin(:,i);
i1 = find(A, 1, 'last');
Vwrite{i} = A(1:i1);

end

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

%final variables are Vwrite (adjusted Vin data) and Iwrite (adjusted Vout
%data)....remember Iwrite can represent either the output voltage across
%the resister or the current if that calculation has already been made.
%This program will not make that calculation

