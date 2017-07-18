function [ Vin, Iout ] = lp_radialscan_format( Vmatrix, Imatrix, sweepsperinterval, minmaxindexstep, averagestep)
%lp_radialscan_format: % This function reads langmuir probe data collected
%while preforming a radial scan and reformats it for analysis. Probe data
%is expected to be collected in continous bias sweeps (from positive to
%negative ideally) while moving from the edge through the center of the
%plasma. The data should be cut into multiple sections, with each column in
%the data matrices representing a section that will be averaged over
%multiple sweeps to smooth out the results. This function will remove the
%transition region between bias sweeps, organize and average the data over
%the number of bias sweeps per section, and output the data in a similar
%format with each column representing the averaged bias voltage and output
%voltage (or current) for each section.
%
%Format: [ Vin, Iout ] =
%lp_radialscan_format( Vmatrix, Imatrix, sweepsperinterval,
%minmaxindexstep, averagestep) with:
%
%Vmatrix - bias voltage (Vin) matrix with each section of bias sweeps in a
%   coloumn
%Imatrix - corresponding voltage (or current, Vout or Iout) matrix with
%   same format and index values as Vmatrix
%sweepsperinterval - approximate number of voltage sweeps in each interval,
%   (error on the high side if necessary)
%minmaxindexstep - the index step size used to locate the minimums and
%   maximums of the bias sweeps (a larger value will make the function run
%   faster but also means it could miss data at the beginning and end of each
%   section, default = 2000)
%averagestep - the bias voltage step that the data is averaged over in each section 
%   once it is seperated into individual bias sweeps (a larger value will make
%   the program run faster and reduce more noise but will also eventually
%   result in less accurate IV data, default = 0.1 V)


%% Idenitify voltage bias sweeps and corresponding current values in each data interval

%calculate the number of sections the data has been cut into (these
%correspond to the radial positions inside the plasma)
[a b] = size(Vmatrix);
numberofintervals = b;

%first loop to seperate data intervals (coloumns of Vmatrix and Imatrix) to
%work on one at a time
%--------------------------------------------------------------------------
for i = 1:numberofintervals

     %create nx1 vectors of data to work with rather than matrix column for
     %the loop iteration
     Vloop = Vmatrix(:,i);
     Iloop = Imatrix(:,i);
     
%Begin parsing through data to locate the minimums and maximums of the
%bias voltage sweeps. A voltage sweep will be indentified using a minimum
%and maximum and then saved in a new matrix. This will eliminate the
%transition regions between voltage sweeps and allow the data to be
%averaged in order smooth out the final IV curves. The corresponding
%current values m will also be identified and save in a matrix
%--------------------------------------------------------------------------
     
     %create initial index step to search for first maximum of bias
     %sweep, the size can be adjusted by changing minmaxindexstep when
     %user is prompted at start of program
     samplestart = 1;
     sampleend = minmaxindexstep;
     
     %identify the last index of the data for use later so the program
     %doesn't try to read data that doesnt exist when it has finished
     %parsing through all the data
     indexlimit = length(Vloop);
     
    
     %create loop that will end when the specified number of bias sweeps in
     %the data have been identified or the final data index is reached
    for j = 1:sweepsperinterval
    
        %set the current index range to the start of the data to look for
        %first maximum
        indexcountLow = samplestart;
        indexcountHigh = sampleend;
        
        %if the current range is larger than the data set, inform the user
        %that they should decrease the index interval size when prompted at
        %the start of the program
        if indexcountHigh > indexlimit
            error('Not enough data points to find minimums and maximums./nUse file with more data points or decrease min/max index interval size')
        end    
        
        %create variable including voltage values in the initial index range
        %find the max value in that range
        Vbefore = Vloop(indexcountLow:indexcountHigh);
        currentmax = max(Vbefore);
        
        %move index range forward by the specificed index step
        indexcountLow = indexcountLow + minmaxindexstep;
        indexcountHigh = indexcountHigh + minmaxindexstep;
        
        %if the end of the index range exceeds the range of the data, set
        %the end of the index to the same value as the final index of the
        %data
        if indexcountHigh > indexlimit
           indexcountHigh = indexlimit;
        end
            
        %create new variable including voltage values in the new index range
        %and calculate the max value in this range
        Vafter = Vloop(indexcountLow:indexcountHigh);
        newmax = max(Vafter);

            %creates loop that will compare old max value with the new one
            %while continuing to move the index range forward until the new
            %max is smaller than the old max. At this point, a local max
            %value has been found which is the start of a voltage bias
            %sweep
            while currentmax <= newmax
        
                %move index range forward
                indexcountLow = indexcountLow + minmaxindexstep;
                indexcountHigh = indexcountHigh + minmaxindexstep;
                
                %set current maximum to new maximum
                currentmax = newmax;
                
                %if the high end of the index range exceeds the data index
                %range, stop the loop and the current max will be recorded
                %as the local max
                if indexcountHigh > indexlimit
                    break
                end
                
                %set the old data range equal to the new range
                Vbefore = Vafter;   
                
                %recalculate the new data range with the index that has
                %been moved forward and recalculate the max for comparison
                Vafter = Vloop(indexcountLow:indexcountHigh);
                newmax = max(Vafter);
                
                %process is repeated until the a local max is indentified
                %or the end of the data range is reached
                
            end
    
            %locate the index value of the current maximum (the loop was exited
            %because the new max was lower) in the original data. This
            %index value is the start of a voltage sweep
            maxindex = find(Vbefore == currentmax,1)+indexcountLow - minmaxindexstep - 1+1;
    
            %set the index range to start one value after the start the
            %bias sweep (the local max)
            indexloopcountLow = maxindex;
            indexloopcountHigh = maxindex + minmaxindexstep;
            
            %the the high end of the index range exceeds the data range, se
            %the end of the index range to the value of the end of the data
            %range, calculate the minimum in this range, find its index
            %value in the data, save this range as the final bias sweep
            %(saves both current and voltage data for this sweep) and then
            %exit the loop because this was the end of the data has been
            %reached (similar code is explained further down in this loop)
            if indexcountHigh > indexlimit
                   indexcountHigh = indexlimit;
                   Vbefore = Vloop(indexcountLow:indexcountHigh);
                   currentmin = min(Vbefore);
                   minindex = find(Vbefore == currentmin,1) + indexcountLow - minmaxindexstep - 1;
    
                   Vsweep{j} = Vloop(maxindex:minindex);
                   Isweep{j} = Iloop(maxindex:minindex);
                   
                   break
                   
            end    
           
            %create variable for data in current range and caculate the
            %minimum value
            Vbefore = Vloop(indexcountLow:indexcountHigh);
            currentmin = min(Vbefore);
            
            %move the index range forward
            indexcountLow = indexcountLow + minmaxindexstep;
            indexcountHigh = indexcountHigh + minmaxindexstep;
            
            %check to see if the end of the data has been reached
            if indexcountHigh > indexlimit
                indexcountHigh = indexlimit;
                
            end
    
            %create variable for data in new range and calculate new minimum
            Vafter = Vloop(indexcountLow:indexcountHigh);
            newmin = min(Vafter);
    
            %similar to maximum loop above: compare the minium values in
            %the before and after ranges until the min in the after range
            %is larger than the min in the before range. At this point the
            %local minimum (and end of the bias sweep) has been located and
            %the loop should be exited. The loop will also be exited if the
            %end of the data range is reached
            while currentmin >= newmin
     
                indexcountLow = indexcountLow + minmaxindexstep;
                indexcountHigh = indexcountHigh + minmaxindexstep;
                
                currentmin = newmin;
        
                if indexcountHigh > indexlimit
                    break
                end    
                
                Vbefore = Vafter;
                
                Vafter = Vloop(indexcountLow:indexcountHigh);
                newmin = min(Vafter);
        
            end    
    
            %locate the minimum index in the orginal data and record this
            %value
            minindex = find(Vbefore == currentmin,1) + indexcountLow - minmaxindexstep - 1-1;
    
            %a sucessive local max and min have now been located, the data
            %between these points constitute a voltage sweep. This data is
            %recorded for each sweep found in the Vsweep and Isweep
            %variables. this process has cut out the transition regions
            %between bias sweeps
            Vsweep{j} = Vloop(maxindex:minindex);
            Isweep{j} = Iloop(maxindex:minindex);
            
            %set the index range to start one value after the current local
            %min to start the loop all over again to look for another
            %voltage sweep
            samplestart = minindex;
            sampleend = minindex + minmaxindexstep;
            
            %if the end of the data range has been reached, exit the loop
            %because the final bias sweep has already been located
            if sampleend > indexlimit
                break
            end    
                
    end
    
    
 %% Average the voltage bias sweeps in each interval to smooth out the data for analysis
    
 
 %SORTING CODE
 
 %the voltage and current data for each data interval must be combined into
 %one nx1 vector in order to be averaged. Each set of data will be set end
 %to end in an unsorted vector and then the voltage will be sorted in ascending
 %order and the current will be reorganized so the values correspond
 %correctly
 
    %set a low and high variable for use in loop to record index values so
    %the data can be recorded in the new vectors end to end
    low = 1;
    high = 0;
    
    for k = 1:9
     %identify the length of the current voltage sweep and set index so the
     %values can be recored from the current index location
     high = length(Vsweep{k})+ high;
     
     %record voltage and current data end to end
     VaveragingvectorUNSORTED(low:high) = Vsweep{k};
     IaveragingvectorUNSORTED(low:high) = Isweep{k};
     
     %adjust loop variables and repeat
     low = high + 1;
        
    end
    
    %sort the bias voltage vector in ascending order, preserve the index
    %change and sort of the current vector accordingly
    [VaveragingvectorSORTED, sortindex] = sort(VaveragingvectorUNSORTED, 'ascend');
    IaveragingvectorSORTED = IaveragingvectorUNSORTED(sortindex);
    
%AVERAGING CODE    
    
%the data for each interval is now combined in one long nx1 vector in
%order of ascending voltage. The data can now be averaged over a
%specified voltage step (can be set by the user when prompted as the
%program is running). 
    
    %find the minimum and maximum value of the sorted voltage data
    maxbias = max(VaveragingvectorSORTED);
    minbias = min(VaveragingvectorSORTED);
    
    %set the initial averaging range based on the minimum of the sorted
    %data and the specified averaging step
    averagelow = minbias;
    averagehigh = minbias + averagestep;
    
    %set loop variables
    indexcountlow = 1;
    loopcount = 1; 
    
    %run averaging loop until the end of the data is reached
    while averagehigh <= maxbias
        
        %find the voltage values inside the current averaging range
        Vcurrentaverageinterval = find(VaveragingvectorSORTED >= averagelow & VaveragingvectorSORTED < averagehigh); 
        
        %if there are no values in the voltage step range, do not complete
        %averaging step and move on to the next voltage interval
        if isempty(Vcurrentaverageinterval) == 1
            averagelow = averagehigh;
            averagehigh = averagelow + averagestep;
            
            continue
        else    
        
        %find the corresponding index range so the current values can be
        %averaged
        indexaddition = length(Vcurrentaverageinterval);
        indexcounthigh = indexcountlow + indexaddition - 1;
        %located corresponding current values
        Icurrentaverageinterval = IaveragingvectorSORTED(indexcountlow:indexcounthigh);
        
        %find the average of the current values in the range and set the
        %corresponding voltage value to the middle of the voltage step
        Vaveraged(loopcount, i) = averagelow + averagestep/2;
        Iaveraged(loopcount, i) = mean(Icurrentaverageinterval);
        
        
        %move average step forward
        averagelow = averagehigh;
        averagehigh = averagelow + averagestep;
        
        %adjust loop variables
        indexcountlow = indexcounthigh + 1;
        loopcount = loopcount + 1;
        
        %repeat until the current and data in the interval has been
        %averaged over the bias sweeps
        
        end
        
    end
    
    %repeat this process for every interval in the .h5 file
end

%set results to the output variables
Vin = Vaveraged;
Iout = Iaveraged;


end

