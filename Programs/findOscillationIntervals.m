function [voltages,ranges] = findOscillationIntervals(filename)
%This program will take a range of input volatages and
%generate a new array containing only the increasing segments of
%the array,


%This variable should be adjusted when we get real data
dt = 5000;
v=h5read(filename,'/V');
v = reshape(v,[size(v,1)*size(v,2),1]);

voltages = zeros(1);
ranges = zeros(0);
first = true;
j =1;
k=1;
for i = 1+dt:(size(v))
    if (v(i)<v(i-dt))
        voltages(k,j) = v(i);
        first = true;
        k=k+1;
    elseif (first == true)
        ranges = [ranges i];
        first = false;
        k = 1;
        j=j+1;
    end
end
plot(voltages,'.');
end



