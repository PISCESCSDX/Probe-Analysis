function [volatages,ranges] = findOscillationIntervals(filename)
%This program will take a range of input volatages and
%generate a new array containing only the increasing segments of
%the array


%This variable should be adjusted when we get real data
dt = 400;

v = makeivdata();

voltages = zeros(1);
ranges = zeros(0);
first = true;
for i = 1:(size(v)-dt)
    if v(i)<v(i+dt)
        voltages = [voltages v(i)];
        first = true;
    elseif (first == true)
        ranges = [ranges i];
        first = false;
    end
end
plot(voltages,'.');


end

