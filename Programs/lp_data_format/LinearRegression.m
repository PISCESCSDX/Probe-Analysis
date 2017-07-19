function [m,b] = LinearRegression(x,y)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    m=0;
    b=0;
    sumx=0;
    sumy=0;
    sumxy=0;
    sumx2=0;
    n=size(x,1);
    
    for i=1:n
        sumx=sumx+x(i);
        sumy=sumy+y(i);
        sumxy=sumxy+x(i)*y(i);
        sumx2=sumx2+x(i)*x(i);
    end
    m=(n*sumxy-sumx*sumy)/(n*sumx2-sumx^2);
    b=(sumy-m*sumx)/n;

    
    



end

