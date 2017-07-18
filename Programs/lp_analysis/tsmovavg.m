function output = tsmovavg(vector,form,numperiod,dim)
%Moving average function

%The first if is just to match the actual function
if(form=='t')
    vectorh=vector;
    size = length(vectorh);
%dim determines the orientation of the output
    if (dim==2)
        output = zeros(1,size);
    elseif (dim==1)
        output = zeros(size, 1);
    end
    halfnumperiod=ceil((numperiod+1)/2);
%Calculating values
    for i=1:size
%For finding low values
        if (i<halfnumperiod)
            sum=0;
            for(j=1:2*i-1)
                sum=sum+vectorh(j);
            end
            output(i)=sum/(2*i-1);   
%For finding high values
        elseif(i>size-halfnumperiod+1)
            sum=0;
            for(j=2*i-size:size)
                sum=sum+vectorh(j);
            end
            output(i)=sum/(2*(size-i)+1);
%For finding middle values
        else
            if(mod(numperiod, 2)==0)
                numperiod=numperiod+1;
            end
            sum=0;
            for(j=i-halfnumperiod+1:i+halfnumperiod-1)
               sum=sum+vectorh(j);
            end
            output(i)=sum/numperiod;
        end
    end
%Do it again
    vectorh=output;
    for i=1:size
        if (i<halfnumperiod)
            sum=0;
            for(j=1:2*i-1)
                sum=sum+vectorh(j);
            end
            output(i)=sum/(2*i-1);   
        elseif(i>size-halfnumperiod+1)
            sum=0;
            for(j=2*i-size:size)
                sum=sum+vectorh(j);
            end
            output(i)=sum/(2*(size-i)+1);
        else
            if(mod(numperiod, 2)==0)
                numperiod=numperiod-1;
            end
            sum=0;
            for(j=i-halfnumperiod+1:i+halfnumperiod-1)
               sum=sum+vectorh(j);
            end
            output(i)=sum/numperiod;
        end
    end
end
end

