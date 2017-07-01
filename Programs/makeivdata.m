function vf = makeivdata()
v = zeros(34500,1);
v(1)=-35;

for i=0:15:34486
    for j=1:10
       v(i+j+1)=v(i+j) + .01;
    end
    for k=1:5
        v(i+k+11)=v(i+k+10) - .01;
    end
   
end
    
for i=34501:34550
    v(i+1)=v(i)-2.3;
end

vf=[v;v;v];    



end

