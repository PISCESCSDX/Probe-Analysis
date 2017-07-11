function [Vavg,Iavg] = formatAndAverageLangmuirData(filebase)
%This program will take a range of input volatages and
%generate a new array containing only the increasing segments of
%the array (eliminates sharp turns)


A=3000;
v=h5read([filebase '.h5'],'/V');
v = reshape(v,[size(v,1)*size(v,2),1]);
I=h5read([filebase '.h5'],'/I');
I = reshape(I,[size(I,1)*size(I,2),1]);

%Read in the current noise and subtract all values in I by the noise



voltages = zeros(1);
currents = zeros(1);
first = true;
j =1;
k=1;
i=1;
while i<(length(v)-2*A)
    sumb=0;
    suma=0;
    for l=i:i+A
        sumb=sumb+v(l);
        suma=suma+v(l+A);
    end
    if sumb<suma-A/75
        voltages(k,j) = v(i);
        currents(k,j) = I(i);
        k=k+1;
        first=true;
    elseif(first==true)
        k = 1;
        j=j+1;
        first=false;
    end
    i=i+1;
end
voltages(:,1)=[];
currents(:,1)=[];
voltages(:,size(voltages,2))=[];
currents(:,size(currents,2))=[];
voltages(4203:size(voltages,1),:)=[];
currents(4203:size(currents,1),:)=[];

Vavg=zeros(size(voltages,1),floor(size(voltages,2)/20));
Iavg=zeros(size(voltages,1),floor(size(currents,2)/20));

for i=1:length(Vavg(1,:))
    for j=1:length(Vavg(:,1))
        for k=1:20
            Vavg(j,i)= Vavg(j,i)+ voltages(j,20*(i-1)+k);
            Iavg(j,i)= Iavg(j,i)+ currents(j,20*(i-1)+k);
        end
        Vavg(j,i) = Vavg(j,i)/20;
        Iavg(j,i) = Iavg(j,i)/20;
    end
end

end



