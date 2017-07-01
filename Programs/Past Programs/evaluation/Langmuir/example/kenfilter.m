% TRY TO REMOVE THE BAD OSCILLATION ON THE CHARACTERISTIC

a=load('meas_000010.dat');
ken = a(:,2);

Wp=[0.1 0.2];
Ws=[0.12 0.18];
Rp=3; Rs=30;

[N,Wn] = BUTTORD(Wp, Ws, Rp, Rs);
[B,A] = butter(N,Wn,'stop');
y=filter(B,A,ken);

hold on, plot(y,'r'), plot(ken,'k'), hold off