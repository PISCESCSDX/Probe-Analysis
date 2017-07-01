function impdat
% IMPDAT: Important dates.
% Fill in manually.

j = 0;

j = j+1;
event(j).name = 'Pregnant days (doctor)     ';
event(j).date = '2013/04/27 08:00:00.000.000.000';

j = j+1;
event(j).name = 'Print EPS-Poster           ';
event(j).date = '2013/06/28 16:00:00.000.000.000';

j = j+1;
event(j).name = 'Flight SAN -> HEL          ';
event(j).date = '2013/06/29 07:42:00.000.000.000';

j = j+1;
event(j).name = 'Flight HEL -> HAM          ';
event(j).date = '2013/07/06 06:05:00.000.000.000';

j = j+1;
event(j).name = 'Flight HAM -> SAN          ';
event(j).date = '2013/07/23 09:00:00.000.000.000';

j = j+1;
event(j).name = 'Submit Paper to PSST       ';
event(j).date = '2013/08/31 23:00:00.000.000.000';

t1 = clock_int;
le = length(event);

disp(' ')
disp('Upcoming Events:')
disp('===================================================')
for i = 1:le
  disp([event(i).name ':   ' clockdiff(t1, event(i).date)]);
end
disp('===================================================')


df1 = clockdiff(t1, event(1).date);
df2 = clockdiff(t1, event(2).date);
sek1 = str2num(df1(18:19));
min1 = str2num(df1(12:13));
hou1 = str2num(df1( 8: 9));
day1 = str2num(df1( 1: 5));
switch sign(day1)
  case -1
    t1 = -1*(sek1+60*min1+3600*hou1+86400*abs(day1));
  case 0
    t1 = sek1+60*min1+3600*hou1+86400*abs(day1);
  case +1
    t1 = sek1+60*min1+3600*hou1+86400*abs(day1);
end
sek2 = str2num(df2(18:19));
min2 = str2num(df2(12:13));
hou2 = str2num(df2( 8: 9));
day2 = str2num(df2( 1: 5));
switch sign(day2)
  case -1
    t2 = -1*(sek2+60*min2+3600*hou2+86400*abs(day2));
  case 0
    t2 = sek2+60*min2+3600*hou2+86400*abs(day2);
  case +1
    t2 = sek2+60*min2+3600*hou2+86400*abs(day2);
end
dtpercent = abs(t1)/(t2-t1)*100;
disp(['verstrichene Zeit          :   ' num2str(dtpercent) ' %'])

end