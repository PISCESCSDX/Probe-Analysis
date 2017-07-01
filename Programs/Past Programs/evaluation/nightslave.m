pth{1}='/home/data/measured/20071113';
pth{2}='/home/data/measured/20071115';
pth{3}='/home/data/measured/20071116';
pth{4}='/home/data/measured/20071118';

for j=1:length(pth)
  cd(pth{j});
  a=findfolders('BA*.MDF');
  b=findfolders('_exc*');
  for i=1:size(a,2)
    a{i}
    cd([a{i}]);
    mdfzip;
  end;
  for i=1:size(b,2)
    cd([b{i}]),
    exczip;
  end;
end;