for i=1:25
  fn=my_filename(i-1,6,'meas_','');
  [fit,par] = loadken(fn);
  n(i) = par.n1;
  
  a=load([fn '.dat']);
  is(i) = mean(a(300:400,2));
end