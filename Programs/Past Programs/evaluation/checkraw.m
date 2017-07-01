function checkraw
%CHECKRAW checks couronne data and Board2 data from Frisch.
%EPS-pictures will be saved.

a=findfiles('cou*MDF');  al=length(a);
b=findfiles('B2*MDF'); bl=length(b);

% # INPUT: Set END FREQUENCY of the plots
  fend = 30e3;

% LOAD and PLOT COURONNE DATA
for i=1:al
  [fpath fname] = fnamesplit( cell2mat(a(i)) );
  cd(fpath);
  [AA At] = readmdf(fname);
  [Af, Aa, Aph] = su_fft(At, AA(:,3));
  [A.f, A.a] = fspcut(Af,Aa,fend);
  % FIGURE COURONNE
  figeps(15,10,1);
  clf
  plot(A.f/1e3,20*log10(A.a));
  xlabel('f [kHz]'); ylabel('S [dB]');
  print('-depsc2', [ 'cou_' num2str(i) '.eps']);
end

% LOAD B2*.MDF data
for i=1:bl
  disp_num(i,bl);
  [fpath fname] = fnamesplit( cell2mat(b(i)) );
  cd(fpath);
  [BB Bt] = readmdf(fname);
  [Bfrq, Bamp, Bph] = su_fft(Bt, BB(:,2));
  [Bf{i} Ba{i}] = fspcut(Bfrq,Bamp,fend);
end
B.a = cell2mat(Ba);
B.f  = cell2mat(Bf);
B.f =  B.f(:,1)/1e3;
B.shot=(1:bl);
  % FIGURE COURONNE
  figeps(15,10,1);
  clf
  plot(Bf{100}, 20*log10(Ba{100}));
  xlabel('f [kHz]'); ylabel('S [dB]');
  print('-depsc2', [ 'B2.eps']);
% FIGURE 3
  figeps(12,8,1);
  clf
  axes('position', [0.17 0.2 0.75 0.7]);
  pcolor(B.shot, B.f, 20*log10(B.a)); shading flat
  xlabel('shot #');
  ylabel('f [kHz');
  print_adv([1], '-r200', 'spectrogram_B2ch2.eps', 100);

save cbmain.mat A B

end