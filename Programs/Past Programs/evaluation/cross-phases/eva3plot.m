load ev2eva.mat
al = length(Gamma.sig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT BUSINESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD THE FOLNAME (for prints)
  filestring='eex_sydw_';
  filetree = fnametree(pwd);
  folname = [filetree{2}(1:8) filetree{1}];
% SET FONTSIZE
  fonts = 14;
% SET TIME INTERVALL
  tint = 2;
  % CREATE TIME INTERVALL END-INDEX
  t = tt*1e3;
  tend = max(find(t<t(1)+tint))+1;
% SET UPPER FREQUENCY
  f_up = 30e3;
% PARAMETER for Axes-Position
  y0=0.25; yd=0.98;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSPORT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:al
fftvar.f = Gamma.f;
fftvar.A = 20*log10( Gamma.fft{i} );
pdf.x    = Gamma.pdfx{i};
pdf.y    = Gamma.pdfy{i};
pdf.xlim = [-10 10];
pdf.ylim = [0.0001 1];

plotfftpdf(fftvar, pdf);
  print('-depsc2', my_filename(i, 1, 'gamma_', '.eps'));
  close
  %
%  [pdf_x{i} pdf_y{i}] = pdf(100*Gamma{i});
%   ix_lo{i} = find(pdf_x{i}<0);
%   ix_hi{i} = find(pdf_x{i}>0);  
%   int_lo{i} = int_discrete(pdf_x{i}(ix_lo{i}), pdf_y{i}(ix_lo{i}));
%   int_hi{i} = int_discrete(pdf_x{i}(ix_hi{i}), pdf_y{i}(ix_hi{i}));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:al
  % *** INSERT YLIM MANUALLY
  [h1 h2] = plotfftphase(Phase.f{i}, 20*log10( abs(Phase.A{i}) ), ...
  Phase.ph{i}, Phase.phstd{i}, [-170 -60]);
  print('-depsc2', '-r300', my_filename(i, 1, 'cpsd_', '.eps'));
  close
end