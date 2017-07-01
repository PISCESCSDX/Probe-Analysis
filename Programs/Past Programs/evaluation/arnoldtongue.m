function [le ri] = arnoldtongue(numlist, fdw, bool_fcut)
%20080225-Mo-01:54 Brandt
%function [le ri] = arnoldtongue(numlist, fdw)
% IN: numlist: no. of cou0007.MDF files (e.g. = 7)
%     fdw: self sustained oscillator frequency [Hz]8  (ca!!!)
%     Pt: pulse time trace
%     Ps: pulse signal (for getting the pulse times)
%OUT: arn_le: Al: Al.f, Al.A  -  left Arnold branch
%     arn_ri: Ar: Ar.f, Ar.A  -  right Arnold branch
% EX: [Al Ar] = arnoldtongue([1:11], 4200, [Bt B(:,8)])

if nargin < 3; bool_fcut = 0; end

  a = dir('B1*.MDF');
  b = dir('B2*.MDF');


    % (2) GET EXCITER AMPLITUDE ( =CURRENT)
    %     Get Currents from Current Monitor signals (look at fdw)
    for i=1:length(numlist)
    %-------FIND PULSE TIMES
        [B2 B2t] = readmdf(b(numlist(i)).name);
        [t0(i) t1(i)] = pulse_limits(B2t, B2(:,8), 0.1);
    %
        [Ex Et] = readmdf(a(numlist(i)).name);
    %-------Get Start and End Index
        i_a = find(Et>t0(i)); i_a = i_a(1); i_b = find(Et<t1(i)); i_b = i_b(end);
    %-------MEAN over 8 Rogowski coils
        for j=1:8
            Acoil(j) = std(Ex(i_a:i_b, j));
            [A_out(j), ph_out] = calibration_cm(j, fdw, Acoil(j), 0);
        end
    %-------CURRENT AMPLITUDE for item numlist(i):
        Curr(i) = mean(A_out);
    end


% (3) STEP THROUGH COURONNE PLOTS - COMPARE EXCITER PEAK AND DW PEAK
% (3.1) CALCULATE MEAN PEAK HEIGHT
% (3.2) SMALL WINDOWED 2D-FFT
a = dir('cou*.MDF');
  for i=1:length(numlist)
      savename = ['arnold_dat' num_prezero(i,length(numlist)) '.mat'];
      savename2= ['arnold_kf' num_prezero(i,length(numlist)) '.mat'];
      savename3= ['arnold_syn' num_prezero(i,length(numlist)) '.mat'];
      % TEST WETHER SAVE ALREADY EXISTS
    if ~exist(savename, 'file')
      disp(' ');
      disp(['evaluation of: ' a(numlist(i)).name]);
      %
      %
      [A At] = readmdf(a(numlist(i)).name);
      % Points per window: f_abtast/fdw
      fwinpts = 100*round( (1/fdw)/(At(2)-At(1)) );
      % %-------METHOD A: LOOK ONLY AT f-Spectrum
      %       [tf frq Amat] = tt2spectrogram(tt, sig, winpts)
      %
      %-------METHOD B: LOOK AT kf-Spectrum
      [tvec mvec fvec mmat fmat] = kfspectrogram(At, A', fwinpts);
      save(savename2, 'tvec', 'mvec', 'fvec', 'mmat', 'fmat');
      % load(savename2);

      %PLOT ALSO WITH BOUNDARIES WHERE EXCITER IS OFF
            % PLOT AND PRINT
              figeps(8,10,1)
              subplot(211); subplotfspec(tvec, fvec, 20*log10(fmat)); freezeColors
              subplot(212); subplotkspec(tvec, mvec, mmat');
              print_adv([1 1], '-r300', ['mode-f_' num_prezero(i,length(numlist)) '.eps']);
              clf
      %
      % FIND SYNCHRONIZATION RANGE (tvec mvec fvec mmat fmat)
          % FIND t0 and t1 in fmat and mmat
            ta = find(tvec>t0(i));     i_t0 = ta(1);
            tb = find(tvec<t1(i));     i_t1 = tb(end);
              tind = i_t0:i_t1;
              tvec_cut = tvec(tind);
              mmat_cut = mmat(:, tind);
              fmat_cut = fmat(:, tind);
          % GET EXCITER FREQUENCY - Calculate fsweep: fex(tvec_cut)
            excdat = load('exc.mat');
            Aex   = excdat.cl(numlist(i), 3);
            mex   = excdat.cl(numlist(i), 5);
            fex_0 = excdat.cl(numlist(i), 6);
            fex_1 = excdat.cl(numlist(i), 7);
            fexvec =  fex_0+(0:length(tvec_cut)-1)*(fex_1-fex_0)/(length(tvec_cut)-1);
          % FIND SYNCH LIMITS (fsy_start, fsy_end)
            [f_synch Maratio Aratio] = arnold_fsynch(tvec_cut,mvec,fvec,mmat_cut,fmat_cut,fdw,fexvec,mex);
            save(savename3, 'f_synch', 'Maratio' ,'Aratio');
            %
          % 3. Zusammenhängenden Bereich raussuchen (es kann sein, dass am Rand noch
          %    einmal synchronisiert vorkommt.
          % CUT YES NO ... fcut
          if bool_fcut == 1
            [istart fsy_start iend fsy_end] = arnold_f_synch_cut(f_synch, fdw);
            tvec_sy = tvec_cut([istart iend]);
          end

% PLOT & PRINT - FOR SYNCHRONIZATION-CHECK
  figeps(10,15,1); clf
    subplot(311); subplotfspec(tvec_cut, fvec, 20*log10(fmat_cut)); freezeColors;
    subplot(312); subplotkspec(tvec_cut, mvec, mmat_cut');
    subplot(313); hold on;
      plot(tvec_cut*1000, f_synch./f_synch, 'o'); ylabel('synch', 'fontsize', 12);
    if bool_fcut == 1
        plot(tvec_sy*1000, [0.9 0.9], 'r*'); set(gca, 'ylim', [0 1]);
          set(gca, 'fontsize', 12)
          set(gca, 'xlim', 1e3*[tvec_cut(1) tvec_cut(end)]);
    end
    print_adv([0 1 1], '-r300', ['synch_mode-f_' num_prezero(i,length(numlist)) '.eps']);


% OUTPUT PARAMETERS
    %
    if bool_fcut ~= 1
        nnan = find ( ~isnan(f_synch)==1 );
        fsy_start = f_synch( nnan(1) );
        fsy_end = f_synch( nnan(end) );
    end
    le.f(i) = fsy_start;
    ri.f(i) = fsy_end;
    le.currA(i) = Curr(i);
    ri.currA(i) = Curr(i);
    le.Ae(i) = Aex;
    ri.Ae(i) = Aex;
    le.arnold(i,:) = f_synch;

        % SAVE ARNOLD TOUNGUE PARTS
          save(savename, 'le', 'ri');
    end         % ----- if ~exist
  end           % ----- for

end