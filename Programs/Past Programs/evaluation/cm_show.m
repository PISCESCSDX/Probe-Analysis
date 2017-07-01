function cm_show(fname)
%20080225-Mo-00:20 Brandt
%function cm_show(fname)
% Plot all 8 Current Monitor Signals.

[B t] = readmdf(fname);
t=t*1e3;
figeps(8,24,1)

i=0
i=i+1; subplot(811); plot(t, B(:,i)); set(gca, 'xlim', [0 t(end)])
i=i+1; subplot(812); plot(t, B(:,i), 'r'); set(gca, 'xlim', [0 t(end)])
i=i+1; subplot(813); plot(t, B(:,i)); set(gca, 'xlim', [0 t(end)])
i=i+1; subplot(814); plot(t, B(:,i), 'r'); set(gca, 'xlim', [0 t(end)])
i=i+1; subplot(815); plot(t, B(:,i)); set(gca, 'xlim', [0 t(end)])
i=i+1; subplot(816); plot(t, B(:,i), 'r'); set(gca, 'xlim', [0 t(end)])
i=i+1; subplot(817); plot(t, B(:,i)); set(gca, 'xlim', [0 t(end)])
i=i+1; subplot(818); plot(t, B(:,i), 'r'); set(gca, 'xlim', [0 t(end)])

end