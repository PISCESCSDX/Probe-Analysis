function [] = plotkf(freq, mvec, kfmat, kfaxis, z_max);
%function [] = plotkf(freq, mvec, kfmat, kfaxis);
%
% m-files needed:   subplotkf
% input:    freq    frequency vector
%           mvec    m-vector
%           kfmat   2d kf-data 
%           kfaxis  [fa fb ma mb 0 kfmax]
% output:   plot of the kf-spectrum
% EXAMPLE:  plotkf(freq, mvec, kfmat, kfaxis);

fonts = 14;

if nargin<5; z_max = 0; end;
if nargin<4; error('Input arguments are missing!'); end;


set(gcf,'PaperUnits','centimeters','PaperType','a4letter',...
        'PaperPosition',[1 200 12 8], 'Color', [1.0 1.0 1.0]);
%-- make pictures look the same on printer and on screen
    wysiwyg;

    axes('Position', [0.15 0.12 0.7 0.90]);

    subplotkf(freq, mvec, kfmat, kfaxis, z_max);
    
end