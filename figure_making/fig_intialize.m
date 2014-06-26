%FIG_INITIALIZE Set figure defaults for x, y axes and figure size
%
%   fig_intialize(fh, xlimits, xticks, ylimits, yticks, [figsize])
%   
%   fh is the figure handle
%   figsize optionally specifies the figure size (usually in cm)
%
% Tara Martin, 2014-03-10


function fig_intialize(fh, xlimits, xticks, ylimits, yticks, figsize)

if nargin<6
    figsize = [0.1 0.1 5 10];
end
set(fh, 'DefaultAxesXLim', xlimits);
set(fh, 'DefaultAxesXTick', xticks);
set(fh, 'DefaultAxesYLim', ylimits);
set(fh, 'DefaultAxesYTick', yticks);
set(fh,'PaperPosition', figsize)


% hold all
% colormap(ct);
% caxis([0 25]);
% set(gcf,'color',background_color, 'InvertHardCopy', paper);
% set(gcf,'PaperUnits','inches')
% set(gca, 'Box', 'Off', 'FontSize', fontsize, 'FontName', fontname,...
%     'LineWidth', 1, 'color', background_color, 'XColor', axis_color,...
%     'YColor', axis_color, 'ZColor', axis_color)
