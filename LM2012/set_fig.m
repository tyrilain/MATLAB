% This script specifies some variables common to all figures

fontsize = 10;
fontname = 'Arial';

% commonly used colors
wh = [1 1 1];
bl = [0 0 0];
dark_grey = [0.4 0.4 0.4];
light_grey = [0.7 0.7 0.7];
steel = [70 130 180]/255;
dark_blue = [35 65 90]/255;

% set_fig used to set figure settings
set(findobj(gca, 'Type', 'line'), 'LineWidth', 1);

set(findobj(gcf, 'Tag', 'PF'), 'Color', steel);

set(findobj(gcf, 'Tag', 'DF'), 'Color', dark_blue);
% set(findobj(gcf, 'Type', 'line', 'Tag', 'DF'), 'LineStyle', '--');
% temp = get(findobj(gcf, 'Tag', 'DF', 'Type', 'hggroup'), 'Children');
% set(temp(2), 'LineStyle', '--');

set(findobj(gcf, 'Tag', 'PR'), 'Color', light_grey);

set(findobj(gcf, 'Tag', 'DR'), 'Color', bl);
% set(findobj(gcf, 'Type', 'line', 'Tag', 'DR'), 'LineStyle', '--');
% temp = get(findobj(gcf, 'Tag', 'DR', 'Type', 'hggroup'), 'Children');
% set(temp(2), 'LineStyle', '--');

set(findobj(gcf, 'Tag', 'fusion'), 'Color', 'r');



set(gca, 'Box', 'off');
set(gca, 'FontSize', fontsize);
set(gca, 'FontName', fontname);
set(gca, 'LineWidth', 1);
%set(gca, 'YTick', [0 1]);


% xSize = 5.5; ySize = 6.3; xLeft = 0.5; yTop = 0.5;
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0.5 0.5 xSize ySize]);
set(gcf, 'Color', 'w');



% Other useful properties that can be set -----------------------
% 
% set(gca, 'XTick', [0 1]);
% set(gca, 'XTickLabel', ['A'; 'P']);
% set(gca, 'YTick', [0 1]);
% set(gca, 'YTickLabel', ['V' ;'D']);
% 
