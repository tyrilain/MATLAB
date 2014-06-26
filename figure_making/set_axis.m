%set axis properties

% set(gca, 'Box', 'Off', 'FontSize', fontsize, 'FontName', fontname,...
%     'LineWidth', 2, 'color', background_color, 'XColor', axis_color,...
%     'YColor', axis_color, 'ZColor', axis_color);

set(gca, 'XLim', [15 90], 'YLim', ylimits);
set(gca, 'XTick', [15 50 90], 'YTick', yticks);