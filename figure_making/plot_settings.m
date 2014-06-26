% This is for figuring out plot settings

% figure settings
set(gcf, 'Resize', 'off');  % Disables manual re-sizing of figure
set(gcf, 'Color', 'w');  % Sets background color

% Sets size and position of figure display, does not affect printing
set(gcf, 'Position', [96   306   700   500]);  

% Sets print size: [left_margin bottom_margin width height], dimensions 
% should sum appropriately to paper size (ie paperwidth=width+2*left)
set(gcf, 'PaperPosition', [0.25 2.5 8 6]);


% axis settings
set(gca, 'XTick', [0 1]); % Sets x-axis to only have ticks at [0 1]
set(gca, 'YTick', [-pi 0 pi]); % Sets y-axis ticks
set(gca, 'YTickLabel', {'D', 'V', 'D'}); % Sets y-axis tick labels
set(gca, 'FontSize', 12, 'LineWidth', 1)

% This command finds a blue line in the current figure and sets it's color
% to purple and width to 2
set(findobj(gca, 'Type', 'line', 'Color', [0 0 1]), ...
    'Color', [0.5 0 0.5], 'LineWidth', 2);


% For peak and trace plots, standard coloring
set(findobj(gcf, 'Type', 'line'), 'LineWidth', 2);

set(findobj(gcf, 'Type', 'line', 'Tag', 'P-F'), 'Color', steel);
set(findobj(gcf, 'Type', 'line', 'Tag', 'D-F'), 'Color', steel, 'LineStyle', '--');
set(findobj(gcf, 'Type', 'line', 'Tag', 'P-R'), 'Color', 'k');
set(findobj(gcf, 'Type', 'line', 'Tag', 'D-R'), 'Color', 'k', 'LineStyle', '--');
