%PLOTLINETRACE  Plot a line trace with filled area representing error.
%
%   plotlinetrace(xx, line_mean, line_sem, [color])
%
%   xx = index of A/P values to plot
%   line_mean = mean line trace values
%   line_sem = standard error of the mean for each A/P point
%   color = color of line, defaults to black
%
% Tara Martin, 2014-02-11


function plotlinetrace(xx, line_mean, line_sem, line_color)

if nargin <4
    line_color = 'black'; % set default line color if not provided
end

hold_state = ishold;
xx_reverse = xx(end:-1:1);

plot(xx,line_mean(xx), 'Color', line_color, 'LineWidth', 2)
hold on;
fill([xx xx_reverse], [line_mean(xx)-line_sem(xx) ...
    line_mean(xx_reverse)+line_sem(xx_reverse)],... 
    line_color, 'EdgeColor', 'none')

if ~hold_state  % if hold was off, return to that state
    hold off;
end
    

