% Have to extract trace info from atlases

% minimal enhancers
atlas_dir = 'atlas20120509.vpc';
atlas_traces;
min_channels = channels;
min_expr = LacZ_expr;

%% fusion enhancers
atlas_dir = 'atlas20120222.vpc';
atlas_traces;
fusion_channels = channels(2:5);
fusion_expr = LacZ_expr(2:5,:);


%% Plot settings

% x-axis index
xaxis = 30:100;


%% Line 207
f3 = figure('Position', [96   306   700   500], 'Resize', 'off', 'Color', 'w');
set(gca, 'FontSize', 12, 'LineWidth', 1)

hold on;

% Line 326
h1 = plot(xaxis,Min_expr(4,xaxis),'-','linewidth',2, 'Color', dark_grey);

% Line 214
h2 = plot(xaxis,Min_expr(2,xaxis),'-','linewidth',2, 'Color', steel);

h3 = plot(xaxis,Fusion_expr(1,xaxis),'-','linewidth',2, 'Color', 'r');

% Add figure info
legend([h1(1), h2(1), h3(1)],...
    {'3+7 minimal', '4+6 minimal', '37F+46F'}, 'Location', 'Best');
xlabel('% Egg Length')
ylabel('Relative intensity')
title('37F+46F');

saveas(gcf, '~/Desktop/line207_trace.eps', 'epsc')


%% Line 208
f3 = figure('Position', [96   306   700   500], 'Resize', 'off', 'Color', 'w');
set(gca, 'FontSize', 12, 'LineWidth', 1)

hold on;

% Line 204
h1 = plot(xaxis,Min_expr(1,xaxis),'-','linewidth',2, 'Color', dark_grey);

% Line 329
h2 = plot(xaxis,Min_expr(7,xaxis),'-','linewidth',2, 'Color', light_grey);

h4 = plot(xaxis,Fusion_expr(2,xaxis),'-','linewidth',2, 'Color', 'r');

% Add figure info
legend([h1(1), h2(1), h4(1)],...
    {'3+7 minimal', '4+6 minimal', '46F+37F'}, 'Location', 'Best');
xlabel('% Egg Length')
ylabel('Relative intensity')
title('46F+37F');

saveas(gcf, '~/Desktop/line208_trace.eps', 'epsc')


%% Line 209
f3 = figure('Position', [96   306   700   500], 'Resize', 'off', 'Color', 'w');
set(gca, 'FontSize', 12, 'LineWidth', 1)

hold on;

% Line 327
h1 = plot(xaxis,Min_expr(5,xaxis),'-','linewidth',2, 'Color', dark_grey);

% Line 214
h2 = plot(xaxis,Min_expr(2,xaxis),'-','linewidth',2, 'Color', light_grey);

h5 = plot(xaxis,Fusion_expr(3,xaxis),'-','linewidth',2, 'Color', 'r');

% Add figure info
legend([h1(1), h2(1), h5(1)],...
    {'3+7 minimal', '4+6 minimal', '37R+46F'}, 'Location', 'Best');
xlabel('% Egg Length');
ylabel('Relative intensity');
title('37R+46F');

saveas(gcf, '~/Desktop/line209_trace.eps', 'epsc')

%% Line 210
f3 = figure('Position', [96   306   700   500], 'Resize', 'off', 'Color', 'w');
set(gca, 'FontSize', 12, 'LineWidth', 1)

hold on;

% Line 326
h1 = plot(xaxis,Min_expr(4,xaxis),'-','linewidth',2, 'Color', dark_grey);

% Line 328
h2 = plot(xaxis,Min_expr(6,xaxis),'-','linewidth',2, 'Color', light_grey);

h6 = plot(xaxis,Fusion_expr(4,xaxis),'-','linewidth',2, 'Color', 'r');

% Add figure info
legend([h1(1), h2(1), h6(1)],...
    {'3+7 minimal', '4+6 minimal', '37F+46R'}, 'Location', 'Best');
xlabel('% Egg Length');
ylabel('Relative intensity');
title('37F+46R');

saveas(gcf, '~/Desktop/line210_trace.eps', 'epsc')

%%  Plot 3+7 traces
f1 = figure();
set(gcf, 'PaperPosition', [0.25 2.5 3 1.5]) % [left_margin bottom_margin width height]
% set(gcf, 'Position', [96   306   700   500]);  
% set(gcf, 'Resize', 'off');
set(gcf, 'Color', 'w');

set(gca, 'FontSize', 12, 'LineWidth', 1)
% set(gca, 'XTick', [0.5 0. 1]); % Sets x-axis to only have ticks in []
% set(gca, 'YTick', [-pi 0 pi]); % Sets y-axis ticks
% set(gca, 'YTickLabel', {'D', 'V', 'D'}); % Sets y-axis tick labels

hold on;

h1 = plot(xaxis, Min_expr(4,xaxis), 'Tag', 'PF');
[x y] = horizontal_errorbar(line204.avg(:,:,cohort), phi2, line204.sem(:,:,cohort), 0.05);
plot(x,y, 'Tag', 'PF')

h2 = plot(line325.avg(:,:,cohort),phi2,'Tag', 'PR');
[x y] = horizontal_errorbar(line325.avg(:,:,cohort), phi2, line325.sem(:,:,cohort), 0.05);
plot(x,y, 'Tag', 'PR')

h3 = plot(line326.avg(:,:,cohort),phi2,'Tag', 'DF');
[x y] = horizontal_errorbar(line326.avg(:,:,cohort), phi2, line326.sem(:,:,cohort), 0.05);
plot(x,y, 'Tag', 'DF')

h4 = plot(line327.avg(:,:,cohort),phi2,'Tag', 'DR');
[x y] = horizontal_errorbar(line327.avg(:,:,cohort), phi2, line327.sem(:,:,cohort), 0.05);
plot(x,y, 'Tag', 'DR')

% For peak and trace plots, standard coloring
set(findobj(gca, 'Type', 'line'), 'LineWidth', 1);
set(findobj(gcf, 'Type', 'line', 'Tag', 'PF'), 'Color', steel);
set(findobj(gcf, 'Type', 'line', 'Tag', 'DF'), 'Color', steel, 'LineStyle', '--');
set(findobj(gcf, 'Type', 'line', 'Tag', 'PR'), 'Color', 'k');
set(findobj(gcf, 'Type', 'line', 'Tag', 'DR'), 'Color', 'k', 'LineStyle', '--');

legend([h1(1), h3(1), h2(1), h4(1)], {'Proximal-Forward', 'Distal-Forward', 'Proximal-Reverse', 'Distal-Reverse'}, 'Location', 'Best');
xlabel('% Egg Length')
ylabel('DV angle (ventral = 0)')
title('eve3+7')
xlim([0.4 1])
ylim([-pi-0.2 pi+0.2])
saveas(gcf, '~/Dropbox/Lydiard-Martin2012/figures/eve37_sem.eps', 'epsc')

