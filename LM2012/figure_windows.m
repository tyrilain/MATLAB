%FIGURE_WINDOWS    Makes plots of gemstat-gl windows for Lydiard-Martin 2013
% 
% This script plots lines representing the fusion constructs and windows
% used for gemstat-gl.
% 
% Tara Martin
% Jan 3, 2013
%

base_dir = '/Users/tmartin/'; % set manually because changes depending on computer
save_dir = [base_dir 'Dropbox/results/'];

color37 = 'r';
color46 = 'b';

% figure_params;
figure('Position', [32 400 415 1200]);
xSize = 8;  ySize = 20;  % Set width and height of figure in cm

% Line 207 (tlm01)
subplot(5,1,2)
plot([0,510], [2 2], color37, [511, 511+800], [2 2], color46);
hold on;
plot([0,0+390], [1 1], color37, [511+320, 511+320+400], [1 1], color46);
set_fig;

% Line 208 (tlm02)
subplot(5,1,4)
plot([0,800], [2 2], color46, [801, 511+800], [2 2], color37);
hold on;
plot([80,80+650], [1 1], color46, [801+240, 801+240+270], [1 1], color37);
set_fig;

% Line 209 (tlm03)
subplot(5,1,1)
plot([0,510], [2 2], color37, [511, 511+800], [2 2], color46);
hold on;
plot([30,30+440], [1 1], color37, [511+320, 511+320+410], [1 1], color46);
set_fig;

% Line 210 (tlm04)
subplot(5,1,3)
plot([0,510], [2 2], color37, [511, 511+800], [2 2], color46);
hold on;
plot([0,0+500], [1 1], color37, [511+30, 511+30+700], [1 1], color46);
set_fig;

% Line 321 (tlm17)
subplot(5,1,5)
plot([0,800], [2 2], color46, [800, 511+800], [2 2], color37);
hold on;
plot([800,800+490], [1 1], color37, [200, 200+590], [1 1], color46);
ylabel('tlm17')
set_fig;

% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1311], 'YLim', [0 3]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 510 800 1311], 'YTick', []);


%% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_windows.eps'), 'epsc')


%% Making alignment of windows to display with heatmaps

% stripe 3/7 windows
win207_37 = [0 390];
win208_37 = [240 510];
win209_37 = [40 480];
win210_37 = [0 500];
win321_37 = [20 510];

figure(1);
plot(win207_37, [3 3], win208_37, [1 1], win209_37, [4 4], win210_37, [2 2], win321_37, [0.5 0.5]);
ylim([0 5]); xlim([0 510]);
title('3/7')
saveas(gcf, fullfile(save_dir,'gemstat_windows37.eps'), 'epsc')

% stripe 4/6 windows
win207_46 = [320 720];
win208_46 = [80 730];
win209_46 = [320 730];
win210_46 = [70 770];
win321_46 = [10 600];

figure(2);
plot(win207_46, [3 3], win208_46, [1 1], win209_46, [4 4], win210_46, [2 2], win321_46, [0.5 0.5]);
ylim([0 5]); xlim([0 800]);
title('4/6')
saveas(gcf, fullfile(save_dir,'gemstat_windows46.eps'), 'epsc')


