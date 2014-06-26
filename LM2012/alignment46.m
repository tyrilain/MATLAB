%ALIGNMENT46  Makes plot of stripe 4/6 fragments that have been cloned or
%predicted by gemstat-gl

base_dir = '/Users/tlm8/'; % set manually because changes depending on computer
save_dir = [base_dir 'Dropbox/results/'];

% cloned fragments
clyde = [0 1156];
prettygood = [767 1702];
stripe46 = [564 1364];
minimal = [564 1156];
hassan = [564 914];
toosmall = [767 1156];

% gemstat-gl stripe 4/6 windows
offset46 = 564;
win207_46 = [320 720];
win208_46 = [80 730];
win209_46 = [320 730];
win210_46 = [70 770];
win321_46 = [10 600];

figure(1);
plot(clyde, [10 10]);
hold on;
plot(prettygood, [9 9]);
plot(toosmall, [8 8]);
plot(minimal, [7 7]);


plot(stripe46, [6 6]);
plot(hassan, [5 5]);

plot(win207_46+offset46, [2.5 2.5]);
plot(win208_46+offset46, [1.5 1.5]);
plot(win209_46+offset46, [3 3]);
plot(win210_46+offset46, [2 2]);
plot(win321_46+offset46, [1 1]);
hold off;

fontsize = 12;
fontname = 'Arial';
set(gca, 'Box', 'off');
set(gca, 'FontSize', fontsize);
set(gca, 'FontName', fontname);

ylim([0 11]); xlim([-100 1800]);
title('4/6')

set(findobj(gcf, 'Type', 'line'), 'LineWidth', 2);
set(gca, 'YTick', []); % Sets y-axis ticks


saveas(gcf, fullfile(save_dir,'alignment_46_fragments.eps'), 'epsc')
