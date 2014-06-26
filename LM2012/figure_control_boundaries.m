% Making figure comparing minimal elements
data_dir = '~/Desktop/pointclouds';
save_dir = '~/Dropbox/results/';
figure_params;

%% Generate summary.mat files (only once)
SummaryStripes(fullfile(data_dir, 'line204', 'measure'), 'lacZ', 2);
SummaryStripes(fullfile(data_dir, 'line214', 'measure'), 'lacZ', 2);
SummaryStripes(fullfile(data_dir, 'line325', 'measure'), 'lacZ', 2);
SummaryStripes(fullfile(data_dir, 'line326', 'measure'), 'lacZ', 2);
SummaryStripes(fullfile(data_dir, 'line327', 'measure'), 'lacZ', 2);
SummaryStripes(fullfile(data_dir, 'line328', 'measure'), 'lacZ', 2);
SummaryStripes(fullfile(data_dir, 'line329', 'measure'), 'lacZ', 2);
SummaryStripes(fullfile(data_dir, 'line330', 'measure'), 'lacZ', 2);

%% Load summary.mat files
line204 = load(fullfile(data_dir,'line204/measure/summary.mat'), '-mat');
line214 = load(fullfile(data_dir,'line214/measure/summary.mat'), '-mat');
line325 = load(fullfile(data_dir,'line325/measure/summary.mat'), '-mat');
line326 = load(fullfile(data_dir,'line326/measure/summary.mat'), '-mat');
line327 = load(fullfile(data_dir,'line327/measure/summary.mat'), '-mat');
line328 = load(fullfile(data_dir,'line328/measure/summary.mat'), '-mat');
line329 = load(fullfile(data_dir,'line329/measure/summary.mat'), '-mat');
line330 = load(fullfile(data_dir,'line330/measure/summary.mat'), '-mat');

%% Parameter settings
cohort = 3;

% DV angle index
phi = [-pi:(pi/8):(7*pi/8)]';
phi2 = repmat(phi, 1, 2); % use for 2 stripes
phi4 = repmat(phi, 1, 4); % use for 4 stripes

%% Pull out mean and standard error of mean (95% conf interval) for lines
line204.avg = squeeze(line204.boundary_mean(cohort,:,:));
line204.sem = 1.96*squeeze(line204.boundary_std(cohort,:,:))/sqrt(line204.N(cohort));

line214.avg = squeeze(line214.boundary_mean(cohort,:,:));
line214.sem = 1.96*squeeze(line214.boundary_std(cohort,:,:))/sqrt(line214.N(cohort));

line325.avg = squeeze(line325.boundary_mean(cohort,:,:));
line325.sem = 1.96*squeeze(line325.boundary_std(cohort,:,:))/sqrt(line325.N(cohort));

line326.avg = squeeze(line326.boundary_mean(cohort,:,:));
line326.sem = 1.96*squeeze(line326.boundary_std(cohort,:,:))/sqrt(line326.N(cohort));

line327.avg = squeeze(line327.boundary_mean(cohort,:,:));
line327.sem = 1.96*squeeze(line327.boundary_std(cohort,:,:))/sqrt(line327.N(cohort));

line328.avg = squeeze(line328.boundary_mean(cohort,:,:));
line328.sem = 1.96*squeeze(line328.boundary_std(cohort,:,:))/sqrt(line328.N(cohort));

line329.avg = squeeze(line329.boundary_mean(cohort,:,:));
line329.sem = 1.96*squeeze(line329.boundary_std(cohort,:,:))/sqrt(line329.N(cohort));

line330.avg = squeeze(line330.boundary_mean(cohort,:,:));
line330.sem = 1.96*squeeze(line330.boundary_std(cohort,:,:))/sqrt(line330.N(cohort));

%% Plot 3+7 construct boundaries
figure;
hold all;

xSize = 7;  ySize = 7;  % Set width and height of figure in cm

h1 = plot(line204.avg, phi4, 'Tag', 'PF');
[x y] = horizontal_errorbar(line204.avg, phi4, line204.sem, 0.05);
plot(x,y, 'Tag', 'PF')

h2 = plot(line325.avg,phi4,'Tag', 'PR');
[x y] = horizontal_errorbar(line325.avg, phi4, line325.sem, 0.05);
plot(x,y, 'Tag', 'PR')

h3 = plot(line326.avg,phi4,'Tag', 'DF');
[x y] = horizontal_errorbar(line326.avg, phi4, line326.sem, 0.05);
plot(x,y, 'Tag', 'PF') % don't want to dash line error bars

h4 = plot(line327.avg,phi4,'Tag', 'DR');
[x y] = horizontal_errorbar(line327.avg, phi4, line327.sem, 0.05);
plot(x,y, 'Tag', 'PR') % don't want to dash line error bars

set(gca, 'YTick', [-pi 0 pi]); % Sets y-axis ticks
set(gca, 'YTickLabel', {'D', 'V', 'D'}); % Sets y-axis tick labels

set_fig;  % use standard figure settings and line color/type/widths
set(gca, 'Box', 'on');

xlabel('Fraction Egg Length');
xlim([0.4 1]);
ylim([-pi-0.2 pi+0.2]);

% save figure
saveas(gcf, fullfile(save_dir, 'eve37_sem.eps'), 'epsc')

%% Plot 4+6 construct boundaries
figure;
hold all;

xSize = 7;  ySize = 7;  % Set width and height of figure in cm

h1 = plot(line214.avg, phi4, 'Tag', 'PF');
[x y] = horizontal_errorbar(line214.avg, phi4, line214.sem, 0.05);
plot(x,y, 'Tag', 'PF')

h2 = plot(line328.avg,phi4,'Tag', 'PR');
[x y] = horizontal_errorbar(line328.avg, phi4, line328.sem, 0.05);
plot(x,y, 'Tag', 'PR')

h3 = plot(line329.avg,phi4,'Tag', 'DF');
[x y] = horizontal_errorbar(line329.avg, phi4, line329.sem, 0.05);
plot(x,y, 'Tag', 'PF') % don't want to dash line error bars

h4 = plot(line330.avg,phi4,'Tag', 'DR');
[x y] = horizontal_errorbar(line330.avg, phi4, line330.sem, 0.05);
plot(x,y, 'Tag', 'PR') % don't want to dash line error bars

set(gca, 'YTick', [-pi 0 pi]); % Sets y-axis ticks
set(gca, 'YTickLabel', {'D', 'V', 'D'}); % Sets y-axis tick labels

set_fig;  % use standard figure settings and line color/type/widths
set(gca, 'Box', 'on');

xlabel('Fraction Egg Length');
xlim([0.4 1]);
ylim([-pi-0.2 pi+0.2]);

% save figure
saveas(gcf, fullfile(save_dir, 'eve46_sem.eps'), 'epsc')

