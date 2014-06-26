%FIGURE_CONTROLS    Makes Figures 2 of Lydiard-Martin 2013
% 
% This script makes line trace figures for the control constructs
% 
% Tara Martin
% Nov 6, 2012
%

%% Prep traces
% creates LacZ_mean(channel,nucleus,cohort) and LacZ_std containing 
% expression trace and std dev respectively.  
% Also 'channels' contains a key to the construct index.

base_dir = '/Users/tlm8/'; % set manually because changes depending on computer
data_dir = [base_dir 'Dropbox/data/atlas_20120810_gemstat/'];
save_dir = [base_dir 'Dropbox/results/'];

%% Load data, extract and reformat relevant channels
load([data_dir 'Data.mat']);

% get list of genes (only LacZ channels)
channels = {'LacZ:0204', 'LacZ:0325', 'LacZ:0326', 'LacZ:0327',...
    'LacZ:0214', 'LacZ:0328', 'LacZ:0329', 'LacZ:0330'};

% extract only relevant time points and channels from scaled data structure
% (assumes only stage 5 embyros and no protein stains)
I = zeros(1,length(gene_names));
for i = 1:length(channels)
    I = I+strcmp(channels{i}, gene_names);
    I = logical(I);
end

channels = gene_names(I); % changes order in user defined list to atlas order

cohorts = [3];
atlas_data = Edata(I, cohorts+3,1);

search_range = 40:90; % where to look for stripe peaks (within 1:100)
num_stripes = 2;  % how many stripes to search for

%% extract line traces and peak heights
figure_prep_traces2;

% Get index number for relevant lines
line204 = strcmp('LacZ:0204', channels);
line325 = strcmp('LacZ:0325', channels);
line326 = strcmp('LacZ:0326', channels);
line327 = strcmp('LacZ:0327', channels);
line214 = strcmp('LacZ:0214', channels);
line328 = strcmp('LacZ:0328', channels);
line329 = strcmp('LacZ:0329', channels);
line330 = strcmp('LacZ:0330', channels);

% Set portion of AP axis to use
xx = 1:100;
xaxis = xx/100;

%% find peaks on x-axis, calculate std from individual traces
for i = 1:length(channels)
    [temp ~] = peakdet(LacZ_mean(i,xx), 0.1);
    stripe_coords(i,:) = temp(:,1);
    norm_factor(i) = temp(1,2);
    stripe_peak(i,:) = temp(:,2)./norm_factor(i);
    stripe_std(i,:) = std(LacZ_peaks{i})./norm_factor(i);
end


%% Make 3+7 trace figure
% Note: all traces are normalized to either stripe 3 or stripe 4 peak value
figure;
hold all;
xSize = 9;  ySize = 6;  % Set width and height of figure in cm


% Line 204
% plot normalized trace
plot(xaxis,LacZ_mean(line204,xx)./norm_factor(line204),'Tag', 'PF');
errorbar(stripe_coords(line204,:)/100, stripe_peak(line204,:), stripe_std(line204,:), '.','Tag', 'PF');

% Line 325
plot(xaxis,LacZ_mean(line325,xx)./norm_factor(line325),'Tag', 'PR');
errorbar(stripe_coords(line325,:)/100, stripe_peak(line325,:), stripe_std(line325,:), '.','Tag', 'PR');

% Line 326
plot(xaxis,LacZ_mean(line326,xx)./norm_factor(line326),'Tag', 'DF');
errorbar(stripe_coords(line326,:)/100, stripe_peak(line326,:), stripe_std(line326,:), '.','Tag', 'DF');

% Line 327
plot(xaxis,LacZ_mean(line327,xx)./norm_factor(line327),'Tag', 'DR');
errorbar(stripe_coords(line327,:)/100, stripe_peak(line327,:), stripe_std(line327,:), '.','Tag', 'DR');


set_fig;  % use standard figure settings and line color/type/widths

xlabel('Fraction Egg Length');
ylabel('Normalized expression');
xlim([0 1]);  
% ylim([0 1]);

% save figure
saveas(gcf, fullfile(save_dir, 'eve37_trace.eps'), 'epsc')

%% Make 4+6 trace figure
figure;
hold all;
% xSize = 7;  ySize = 4;  % Set width and height of figure in cm

% Line 214
plot(xaxis,LacZ_mean(line214,xx)./norm_factor(line214),'Tag', 'PF');
errorbar(stripe_coords(line214,:)/100, stripe_peak(line214,:), stripe_std(line214,:), '.','Tag', 'PF');

% Line 328
plot(xaxis,LacZ_mean(line328,xx)./norm_factor(line328),'Tag', 'PR');
errorbar(stripe_coords(line328,:)/100, stripe_peak(line328,:), stripe_std(line328,:), '.','Tag', 'PR');

% Line 329
plot(xaxis,LacZ_mean(line329,xx)./norm_factor(line329),'Tag', 'DF');
errorbar(stripe_coords(line329,:)/100, stripe_peak(line329,:), stripe_std(line329,:), '.','Tag', 'DF');

% Line 330
plot(xaxis,LacZ_mean(line330,xx)./norm_factor(line330),'Tag', 'DR');
errorbar(stripe_coords(line330,:)/100, stripe_peak(line330,:), stripe_std(line330,:), '.','Tag', 'DR');


set_fig;  % use standard figure settings and line color/type/widths

xlabel('Fraction Egg Length');
ylabel('Normalized expression');
xlim([0 1]);  
% ylim([0 3.5]);

% save figure
saveas(gcf, fullfile(save_dir, 'eve46_trace.eps'), 'epsc')


