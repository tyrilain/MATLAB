%FIGURE_FUSIONS    Makes Figures 3 and 5 of Lydiard-Martin 2013
% 
% This script makes line trace figures for the fusion constructs as well as
% fusion plus spacer constructs
% 
% Tara Martin
% Aug 21, 2012
%


%% Prep traces
% creates LacZ_mean(channel,nucleus,cohort) and LacZ_std containing 
% expression trace and std dev respectively.  
% Also 'channels' contains a key to the construct index.

base_dir = '/Users/tmartin/'; % set manually because changes depending on computer
data_dir = [base_dir 'Desktop/atlas_20130125_freeze/'];
save_dir = [base_dir 'Dropbox/results/'];

%% Load data, extract and reformat relevant channels
load([data_dir 'Data.mat']);
cohorts = [3];

%% get list of genes (only LacZ channels)
channels37 = {'LacZ:0204', 'LacZ:0207','LacZ:0209', 'LacZ:0325','LacZ:0326',...
    'LacZ:0327', 'LacZ:0356', 'LacZ:0323'};
channels46 = {'LacZ:0208','LacZ:0210','LacZ:0214', 'LacZ:0240', 'LacZ:0328',...
    'LacZ:0329','LacZ:0330', 'LacZ:0321', 'LacZ:0322', 'LacZ:0324'};


%% Pull out traces and peak data for 37 dominant lines
channels = channels37;
% extract only relevant time points and channels from scaled data structure
% (assumes only stage 5 embyros and no protein stains)
I = zeros(1,length(gene_names));
for i = 1:length(channels)
    I = I+strcmp(channels{i}, gene_names);
    I = logical(I);
end
channels37 = gene_names(I); % changes order in user defined list to atlas order
atlas_data = Edata(I, cohorts+3,1);
search_range = 40:90; % where to look for stripe peaks (within 1:100)
num_stripes = 2;  % how many stripes to search for
figure_prep_traces2;

%% find peaks on x-axis, calculate normalization factor
for i = 1:length(channels)
    [temp ~] = peakdet(LacZ_mean(i,45:55), 0.1);
%     LacZ_peak{i} = temp;
    norm_factor = max(temp(:,2)); % use stripe 3 peak to normalize
    LacZ_mean(i,:) = LacZ_mean(i,:)/norm_factor;
    LacZ_std(i,:) = LacZ_std(i,:)/norm_factor;
end

temp_mean = LacZ_mean;
temp_std = LacZ_std;
temp_numPC = LacZ_numPC;
% temp_peak = LacZ_peak;

%% Pull out traces and peak data for 46 dominant lines
channels = channels46;
% extract only relevant time points and channels from scaled data structure
% (assumes only stage 5 embyros and no protein stains)
I = zeros(1,length(gene_names));
for i = 1:length(channels)
    I = I+strcmp(channels{i}, gene_names);
    I = logical(I);
end

channels46 = gene_names(I); % changes order in user defined list to atlas order
atlas_data = Edata(I, cohorts+3,1);
search_range = 55:85; % where to look for stripe peaks (within 1:100)
num_stripes = 2;  % how many stripes to search for
figure_prep_traces2;

%% find peaks on x-axis, calculate normalization factor
for i = 1:length(channels)
    [temp ~] = peakdet(LacZ_mean(i,55:65), 0.1);
%     LacZ_peak{i} = temp;
    norm_factor = max(temp(:,2)); % use stripe 4 peak to normalize
    LacZ_mean(i,:) = LacZ_mean(i,:)/norm_factor;
    LacZ_std(i,:) = LacZ_std(i,:)/norm_factor;
end

%% Combine data from all lines
channels = [channels37 channels46];
LacZ_mean = [temp_mean; LacZ_mean];
LacZ_std = [temp_std; LacZ_std];
LacZ_numPC = [temp_numPC; LacZ_numPC];
% LacZ_peak = [temp_peak, LacZ_peak];


% Get index number for relevant lines
line204 = strcmp('LacZ:0204', channels);
line207 = strcmp('LacZ:0207', channels);
line208 = strcmp('LacZ:0208', channels);
line209 = strcmp('LacZ:0209', channels);
line210 = strcmp('LacZ:0210', channels);
line214 = strcmp('LacZ:0214', channels);
line326 = strcmp('LacZ:0326', channels);
line327 = strcmp('LacZ:0327', channels);
line328 = strcmp('LacZ:0328', channels);
line329 = strcmp('LacZ:0329', channels);

%% Set portion of AP axis to use
xx = 1:100;
xaxis = xx/100;


%% Plotting all fusion lines (Fig 3)

% figure_params;
figure('Position', [32 400 415 1200]);
xSize = 8;  ySize = 20;  % Set width and height of figure in cm

% Line 207, Fusion B
subplot(4,1,2)
hold on;
plot(xaxis,LacZ_mean(line326,xx)+LacZ_mean(line214,xx),'g'); % Line 326+214
plot(xaxis,LacZ_mean(line207,xx),'Tag', 'fusion'); % Line 207
plot(xaxis,LacZ_mean(line207,xx)+LacZ_std(line207,xx), 'b', xaxis,LacZ_mean(line207,xx)-LacZ_std(line207,xx), 'b');
ylabel('207')
set_fig;

% Line 208, Fusion D
subplot(4,1,4)
hold on;
plot(xaxis,LacZ_mean(line204,xx)+LacZ_mean(line329,xx),'g'); % Line 204+329
plot(xaxis,LacZ_mean(line208,xx),'Tag', 'fusion'); % Line 208
plot(xaxis,LacZ_mean(line208,xx)+LacZ_std(line208,xx), 'b', xaxis,LacZ_mean(line208,xx)-LacZ_std(line208,xx), 'b');
ylabel('208')
set_fig;

% Line 209, Fusion A
subplot(4,1,1)
hold on;
plot(xaxis,LacZ_mean(line327,xx)+LacZ_mean(line214,xx),'g'); % Line 327+214
plot(xaxis,LacZ_mean(line209,xx),'Tag', 'fusion'); % Line 209
plot(xaxis,LacZ_mean(line209,xx)+LacZ_std(line209,xx), 'b', xaxis,LacZ_mean(line209,xx)-LacZ_std(line209,xx), 'b');
ylabel('209')
set_fig;

% Line 210, Fusion C
subplot(4,1,3)
hold on;
plot(xaxis,LacZ_mean(line326,xx)+LacZ_mean(line328,xx),'g'); % Line 326+328
plot(xaxis,LacZ_mean(line210,xx),'Tag', 'fusion'); % Line 210
plot(xaxis,LacZ_mean(line210,xx)+LacZ_std(line210,xx), 'b', xaxis,LacZ_mean(line210,xx)-LacZ_std(line210,xx), 'b');
ylabel('210')
set_fig;


% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 3]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 0.5 1], 'YTick', [0 1 2 3]);

% Save figure
saveas(gcf, fullfile(save_dir,'fusion_traces.eps'), 'epsc')




%% Fusions with spacers (Figure 5)

figure('Position', [32 400 415 1200]);

% Get index number for relevant lines
line356 = strcmp('LacZ:0356', channels);
line322 = strcmp('LacZ:0322', channels);
line323 = strcmp('LacZ:0323', channels);
line324 = strcmp('LacZ:0324', channels);

% Line 356 (corresponds to 207); Spacer B
subplot(4,1,2)
hold on;
plot(xaxis,LacZ_mean(line207,xx),'g'); % Line 207
plot(xaxis,LacZ_mean(line326,xx)+LacZ_mean(line214,xx),'k'); % Line 326+214
plot(xaxis,LacZ_mean(line356,xx),'Tag', 'fusion');
plot(xaxis,LacZ_mean(line356,xx)+LacZ_std(line356,xx), 'b', xaxis,LacZ_mean(line356,xx)-LacZ_std(line356,xx), 'b');
ylabel('356 (207)')
set_fig;

% Line 322 (corresponds to 208); Spacer D
subplot(4,1,4)
hold on;
plot(xaxis,LacZ_mean(line208,xx),'g'); % Line 208
plot(xaxis,LacZ_mean(line204,xx)+LacZ_mean(line329,xx),'k'); % Line 204+329
plot(xaxis,LacZ_mean(line322,xx),'Tag', 'fusion');
plot(xaxis,LacZ_mean(line322,xx)+LacZ_std(line322,xx), 'b', xaxis,LacZ_mean(line322,xx)-LacZ_std(line322,xx), 'b');
ylabel('322 (208)')
set_fig;

% Line 323 (corresponds to 209); Spacer A
subplot(4,1,1)
hold on;
plot(xaxis,LacZ_mean(line209,xx),'g'); % Line 209
plot(xaxis,LacZ_mean(line327,xx)+LacZ_mean(line214,xx),'k'); % Line 327+214
plot(xaxis,LacZ_mean(line323,xx),'Tag', 'fusion'); 
plot(xaxis,LacZ_mean(line323,xx)+LacZ_std(line323,xx), 'b', xaxis,LacZ_mean(line323,xx)-LacZ_std(line323,xx), 'b');
ylabel('323 (209)')
set_fig;

% Line 324 (corresponds to 210); Spacer C
subplot(4,1,3)
hold on;
plot(xaxis,LacZ_mean(line210,xx),'g'); % Line 210
plot(xaxis,LacZ_mean(line326,xx)+LacZ_mean(line328,xx),'k'); % Line 326+328
plot(xaxis,LacZ_mean(line324,xx),'Tag', 'fusion'); 
plot(xaxis,LacZ_mean(line324,xx)+LacZ_std(line324,xx), 'b', xaxis,LacZ_mean(line324,xx)-LacZ_std(line324,xx), 'b');
ylabel('324 (210)')
set_fig;


% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 3]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 0.5 1], 'YTick', [0 1 2 3]);

% Save figure
saveas(gcf, fullfile(save_dir,'fusion_spacers.eps'), 'epsc')

%% Fusions with 1000bp spacer (Figure S6)

figure('Position', [32 400 415 300]);
xSize = 8;  ySize = 5;  % Set width and height of figure in cm

% Get index number for relevant lines
line240 = strcmp('LacZ:0240', channels);

% Line 240 (corresponds to 207); 1000bp Spacer
hold on;
plot(xaxis,LacZ_mean(line240,xx),'g'); % Line 240
plot(xaxis,LacZ_mean(line326,xx)+LacZ_mean(line214,xx),'k'); % Line 326+214
plot(xaxis,LacZ_mean(line240,xx),'Tag', 'fusion');
plot(xaxis,LacZ_mean(line240,xx)+LacZ_std(line240,xx), 'b', xaxis,LacZ_mean(line240,xx)-LacZ_std(line240,xx), 'b');
ylabel('line240 (207)')
ylim([0 3]);
set(gca, 'XTick', [0 0.5 1])
set_fig;

% Save figure
saveas(gcf, fullfile(save_dir,'fusion_bigspacer.eps'), 'epsc')

%% Fusion E for Angela's talk 20130619

line321 = strcmp('LacZ:0321', channels);

figure('Position', [32 400 415 300]);
xSize = 8;  ySize = 5;  % Set width and height of figure in cm

% line 321 is Fusion E, line 208 is Fusion D (use as comparison)
hold on;
plot(xaxis, LacZ_mean(line208, xx), 'g'); % line 208 for comparison
plot(xaxis, LacZ_mean(line321,xx), 'Tag', 'fusion'); % line 321 mean
plot(xaxis, LacZ_mean(line321,xx)+LacZ_std(line321,xx), 'b', xaxis, LacZ_mean(line321,xx)-LacZ_std(line321,xx), 'b');

ylabel('Relative expression')
ylim([0 3]);
set(gca, 'XTick', [0 0.5 1])
set_fig;

% Save figure
saveas(gcf, fullfile(save_dir,'fusionE_angela.eps'), 'epsc')

