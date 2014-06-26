%FIGURE_GEMSTAT    Makes Gemstat Figures of Lydiard-Martin 2013
% 
% This script makes line trace figures for the gemstat modeling results
% 
% Tara Martin
% Jan 2, 2013
%

%% Initial stuff
% Pre-req for this script:
% Manually import data from model-predictions.xlsx and save to the variables
% "gemstat", "twogrf", and "gemstat_gl"; also save gemstat_gl row names as
% "gemstat_gl_rows" cell array for reference.

base_dir = '/Users/tmartin/'; % set manually because changes depending on computer
data_dir = [base_dir 'Dropbox/data/atlas_20120810_gemstat/'];
save_dir = [base_dir 'Dropbox/results/'];

%% Load data, extract and reformat relevant channels
load([data_dir 'Data.mat']);
cohorts = [3];

%% get list of genes (only LacZ channels)
channels37 = {'LacZ:0204', 'LacZ:0207', 'LacZ:0209', 'LacZ:0321', ...
    'LacZ:0325', 'LacZ:0326', 'LacZ:0327', 'LacZ:0356', 'LacZ:0323'};
channels46 = {'LacZ:0208', 'LacZ:0210', 'LacZ:0214', 'LacZ:0240', ...
    'LacZ:0328', 'LacZ:0329', 'LacZ:0330', 'LacZ:0322', 'LacZ:0324'};


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

temp_mean = LacZ_mean;

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

%% Combine data from all lines
channels = [channels37 channels46];
LacZ_mean = [temp_mean; LacZ_mean];

%% Get index number for relevant lines
line204 = strcmp('LacZ:0204', channels);
line214 = strcmp('LacZ:0214', channels);

line207 = strcmp('LacZ:0207', channels);
line208 = strcmp('LacZ:0208', channels);
line209 = strcmp('LacZ:0209', channels);
line210 = strcmp('LacZ:0210', channels);
line321 = strcmp('LacZ:0321', channels);

line356 = strcmp('LacZ:0356', channels);
line322 = strcmp('LacZ:0322', channels);
line323 = strcmp('LacZ:0323', channels);
line324 = strcmp('LacZ:0324', channels);

line240 = strcmp('LacZ:0240', channels);

%% Set portion of AP axis to use
xx = 1:100;
xaxis = xx/100;

%% normalize to 1 within each channel
norm_mat = max(LacZ_mean,[],2);
% LacZ_mean = LacZ_mean./repmat(norm_mat,1,100);

% also normalize gemstat fits to same factors as measured
% norm_mat_gem = [norm_mat(line207); ...
%     norm_mat(line208);...
%     norm_mat(line209);...
%     norm_mat(line210);...
%     norm_mat(line356);...
%     norm_mat(line240);...
%     norm_mat(line322);...
%     norm_mat(line323);...
%     norm_mat(line324);...
%     max(gemstat(10,:));...
%     max(gemstat(11,:))];

%% actually normalizing to 1
gemstat = gemstat_tmp./repmat(max(gemstat_tmp,[],2),1,100);
twogrf = twogrf_tmp./repmat(max(twogrf_tmp,[],2),1,100);


%% Plotting gemstat fits to components (done separately from the rest so variable names are weird)

% figure_params;
figure('Position', [32 400 415 1200]);
xSize = 8;  ySize = 10;  % Set width and height of figure in cm

subplot(2,1,1)
hold on;
plot(xaxis,gemstat(10,:), 'r');
plot(xaxis, LacZ_mean(line204,xx), 'b');
ylabel('eve37')
set_fig;

subplot(2,1,2)
hold on;
plot(xaxis,gemstat(11,:), 'r');
plot(xaxis, LacZ_mean(line214,xx), 'b');
ylabel('eve46')
set_fig;

% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 1]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 0.5 1], 'YTick', [0 1]);

%% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_components.eps'), 'epsc')


%% Plotting single GRF fits to fusions

% figure_params;
figure('Position', [32 400 415 1200]);
xSize = 8;  ySize = 20;  % Set width and height of figure in cm

% Line 207 (tlm01)
subplot(4,1,2)
hold on;
plot(xaxis,gemstat(1,:), 'r');
plot(xaxis, LacZ_mean(line207,xx), 'b');
ylabel('tlm01')
set_fig;

% Line 208 (tlm02)
subplot(4,1,4)
hold on;
plot(xaxis,gemstat(2,:), 'r');
plot(xaxis, LacZ_mean(line208,xx), 'b');
ylabel('tlm02')
set_fig;

% Line 209 (tlm03)
subplot(4,1,1)
hold on;
plot(xaxis,gemstat(3,:), 'r');
plot(xaxis, LacZ_mean(line209,xx), 'b');
ylabel('tlm03')
set_fig;

% Line 210 (tlm04)
subplot(4,1,3)
hold on;
plot(xaxis,gemstat(4,:), 'r');
plot(xaxis, LacZ_mean(line210,xx), 'b');
ylabel('tlm04')
set_fig;


% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 1]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 0.5 1], 'YTick', [0 1]);

% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_single.eps'), 'epsc')


%% Plotting double GRF fits to fusions

% figure_params;
figure('Position', [32 400 415 1200]);
xSize = 8;  ySize = 20;  % Set width and height of figure in cm

% Line 207 (tlm01)
subplot(4,1,2)
hold on;
plot(xaxis,twogrf(1,:), 'r');
plot(xaxis, LacZ_mean(line207,xx), 'b');
ylabel('tlm01')
set_fig;

% Line 208 (tlm02)
subplot(4,1,4)
hold on;
plot(xaxis,twogrf(2,:), 'r');
plot(xaxis, LacZ_mean(line208,xx), 'b');
ylabel('tlm02')
set_fig;

% Line 209 (tlm03)
subplot(4,1,1)
hold on;
plot(xaxis,twogrf(3,:), 'r');
plot(xaxis, LacZ_mean(line209,xx), 'b');
ylabel('tlm03')
set_fig;

% Line 210 (tlm04)
subplot(4,1,3)
hold on;
plot(xaxis,twogrf(4,:), 'r');
plot(xaxis, LacZ_mean(line210,xx), 'b');
ylabel('tlm04')
set_fig;


% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 1]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 0.5 1], 'YTick', [0 1]);

% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_twogrf.eps'), 'epsc')


%% Plotting GEMSTAT-GL window fits to fusions

% figure_params;
figure('Position', [32 400 415 1200]);
xSize = 8;  ySize = 20;  % Set width and height of figure in cm

% Line 207 (tlm01)
subplot(4,1,2)
hold on;
plot(xaxis, LacZ_mean(line207,xx), 'b');
plot(xaxis, gemstat_gl(1,:)./norm_mat(line207), 'r', xaxis, gemstat_gl(2,:)./norm_mat(line207), 'g');
ylabel('tlm01')
set_fig;

% Line 208 (tlm02)
subplot(4,1,4)
hold on;
plot(xaxis, LacZ_mean(line208,xx), 'b');
plot(xaxis, gemstat_gl(4,:)./norm_mat(line208), 'r', xaxis, gemstat_gl(5,:)./norm_mat(line208), 'g');
ylabel('tlm02')
set_fig;

% Line 209 (tlm03)
subplot(4,1,1)
hold on;
plot(xaxis, LacZ_mean(line209,xx), 'b');
plot(xaxis, gemstat_gl(7,:)./norm_mat(line209), 'r', xaxis, gemstat_gl(8,:)./norm_mat(line209), 'g');
ylabel('tlm03')
set_fig;

% Line 210 (tlm04)
subplot(4,1,3)
hold on;
plot(xaxis, LacZ_mean(line210,xx), 'b');
plot(xaxis, gemstat_gl(10,:)./norm_mat(line210), 'r', xaxis, gemstat_gl(11,:)./norm_mat(line210), 'g');
ylabel('tlm04')
set_fig;


% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 1.1]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 1], 'YTick', [0 1]);

% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_GL1.eps'), 'epsc')

%% Plotting GEMSTAT-GL fits to fusions

% figure_params;
figure('Position', [32 400 415 1200]);
xSize = 8;  ySize = 20;  % Set width and height of figure in cm

% Line 207 (tlm01)
subplot(4,1,2)
hold on;
plot(xaxis, LacZ_mean(line207,xx), 'b');
plot(xaxis, gemstat_gl(3,:)/max(gemstat_gl(3,:)), 'r');
ylabel('tlm01')
set_fig;

% Line 208 (tlm02)
subplot(4,1,4)
hold on;
plot(xaxis, LacZ_mean(line208,xx), 'b');
plot(xaxis, gemstat_gl(6,:)/max(gemstat_gl(6,:)), 'r');
ylabel('tlm02')
set_fig;

% Line 209 (tlm03)
subplot(4,1,1)
hold on;
plot(xaxis, LacZ_mean(line209,xx), 'b');
plot(xaxis, gemstat_gl(9,:)/max(gemstat_gl(9,:)), 'r');
ylabel('tlm03')
set_fig;

% Line 210 (tlm04)
subplot(4,1,3)
hold on;
plot(xaxis, LacZ_mean(line210,xx), 'b');
plot(xaxis, gemstat_gl(12,:)/max(gemstat_gl(12,:)), 'r');
ylabel('tlm04')
set_fig;


% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 1]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 0.5 1], 'YTick', [0 1]);

%% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_GL2.eps'), 'epsc')



%% -----------SPACERS-----------------------------

% Plotting single GRF fits to spacers

% figure_params;
figure('Position', [32 400 415 1200]);
xSize = 8;  ySize = 20;  % Set width and height of figure in cm

% Line 356 (tlm07; spacer B)
subplot(5,1,2)
hold on;
plot(xaxis,gemstat(5,:), 'r');
plot(xaxis, LacZ_mean(line356,xx), 'b');
ylabel('tlm07')
set_fig;

% Line 322 (tlm18; spacer D)
subplot(5,1,4)
hold on;
plot(xaxis,gemstat(7,:), 'r');
plot(xaxis, LacZ_mean(line322,xx), 'b');
ylabel('tlm18')
set_fig;

% Line 323 (tlm19; spacer A)
subplot(5,1,1)
hold on;
plot(xaxis,gemstat(8,:), 'r');
plot(xaxis, LacZ_mean(line323,xx), 'b');
ylabel('tlm19')
set_fig;

% Line 324 (tlm20; spacer C)
subplot(5,1,3)
hold on;
plot(xaxis,gemstat(9,:), 'r');
plot(xaxis, LacZ_mean(line324,xx), 'b');
ylabel('tlm20')
set_fig;

% Line 240 (tlm11; 1000bp)
subplot(5,1,5)
hold on;
plot(xaxis,gemstat(6,:), 'r');
plot(xaxis, LacZ_mean(line240,xx), 'b');
ylabel('tlm11')
set_fig;


% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 1]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 1], 'YTick', [0 1]);

% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_single_spacer.eps'), 'epsc')


%% Plotting double GRF fits to spacers

% figure_params;
figure('Position', [32 400 415 1200]);
xSize = 8;  ySize = 20;  % Set width and height of figure in cm

% Line 356 (tlm07; spacer B)
subplot(5,1,2)
hold on;
plot(xaxis,twogrf(5,:), 'r');
plot(xaxis, LacZ_mean(line356,xx), 'b');
ylabel('tlm07')
set_fig;

% Line 322 (tlm18; spacer D)
subplot(5,1,4)
hold on;
plot(xaxis,twogrf(7,:), 'r');
plot(xaxis, LacZ_mean(line322,xx), 'b');
ylabel('tlm18')
set_fig;

% Line 323 (tlm19; spacer A)
subplot(5,1,1)
hold on;
plot(xaxis,twogrf(8,:), 'r');
plot(xaxis, LacZ_mean(line323,xx), 'b');
ylabel('tlm19')
set_fig;

% Line 324 (tlm20; spacer C)
subplot(5,1,3)
hold on;
plot(xaxis,twogrf(9,:), 'r');
plot(xaxis, LacZ_mean(line324,xx), 'b');
ylabel('tlm20')
set_fig;

% Line 240 (tlm11; 1000bp)
subplot(5,1,5)
hold on;
plot(xaxis,twogrf(6,:), 'r');
plot(xaxis, LacZ_mean(line240,xx), 'b');
ylabel('tlm11')
set_fig;


% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 1]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 1], 'YTick', [0 1]);

% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_twogrf_spacer.eps'), 'epsc')


%% Plotting GEMSTAT-GL window fits to spacers

norm_mat_gl = max(gemstat_gl,[],2);

%%
% figure_params;
figure('Position', [32 400 415 1200]);
xSize = 8;  ySize = 20;  % Set width and height of figure in cm

% Line 356 (tlm07; spacer B)
subplot(5,1,2)
hold on;
plot(xaxis, gemstat_gl(21,:)./norm_mat_gl(21), 'b');
plot(xaxis, gemstat_gl(22,:)./norm_mat_gl(21), 'r', xaxis, gemstat_gl(23,:)./norm_mat_gl(21), 'g');
ylabel('tlm07')
set_fig;

% Line 322 (tlm18; spacer D)
subplot(5,1,4)
hold on;
plot(xaxis, gemstat_gl(31,:)./norm_mat_gl(31), 'b');
plot(xaxis, gemstat_gl(32,:)./norm_mat_gl(31), 'r', xaxis, gemstat_gl(33,:)./norm_mat_gl(31), 'g');
ylabel('tlm18')
set_fig;

% Line 323 (tlm19; spacer A)
subplot(5,1,1)
hold on;
plot(xaxis, gemstat_gl(36,:)./norm_mat_gl(36), 'b');
plot(xaxis, gemstat_gl(37,:)./norm_mat_gl(36), 'r', xaxis, gemstat_gl(38,:)./norm_mat_gl(36), 'g');
ylabel('tlm19')
set_fig;

% Line 324 (tlm20; spacer C)
subplot(5,1,3)
hold on;
plot(xaxis, gemstat_gl(41,:)./norm_mat_gl(41), 'b');
plot(xaxis, gemstat_gl(42,:)./norm_mat_gl(42), 'r', xaxis, gemstat_gl(43,:)./norm_mat_gl(41), 'g');
ylabel('tlm20')
set_fig;

% Line 240 (tlm11; 1000bp)
subplot(5,1,5)
hold on;
plot(xaxis, gemstat_gl(26,:)./norm_mat_gl(26), 'b');
plot(xaxis, gemstat_gl(27,:), 'r', xaxis, gemstat_gl(28,:), 'g');
ylabel('tlm11')
set_fig;


% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 1]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 1], 'YTick', [0 1]);

% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_GL1_spacer.eps'), 'epsc')

%% Plotting GEMSTAT-GL fits to spacers

% figure_params;
figure('Position', [32 400 415 1200]);
xSize = 8;  ySize = 20;  % Set width and height of figure in cm


% Line 356 (tlm07; spacer B)
subplot(5,1,2)
hold on;
plot(xaxis, gemstat_gl(21,:)./norm_mat_gl(21), 'b');
plot(xaxis, gemstat_gl(25,:)./norm_mat_gl(25), 'r');
ylabel('tlm07')
set_fig;

% Line 322 (tlm18; spacer D)
subplot(5,1,4)
hold on;
plot(xaxis, gemstat_gl(31,:)./norm_mat_gl(31), 'b');
plot(xaxis, gemstat_gl(35,:)./norm_mat_gl(35), 'r');
ylabel('tlm18')
set_fig;

% Line 323 (tlm19; spacer A)
subplot(5,1,1)
hold on;
plot(xaxis, gemstat_gl(36,:)./norm_mat_gl(36), 'b');
plot(xaxis, gemstat_gl(40,:)./norm_mat_gl(40), 'r');
ylabel('tlm19')
set_fig;

% Line 324 (tlm20; spacer C)
subplot(5,1,3)
hold on;
plot(xaxis, gemstat_gl(41,:)./norm_mat_gl(41), 'b');
plot(xaxis, gemstat_gl(45,:)./norm_mat_gl(45), 'r');
ylabel('tlm20')
set_fig;

% Line 240 (tlm11; 1000bp)
subplot(5,1,5)
hold on;
plot(xaxis, gemstat_gl(26,:)./norm_mat_gl(26), 'b');
plot(xaxis, gemstat_gl(30,:)./norm_mat_gl(30), 'r');
ylabel('tlm11')
set_fig;

% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 1]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 1], 'YTick', [0 1]);

% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_GL2_spacer.eps'), 'epsc')


%% ----------------------------------------
% Plot tlm17


% Line 321 (tlm17; fusion E) windows
% figure_params;
figure('Position', [32 400 415 300]);
xSize = 8;  ySize = 5;  % Set width and height of figure in cm
hold on;
plot(xaxis, tlm17(1,:)./max(tlm17(1,:)), 'b');
plot(xaxis, tlm17(2,:)./max(tlm17(1,:)), 'r', xaxis, tlm17(3,:)./max(tlm17(1,:)), 'g');
ylabel('tlm17')
set_fig;

% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 1]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 0.5 1], 'YTick', [0 1]);

% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_GL1_tlm17.eps'), 'epsc')


% Line 321 (tlm17; fusion E) total model
% figure_params;
figure('Position', [32 400 415 300]);
xSize = 8;  ySize = 5;  % Set width and height of figure in cm
hold on;
plot(xaxis, tlm17(1,:)./max(tlm17(1,:)), 'b');
plot(xaxis, tlm17(5,:)/max(tlm17(5,:)), 'r');
ylabel('tlm17')
set_fig;


% Set general figure properties
set(findobj(gcf, 'Type', 'axes'), 'XLim', [0 1], 'YLim', [0 1]);
set(findobj(gcf, 'Type', 'axes'), 'XTick', [0 0.5 1], 'YTick', [0 1]);

% Save figure
saveas(gcf, fullfile(save_dir,'gemstat_GL2_tlm17.eps'), 'epsc')

