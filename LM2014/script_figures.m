% Wrapper script to make figures for paper using atlas registered hkb
% co-stained data.
%
% Tara Martin, 2014-02-04
% Edits:
% 2014-04-12: Added Figure 2 plots level plots.

%% Load and normalize atlas data (can be shared with other scripts)

script_load_and_normalize_atlas;

%% Initialize variables
ap_bins = 200;
line_mean = zeros(length(line_numbers),ap_bins);
line_std = line_mean;
line_sem = line_mean;

%% For each line extract line traces
for gg = 1:length(line_numbers)
% for gg = 1
    fprintf('Analyzing line %s\n', line_numbers{gg});
    nPces = size(atlas_data{gg},1);
    values_hkb = zeros(nPces,ap_bins);
    pos = zeros(nPces,4);
    
    % extract line trace for each embryo
    for e=1:nPces
        % convert atlas data to pointcloud
        pce = pointcloud(tm36.X{6}, atlas_data{gg}(e,:));
        pce = alignpce(pce);
    
        % extract line trace and normalize by hkb
        pce_trace = getapprojection(pce([0,1],dv_slice), ap_bins,0.005, false);
        values_hkb(e,:) = pce_trace/hkb_mean{gg}(e); % divide traces by mean hkb
    end
    
    %% extract mean and std for each line
    line_mean(gg,:) = mean(values_hkb);
    line_std(gg,:) = std(values_hkb);
    line_sem(gg,:) = std(values_hkb)/sqrt(nPces);
    
end

%% Normalize line traces by line 204
% line_mean_by204 = line_mean/tm36_mean204;
% line_std_by204 = line_std/tm36_mean204;
% line_sem_by204 = line_sem/tm36_mean204;

%% Now process control lines from hybridization tm-35

% Line numbers
line_numbers_controls = {'LacZ:0204','LacZ:0214','LacZ:0325','LacZ:0326',...
    'LacZ:0327','LacZ:0328','LacZ:0329','LacZ:0330','LacZ:0453',...
    'LacZ:0454','LacZ:0455','LacZ:0456'};

%% Initialize variables
controls_mean = zeros(length(line_numbers_controls),ap_bins);
controls_std = zeros(length(line_numbers_controls),ap_bins);
controls_sem = zeros(length(line_numbers_controls),ap_bins);

%% For each line extract line traces
for gg = 1:length(line_numbers_controls)
% for gg = 1
    fprintf('Analyzing line %s\n', line_numbers_controls{gg});
    nPces = size(atlas_data_controls{gg},1);
    values_hkb = zeros(nPces,ap_bins);
    
    % extract line trace for each embryo
    for e=1:nPces
        % convert atlas data to pointcloud
        pce = pointcloud(tm36.X{6}, atlas_data_controls{gg}(e,:));
        pce = alignpce(pce);
    
        % extract line trace and normalize by hkb
        pce_trace = getapprojection(pce([0,1],dv_slice), ap_bins,0.005, false);
        values_hkb(e,:) = pce_trace/hkb_mean_controls{gg}(e); % divide traces by mean hkb
    end
    
    %% extract mean and std for each line
    controls_mean(gg,:) = mean(values_hkb);
    controls_std(gg,:) = std(values_hkb);
    controls_sem(gg,:) = std(values_hkb)/sqrt(nPces);
end

%% Normalize controls by line 204
[controls_mean_by204 controls_sem_by204] = errorprop('multiply', ...
    controls_mean, controls_sem, cross_hybe, cross_hybe_sem);


%% Prep for making figures

close all;
myfiguredefaults;

%set some plot axis properties
xlimits = [30 90]*(ap_bins/100);
xticks = [30 60 90]*(ap_bins/100);
xlabels = {0.3, 0.6, 0.9};
ylimits = [-0.1 3.8];
yticks = [0 1 2 3];

red = [1 0.1 0.1];
bl = [0 0 0];
% cmap = [bl; red];
xx = 30*(ap_bins/100):90*(ap_bins/100);

%optional setting to make each spacer length/type have unique color
% colororder = cbrewer('qual', 'Set1', 5);
additive = [35 31 32]/255;
big_spacer = [190 30 45]/255;
small_spacer = [251 176 64]/255;
fusions = [43 57 144]/255;
inversions = [0 174 239]/255;
distal_fusion = [141 198 63]/255;


%% Figure 2: Comparing expression levels for single enhancer lines 
%            (position plots in script_thresholded_position_controls.m)

panel_size = [0 0 6 2.5]; %dimensions in cm
colors = cbrewer('qual', 'Set1', 5);

%holder for line expression values
expression99 = cell(size(hkb_mean_controls));
mu99 = zeros(size(line_numbers_controls));
sem99 = mu99;
max99 = mu99;

%get the 99% expression value for each construct
for gg = 1:length(line_numbers_controls)
    %first normalize each cell of each embryo
    temp = atlas_data_controls{gg}./repmat(hkb_mean_controls{gg}, [1 6078]);
    %get 99% value in trunk
    temp_trimmed = temp(:,I_cells); %extract only trunk cells
    expression99{gg} = prctile(temp_trimmed, 99, 2); %get 99% values
    
    %save mean and sem values for each line
    mu99(gg) = mean(expression99{gg});
    sem99(gg) = sem(expression99{gg});
    
    %in order to plot nicely, need to check what the max values of
    %expression99 are
    max99(gg) = max(expression99{gg});
end


distances = [-1000 -500 0];
fig2a_lines = [i453c i326c i204c]; %eve 3/7 forward
fig2b_lines = [i454c i327c i325c]; %eve 3/7 reverse
fig2d_lines = [i455c i329c i214c]; %eve 4/6 forward
fig2e_lines = [i456c i330c i328c]; %eve 4/6 reverse

figure(1);
clf;
hold on;
for i = 1:3
    jitter = (rand(size(expression99{fig2a_lines(i)}))-0.5)*20;
    %add a little random "jitter" to aid visibility
    plot(distances(i)+jitter, expression99{fig2a_lines(i)}, '.', 'Color',[.6 .6 .6], 'LineWidth', 1);
end
errorbar(distances, mu99(fig2a_lines), 1.96*sem99(fig2a_lines), '--', 'Color', colors(2,:), 'LineWidth', 1);
hold off;
ylim([0 max(max99)+0.02]); %make sure all data gets plotted
xlim([-1200 200])
set(gca, 'XTick', [-1000 -500 0],'xticklabel',{'','',''}, 'YTick', [0 1 2 3]);
set(gcf,'PaperPosition', panel_size); %figure dimensions in cm
set(gca,'FontSize', 8); % use 8pt font

% Save panel
saveas(gcf, fullfile(save_dir,'fig2A_levels.eps'), 'epsc')


figure(2);
clf;
hold on;
for i = 1:3
    jitter = (rand(size(expression99{fig2b_lines(i)}))-0.5)*20;
    %add a little random "jitter" to aid visibility
    plot(distances(i)+jitter, expression99{fig2b_lines(i)}, '.', 'Color',[.6 .6 .6], 'LineWidth', 1);
end
errorbar(distances, mu99(fig2b_lines), 1.96*sem99(fig2b_lines), '--', 'Color', colors(3,:), 'LineWidth', 1);
hold off;
ylim([0 max(max99)+0.02]); %make sure all data gets plotted
xlim([-1200 200])
set(gca, 'XTick', [-1000 -500 0],'xticklabel',{'','',''}, 'YTick', [0 1 2 3]);
set(gcf,'PaperPosition', panel_size); %figure dimensions in cm
set(gca,'FontSize', 8); % use 8pt font

% Save panel
saveas(gcf, fullfile(save_dir,'fig2B_levels.eps'), 'epsc')

figure(3);
clf;
hold on;
for i = 1:3
    jitter = (rand(size(expression99{fig2d_lines(i)}))-0.5)*20;
    %add a little random "jitter" to aid visibility
    plot(distances(i)+jitter, expression99{fig2d_lines(i)}, '.', 'Color',[.6 .6 .6], 'LineWidth', 1);
end
errorbar(distances, mu99(fig2d_lines), 1.96*sem99(fig2d_lines), '--', 'Color', colors(4,:), 'LineWidth', 1);
hold off;
ylim([0 max(max99)+0.02]); %make sure all data gets plotted
xlim([-1200 200])
set(gca, 'XTick', [-1000 -500 0],'xticklabel',{'','',''}, 'YTick', [0 1 2 3]);
set(gcf,'PaperPosition', panel_size); %figure dimensions in cm
set(gca,'FontSize', 8); % use 8pt font

% Save panel
saveas(gcf, fullfile(save_dir,'fig2D_levels.eps'), 'epsc')

figure(4);
clf;
hold on;
for i = 1:3
    jitter = (rand(size(expression99{fig2e_lines(i)}))-0.5)*20;
    %add a little random "jitter" to aid visibility
    plot(distances(i)+jitter, expression99{fig2e_lines(i)}, '.', 'Color',[.6 .6 .6], 'LineWidth', 1);
end
errorbar(distances, mu99(fig2e_lines), 1.96*sem99(fig2e_lines), '--', 'Color', colors(1,:), 'LineWidth', 1);
hold off;
ylim([0 max(max99)+0.02]); %make sure all data gets plotted
xlim([-1200 200])
set(gca, 'XTick', [-1000 -500 0],'xticklabel',{'','',''}, 'YTick', [0 1 2 3]);
set(gcf,'PaperPosition', panel_size); %figure dimensions in cm
set(gca,'FontSize', 8); % use 8pt font

% Save panel
saveas(gcf, fullfile(save_dir,'fig2E_levels.eps'), 'epsc')

%%
panel_size2 = [0 0 3.8 4]; %dimensions in cm

figure(5);
clf;
hold on;
errorbar(distances, mu99(fig2a_lines), 1.96*sem99(fig2a_lines), '.-', 'Color', colors(2,:), 'LineWidth', 0.5);
errorbar(distances, mu99(fig2b_lines), 1.96*sem99(fig2b_lines), '.-', 'Color', colors(3,:), 'LineWidth', 0.5);
hold off;
ylim([0 max(max99)+0.02]); %make sure all data gets plotted
xlim([-1200 100])
set(gca, 'XTick', [-1000 -500 0],'xticklabel',{'','',''}, 'YTick', [0 1 2 3]);
set(gcf,'PaperPosition', panel_size2); %figure dimensions in cm
set(gca,'FontSize', 8); % use 8pt font

% Save panel
saveas(gcf, fullfile(save_dir,'fig2C_levels.eps'), 'epsc')

figure(6);
clf;
hold on;
errorbar(distances, mu99(fig2d_lines), 1.96*sem99(fig2d_lines), '.-', 'Color', colors(4,:), 'LineWidth', 0.5);
errorbar(distances, mu99(fig2e_lines), 1.96*sem99(fig2e_lines), '.-', 'Color', colors(1,:), 'LineWidth', 0.5);
hold off;
ylim([0 max(max99)+0.02]); %make sure all data gets plotted
xlim([-1200 100])
set(gca, 'XTick', [-1000 -500 0],'xticklabel',{'','',''}, 'YTick', [0 1 2 3]);
set(gcf,'PaperPosition', panel_size2); %figure dimensions in cm
set(gca,'FontSize', 8); % use 8pt font

% Save panel
saveas(gcf, fullfile(save_dir,'fig2F_levels.eps'), 'epsc')


%% Figure 3: 1000bp spacers compared to additive single enhancers
% Edited to plot 1000bp spacer constructs compared to expected additive
% model

fh = figure('Position', [32 400 415 1200]);
set(fh, 'DefaultAxesXLim', xlimits);
set(fh, 'DefaultAxesXTick', xticks);
set(fh, 'DefaultAxesXTickLabel', xlabels)
set(fh, 'DefaultAxesYLim', ylimits);
set(fh, 'DefaultAxesYTick', yticks);

xSize = 4.5;  ySize = 13; xLeft = 0.1; yBottom = 0.1; % Set width and height of figure in cm
set(fh,'PaperPosition', [xLeft yBottom xSize ySize])

set(fh, 'name', 'Figure 3, 1000bp spacers')

% big_spacer = red;
% additive = bl;

subplot(412) % equivalent to line 207
hold on;
plotlinetrace(xx, line_mean(i240,:), line_sem(i240,:), big_spacer);
[c_mean c_sem] = errorprop('add',...
    controls_mean_by204(i453c,:), controls_sem_by204(i453c,:),...
    line_mean(i214,:), line_sem(i214,:) );
plotlinetrace(xx, c_mean, c_sem, additive);
% title('37+1000bp+46');
ylabel({'Normalized', 'expression'})

subplot(414) % equivalent to line 208
hold on;
plotlinetrace(xx, line_mean(i450,:), line_sem(i450,:), big_spacer);
[c_mean c_sem] = errorprop('add',...
    controls_mean_by204(i455c,:), controls_sem_by204(i455c,:),...
    line_mean(i204,:), line_sem(i204,:) );
plotlinetrace(xx, c_mean, c_sem, additive);
% title('46+1000bp+37');
ylabel({'Normalized', 'expression'})

subplot(411) % equivalent to line 209
hold on;
plotlinetrace(xx, line_mean(i451,:), line_sem(i451,:), big_spacer);
[c_mean c_sem] = errorprop('add',...
    controls_mean_by204(i454c,:), controls_sem_by204(i454c,:),...
    line_mean(i214,:), line_sem(i214,:) );
plotlinetrace(xx, c_mean, c_sem, additive);
% title('37R+1000bp+46');
ylabel({'Normalized', 'expression'})

subplot(413) % equivalent to line 210
hold on;
plotlinetrace(xx, line_mean(i452,:), line_sem(i452,:), big_spacer);
[c_mean c_sem] = errorprop('add',...
    controls_mean_by204(i453c,:), controls_sem_by204(i453c,:),...
    controls_mean_by204(i328c,:), controls_sem_by204(i328c,:) );
plotlinetrace(xx, c_mean, c_sem, additive);
% title('37+1000bp+46R');
ylabel({'Normalized', 'expression'})

% Save figure
saveas(gcf, fullfile(save_dir,'Figure3.eps'), 'epsc')


%% Figure 4: Shorter spacer allows greater interference between enhancers 
% (200bp spacers compared to 1000bp)

% big_spacer = bl;
% small_spacer = red;

fh = figure('Position', [32 400 415 1200]);
set(fh, 'DefaultAxesXLim', xlimits);
set(fh, 'DefaultAxesXTick', xticks);
set(fh, 'DefaultAxesXTickLabel', xlabels)
set(fh, 'DefaultAxesYLim', ylimits);
set(fh, 'DefaultAxesYTick', yticks);

xSize = 4.5;  ySize = 13.5; xLeft = 0.1; yBottom = 0.1; % Set width and height of figure in cm
set(fh,'PaperPosition', [xLeft yBottom xSize ySize])
set(fh, 'name', 'Figure 4, 1000bp vs 200bp')

subplot(412) % equivalent to line 207
hold on;
plotlinetrace(xx, line_mean(i356,:), line_sem(i356,:), small_spacer)
[c_mean c_sem] = errorprop('add',...
    controls_mean_by204(i453c,:), controls_sem_by204(i453c,:),...
    line_mean(i214,:), line_sem(i214,:) );
plotlinetrace(xx, c_mean, c_sem, additive);
plotlinetrace(xx, line_mean(i240,:), line_sem(i240,:), big_spacer)
% title('37+200bp+46');
ylabel({'Normalized', 'expression'})

subplot(414) % equivalent to line 208
hold on;
plotlinetrace(xx, line_mean(i322,:), line_sem(i322,:), small_spacer)
[c_mean c_sem] = errorprop('add',...
    controls_mean_by204(i455c,:), controls_sem_by204(i455c,:),...
    line_mean(i204,:), line_sem(i204,:) );
plotlinetrace(xx, c_mean, c_sem, additive);
plotlinetrace(xx, line_mean(i450,:), line_sem(i450,:), big_spacer)
% title('46+200bp+37');
ylabel({'Normalized', 'expression'})

subplot(411) % equivalent to line 209
hold on;
plotlinetrace(xx, line_mean(i323,:), line_sem(i323,:), small_spacer)
[c_mean c_sem] = errorprop('add',...
    controls_mean_by204(i454c,:), controls_sem_by204(i454c,:),...
    line_mean(i214,:), line_sem(i214,:) );
plotlinetrace(xx, c_mean, c_sem, additive);
plotlinetrace(xx, line_mean(i451,:), line_sem(i451,:), big_spacer)
% title('37R+200bp+46');
ylabel({'Normalized', 'expression'})

subplot(413) % equivalent to line 210
hold on;
plotlinetrace(xx, line_mean(i324,:), line_sem(i324,:), small_spacer)
[c_mean c_sem] = errorprop('add',...
    controls_mean_by204(i453c,:), controls_sem_by204(i453c,:),...
    controls_mean_by204(i328c,:), controls_sem_by204(i328c,:) );
plotlinetrace(xx, c_mean, c_sem, additive);
plotlinetrace(xx, line_mean(i452,:), line_sem(i452,:), big_spacer)
% title('37+200bp+46R');
ylabel({'Normalized', 'expression'})

% Save figure
saveas(gcf, fullfile(save_dir,'Figure4.eps'), 'epsc')


%% Figure 5: Local interactions explain changes in pattern but not level 
% (200bp spacers vs fusions)

xSize = 5.5;  ySize = 13; xLeft = 0.1; yBottom = 0.1; % Set width and height of figure in cm

fh = figure('Position', [600 400 415 1200]);
fig_intialize(fh, xlimits, xticks, ylimits, yticks, [xLeft yBottom xSize ySize]);
set(fh, 'name', 'Figure 5, 200bp spacers vs fusions')

% small_spacer = bl;
% fusions = red;

subplot(412)
hold on;
plotlinetrace(xx, line_mean(i207,:), line_sem(i207,:), fusions)
plotlinetrace(xx, line_mean(i356,:), line_sem(i356,:), small_spacer)
title('37+46');
ylabel({'Normalized', 'expression'})

subplot(414)
hold on;
plotlinetrace(xx, line_mean(i208,:), line_sem(i208,:), fusions)
plotlinetrace(xx, line_mean(i322,:), line_sem(i322,:), small_spacer)
title('46+37');
ylabel({'Normalized', 'expression'})

subplot(411)
hold on;
plotlinetrace(xx, line_mean(i209,:), line_sem(i209,:), fusions)
plotlinetrace(xx, line_mean(i323,:), line_sem(i323,:), small_spacer)
title('37R+46');
ylabel({'Normalized', 'expression'})

subplot(413)
hold on;
plotlinetrace(xx, line_mean(i210,:), line_sem(i210,:), fusions)
plotlinetrace(xx, line_mean(i324,:), line_sem(i324,:), small_spacer)
title('37+46R');
ylabel({'Normalized', 'expression'})

% Save figure
saveas(gcf, fullfile(save_dir,'Figure5.eps'), 'epsc')


%% Figure 6: Reversing fusions
xSize = 8;  ySize = 10; xLeft = 0.1; yBottom = 0.1; % Set width and height of figure in cm

fh = figure('Position', [32 400 600 600]);
fig_intialize(fh, xlimits, xticks, ylimits, yticks, [xLeft yBottom xSize ySize]);
set(fh, 'name', 'Figure 6, inverted fusions')

% fusions = bl;
% inversions = red;

% ylimits = [-0.1 3.8];
% yticks = [0 1 2 3];

subplot(211)
hold on;
title('207 vs 320')
plotlinetrace(xx, line_mean(i207,:), line_sem(i207,:), fusions)
plotlinetrace(xx, line_mean(i320,:), line_sem(i320,:), inversions)
ylabel({'Normalized', 'expression'})

subplot(212)
hold on;
title('208 vs 321');
plotlinetrace(xx, line_mean(i208,:), line_sem(i208,:), fusions)
plotlinetrace(xx, line_mean(i321,:), line_sem(i321,:), inversions)
ylabel({'Normalized', 'expression'})


% Save figure
saveas(gcf, fullfile(save_dir,'Figure6.eps'), 'epsc')


%% Figure 7:  1000bp spacer between promoter and fusions
ylimits = [-0.1 3.8];
yticks = [0 1 2 3];
xSize = 5.5;  ySize = 13; xLeft = 0.1; yBottom = 0.1; % Set width and height of figure in cm

% fusions = bl;
% distal_fusion = red;

fh = figure('Position', [32 400 415 1200]);
fig_intialize(fh, xlimits, xticks, ylimits, yticks, [xLeft yBottom xSize ySize]);
set(fh, 'name', 'Fusions at a distance from promoter (red)');

subplot(412)
hold on;
plotlinetrace(xx, line_mean(i207,:), line_sem(i207,:), fusions)
plotlinetrace(xx, line_mean(i458,:), line_sem(i458,:), distal_fusion)
title('37+46');
ylabel({'Normalized', 'expression'})

subplot(414)
hold on;
plotlinetrace(xx, line_mean(i208,:), line_sem(i208,:), fusions)
plotlinetrace(xx, line_mean(i460,:), line_sem(i460,:), distal_fusion)
title('46+37');
ylabel({'Normalized', 'expression'})

subplot(411)
hold on;
plotlinetrace(xx, line_mean(i209,:), line_sem(i209,:), fusions)
plotlinetrace(xx, line_mean(i457,:), line_sem(i457,:), distal_fusion)
title('37R+46');
ylabel({'Normalized', 'expression'})

subplot(413)
hold on;
plotlinetrace(xx, line_mean(i210,:), line_sem(i210,:), fusions)
plotlinetrace(xx, line_mean(i459,:), line_sem(i459,:), distal_fusion)
title('37+46R');
ylabel({'Normalized', 'expression'})

% Save figure
saveas(gcf, fullfile(save_dir,'Figure7.eps'), 'epsc')

%% Double check how 204 and 214 compare

fh = figure('Position', [32 400 415 1200]);
fig_intialize(fh, xlimits, xticks, ylimits, yticks, [xLeft yBottom xSize ySize]);
set(fh, 'name', 'testing hkb normalization');


subplot(211)
hold on;
plotlinetrace(xx, line_mean(1,:), line_sem(1,:), bl);
plotlinetrace(xx, controls_mean_by204(i204c,:), controls_sem_by204(i204c,:), red);

subplot(212)
hold on;
plotlinetrace(xx, line_mean(6,:), line_sem(6,:), bl);
plotlinetrace(xx, line_mean(i214,:), line_sem(i214,:), red);

%% Supplementary figure
% (Stripe 7 in 200bp spacers vs fusions)

%initialize variables
ap_bins = 200;

peak_pos = cell(length(line_numbers),1);
peak_pos_mean = zeros(length(line_numbers),1);
peak_pos_sem = zeros(length(line_numbers),1);
fusion_lines = [i204 i323 i209 i356 i207 i324 i210];

%% For each line extract line traces and stripe 7 peak position
for gg = fusion_lines
% for gg = 1
    fprintf('Analyzing line %s\n', line_numbers{gg});
    nPces = size(atlas_data{gg},1);
    pos = zeros(nPces,1);
    
    % extract line trace for each embryo
    for e=1:nPces
        % convert atlas data to pointcloud
        pce = pointcloud(tm36.X{6}, atlas_data{gg}(e,:));
        pce = alignpce(pce);
    
        % extract line trace
        pce_trace = getapprojection(pce([0,1],dv_slice), ap_bins,0.005, false);
        
        % save stripe peak info
        smooth_trace = double(gaussf(pce_trace,0.5));
        pos(e) = locateappeaks(smooth_trace, [0.82]);        
    end
    
    %% extract mean stripe 7 position and SEM for each line
    peak_pos{gg} = pos;
    peak_pos_mean(gg) = mean(pos);
    peak_pos_sem(gg) = 1.96*std(pos)/sqrt(nPces);

end

%% Now process control lines from hybridization tm-35

% Line numbers
line_numbers_controls = {'LacZ:0204','LacZ:0214','LacZ:0325','LacZ:0326',...
    'LacZ:0327','LacZ:0328','LacZ:0329','LacZ:0330','LacZ:0453',...
    'LacZ:0454','LacZ:0455','LacZ:0456'};

%% Initialize variables
peak_pos_controls = cell(length(line_numbers_controls),1);
peak_pos_mean_controls = zeros(length(line_numbers_controls),1);
peak_pos_sem_controls = zeros(length(line_numbers_controls),1);
control_lines = [i204c i325c i326c i327c i453c i454c];

%% For each line extract line traces
for gg = control_lines
% for gg = 1
    fprintf('Analyzing line %s\n', line_numbers_controls{gg});
    nPces = size(atlas_data_controls{gg},1);
    pos = zeros(nPces,1);

    % extract line trace for each embryo
    for e=1:nPces
        % convert atlas data to pointcloud
        pce = pointcloud(tm36.X{6}, atlas_data_controls{gg}(e,:));
        pce = alignpce(pce);
    
        % extract line trace
        pce_trace = getapprojection(pce([0,1],dv_slice), ap_bins,0.005, false);
        
        % save stripe peak info
        smooth_trace = double(gaussf(pce_trace,0.5));
        pos(e) = locateappeaks(smooth_trace, [0.82]);        
    end
    
    %% extract mean stripe 7 position and SEM for each line
    peak_pos_controls{gg} = pos;
    peak_pos_mean_controls(gg) = mean(pos);
    peak_pos_sem_controls(gg) = 1.96*std(pos)/sqrt(nPces);
end

%% Aggregate line 204 data

peak_204 = [peak_pos{i204}; peak_pos_controls{i204c}];
peak_204_mean = mean(peak_204);
peak_204_sem = 1.96*std(peak_204)/sqrt(length(peak_204));

%% T-test

% %holder variables
% h_controls = zeros(size(control_lines)); %reject null hypothesis (1=yes, 2=no)
% p_controls = h_controls;  %p-value
% 
% for gg = 1:length(control_lines)
%     [h_controls(gg) p_controls(gg)] = ttest2(peak_204, peak_pos_controls{control_lines(gg)});    
% end
% 
% %holder variables
% h = zeros(size(fusion_lines)); %reject null hypothesis (1=yes, 2=no)
% p = h;  %p-value
% 
% for gg = 1:length(fusion_lines)
%     [h(gg) p(gg)] = ttest2(peak_204, peak_pos{fusion_lines(gg)});    
% end
% 
% diff = peak_204_mean - peak_pos_mean(fusion_lines);

diff_A = peak_pos_mean(i323) - peak_pos_mean(i209)
[hA pA] = ttest2(peak_pos{i323}, peak_pos{i209})

diff_B = peak_pos_mean(i356) - peak_pos_mean(i207)
[hB pB] = ttest2(peak_pos{i356}, peak_pos{i207})

diff_C = peak_pos_mean(i324) - peak_pos_mean(i210)
[hC pC] = ttest2(peak_pos{i324}, peak_pos{i210})

%% Plot control lines
panel_size = [0.1 0.1 7.5 10]; %[xLeft yBottom xSize ySize] in cm
close all;
fh = figure('Position', [100 400 415 400]);
fig_intialize(fh, [0.82 0.86], [0.82 0.84 0.86], [0 6.5], [1 2 3 4 5 6], panel_size);
set(fh, 'name', 'Supplement, stripe 7 in 200bp spacers vs fusions')

% %first plot 204
% subplot(211)
% hold on;
% [xb, yb] = horizontal_errorbar(peak_204_mean,6,peak_204_sem, 0.1);
% plot(xb,yb)
% plot([peak_204_mean peak_204_mean],[6 0], '.-')
% 
% %then the rest of controls
% [xb, yb] = horizontal_errorbar(peak_pos_mean_controls(control_lines(2:end)),...
%     [5 4 3 2 1],peak_pos_sem_controls(control_lines(2:end)), 0.1);
% plot(xb,yb)
% set(gca, 'yticklabel', {'reverse+1000bp','forward+1000bp','reverse+500bp',... 
%     'forward+500bp', 'reverse','forward'})
% hold off;

%plot fusion lines
% subplot(212)
hold on;
% [xb, yb] = horizontal_errorbar(peak_204_mean,7,peak_204_sem, 0.1);
% plot(xb,yb)
% plot([peak_204_mean peak_204_mean],[7 0], '.-')

[xb, yb] = horizontal_errorbar(peak_pos_mean(fusion_lines(2:end)),...
    [6 5 4 3 2 1],peak_pos_sem(fusion_lines(2:end)), 0.1);
plot(xb,yb)
ylim([0 6.5])
set(gca, 'ytick', [1:7],'yticklabel', {'C_fusion', 'C_200',... 
    'B_fusion', 'B_200', 'A_fusion', 'A_200'})

plot([0.85 0.85], [6 5], 'k')
text(0.851, 5.5, ['p = ' num2str(round(pA*1000)/1000)])

plot([0.85 0.85], [4 3], 'k')
text(0.851, 3.5, ['p = ' num2str(round(pB*1000)/1000)])

plot([0.85 0.85], [2 1], 'k')
text(0.851, 1.5, ['p = ' num2str(round(pC*1000)/1000)])

hold off;
xlabel('x/L')

% Save figure
saveas(gcf, fullfile(save_dir,'SuppFig_stripe7.eps'), 'epsc')


%%
panel_size = [0.1 0.1 5.5 10]; %[xLeft yBottom xSize ySize] in cm

close all;
fh = figure('Position', [100 400 415 1200]);
fig_intialize(fh, [0.8 0.9], [0.8 0.85 0.9], [0 2.5], [1 2], panel_size);
set(fh, 'name', 'Supplement, stripe 7 in 200bp spacers vs fusions')

% small_spacer = bl;
% fusions = red;

subplot(311)
[xb, yb] = horizontal_errorbar(peak_pos([i323 i209],4),[2 1],peak_pos_sem([i323 i209],4), 0.1);
hold on;
plot(xb,yb)
plot([peak_pos(i323,4) peak_pos(i323,4)],[2 0], '.-')
plot([peak_pos(i209,4) peak_pos(i209,4)],[1 0], 'g.-')
set(gca, 'yticklabel',{'A_200', 'A_fusion'})
text(peak_pos(i323,4)+0.002, 0.5, ['diff = ' num2str(peak_pos(i323,4)-peak_pos(i209,4))])


subplot(312)
[xb, yb] = horizontal_errorbar(peak_pos([i356 i207],4),[2 1],peak_pos_sem([i356 i207],4), 0.1);
hold on;
plot(xb,yb)
% plot(peak_pos([i356 i207],4),[2 1], '.')
plot([peak_pos(i356,4) peak_pos(i356,4)],[2 0], '.-')
plot([peak_pos(i207,4) peak_pos(i207,4)],[1 0], 'g.-')
set(gca, 'yticklabel',{'B_200', 'B_fusion'})
text(peak_pos(i356,4)+0.002, 0.5, ['diff = ' num2str(peak_pos(i356,4)-peak_pos(i207,4))])

subplot(313)
[xb, yb] = horizontal_errorbar(peak_pos([i324 i210],4),[2 1],peak_pos_sem([i324 i210],4), 0.1);
hold on;
plot(xb,yb)
% plot(peak_pos([i324 i210],4),[2 1], '.')
plot([peak_pos(i324,4) peak_pos(i324,4)],[2 0], '.-')
plot([peak_pos(i210,4) peak_pos(i210,4)],[1 0], 'g.-')
set(gca, 'yticklabel',{'C_200', 'C_fusion'})
text(peak_pos(i324,4)+0.002, 0.5, ['diff = ' num2str(peak_pos(i324,4)-peak_pos(i210,4))])

% subplot(414)
% hold on;
% plotlinetrace(xx, line_mean(i208,:), line_sem(i208,:), fusions)
% plotlinetrace(xx, line_mean(i322,:), line_sem(i322,:), small_spacer)
% title('46+37');
% ylabel({'Normalized', 'expression'})
% 
% subplot(411)
% hold on;
% plotlinetrace(xx, line_mean(i209,:), line_sem(i209,:), fusions)
% plotlinetrace(xx, line_mean(i323,:), line_sem(i323,:), small_spacer)
% title('37R+46');
% ylabel({'Normalized', 'expression'})
% 
% subplot(413)
% hold on;
% plotlinetrace(xx, line_mean(i210,:), line_sem(i210,:), fusions)
% plotlinetrace(xx, line_mean(i324,:), line_sem(i324,:), small_spacer)
% title('37+46R');
% ylabel({'Normalized', 'expression'})

% Save figure
saveas(gcf, fullfile(save_dir,'SuppFig_stripe7.eps'), 'epsc')



