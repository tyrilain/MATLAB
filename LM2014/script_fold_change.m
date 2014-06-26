% Extracting mean lacZ expression for each stripe using thresholding of
% stripes in single enhancer controls to define "on" cells.
%
% Tara Martin, 2014-02-09
% Edits:
% 2014-03-07: Made loading and normalizing data into separate script.

%% Load and normalize atlas data (can be shared with other scripts)

% script_load_and_normalize_atlas;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: Define stripe cells using thresholding of lines 204 and 214

% Threshold controls to find stripe cells

% in tm-36 the indices of line 204=1 and line 214=6
% in tm-35 the indices of line 204=1 and line 214=2

% pull out 204 data
data_204 = [atlas_data{1}; atlas_data_controls{1}];
nEmbryos(1) = size(data_204,1);
% threshold each embryo individually
for i = 1:nEmbryos(1)
    embryo = data_204(i,I_cells);   %trunk of embryo
    t = findThreshold(embryo);      %find "on" threshold
    data_204(i,:) = data_204(i,:)>t;  %save on/off values of cells
end

% pull out 214 data
data_214 = [atlas_data{6}; atlas_data_controls{2}];
nEmbryos(2) = size(data_214,1);
% threshold each embryo individually
for i = 1:nEmbryos(2)
    embryo = data_214(i,I_cells);   %trunk of embryo
    t = findThreshold(embryo);      %find "on" threshold
    data_214(i,:) = data_214(i,:)>t;  %save on/off values of cells
end

%% Visualize thresholded stripes
% subplot(2,1,1)
% pce = pointcloud(coords, mean(data_204));
% disp(pce, 'print');
% subplot(2,1,2)
% pce = pointcloud(coords, mean(data_214));
% disp(pce, 'print');

%% now find cells in each stripe (indexed by I_3, I_4, I_6, I_7)

thresh = 0.3;   % threshold on the probability of stripe being expressed

I_204 = find(mean(data_204)>thresh); %indices of "on" cells
I_3 = find(coords(:,1)>0.40 & coords(:,1)<0.6); %limit search to area of embryo
I_3 = intersect(I_204', I_3);   %cells in right region that are "on"
I_7 = find(coords(:,1)>0.60 & coords(:,1)<0.90); %limit search to area of embryo
I_7 = intersect(I_204', I_7);   %cells in right region that are "on"


I_214 = find(mean(data_214)>thresh);
I_4 = find(coords(:,1)>0.45 & coords(:,1)<0.65);
I_4 = intersect(I_214', I_4);
I_6 = find(coords(:,1)>0.65 & coords(:,1)<0.85);
I_6 = intersect(I_214', I_6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2: Calculate mean and 95% of each stripe for each line

% intialize holder variables
stripe_means = zeros(length(line_numbers), 4);
stripe_std = stripe_means;
stripe_sem = stripe_means;
% stripe_95 = stripe_means;
% stripe_95_std = stripe_means;
% stripe_95_sem = stripe_means;

for gg = 1:length(line_numbers)
    fprintf('Analyzing line %s\n', line_numbers{gg});
    n = length(hkb_mean{gg}); % number of embryos for line
    
    % initialize temporary holder variables
    str_mean = zeros(n,4);
    str_95 = str_mean;
    
    % normalize and find mean and 95% for each embryo
    for e = 1:n
        embryo = atlas_data{gg}(e,:)/hkb_mean{gg}(e); %normalized expression
        
        % save mean and 95% lacZ expression for each stripe in each embryo
        str_mean(e,:) = [mean(embryo(I_3)) mean(embryo(I_7))...
                         mean(embryo(I_4)) mean(embryo(I_6))];
        str_95(e,:) = [prctile(embryo(I_3),95) prctile(embryo(I_7),95)...
                       prctile(embryo(I_4),95) prctile(embryo(I_6),95)];
    end
    
    % calculate stripe mean, std, and sem of stripe means and 95%
    stripe_means(gg,:) = mean(str_mean);
    stripe_std(gg,:) = std(str_mean);
    stripe_sem(gg,:) = std(str_mean)/sqrt(n);
    
%     stripe_95(gg,:) = mean(str_95);
%     stripe_95_std(gg,:) = std(str_95);
%     stripe_95_sem(gg,:) = std(str_95)/sqrt(n);
    
end

%% Do it again for the control lines

% intialize holder variables
control_means = zeros(length(line_numbers_controls), 4);
control_std = control_means;
control_sem = control_means;
stripe_95_controls = control_means;
stripe_95_std_controls = control_means;
stripe_95_sem_controls = control_means;

for gg = 1:length(line_numbers_controls)
    fprintf('Analyzing line %s\n', line_numbers_controls{gg});
    n = length(hkb_mean_controls{gg}); % number of embryos for line
    
    % initialize temporary holder variables
    str_mean = zeros(n,4);
    str_95 = str_mean;
    
    % normalize and find mean and 95% for each embryo
    for e = 1:n
        embryo = atlas_data_controls{gg}(e,:)/hkb_mean_controls{gg}(e); %normalized expression
        str_mean(e,:) = [mean(embryo(I_3)) mean(embryo(I_7))...
                         mean(embryo(I_4)) mean(embryo(I_6))];
        str_95(e,:) = [prctile(embryo(I_3),95) prctile(embryo(I_7),95)...
                       prctile(embryo(I_4),95) prctile(embryo(I_6),95)];
    end
    
    % calculate mean, std, and sem of stripe means and 95%
    control_means(gg,:) = mean(str_mean);
    control_std(gg,:) = std(str_mean);
    control_sem(gg,:) = std(str_mean)/sqrt(n);
    
    stripe_95_controls(gg,:) = mean(str_95);
    stripe_95_std_controls(gg,:) = std(str_95);
    stripe_95_sem_controls(gg,:) = std(str_95)/sqrt(n);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2: Normalize all lines by 204

% divide stripe means by line 204 mean
[stripe_means stripe_sem] = errorprop('divide',stripe_means, stripe_sem,... 
    repmat(tm36_mean204, size(stripe_means)), repmat(tm36_sem204, size(stripe_means)));

% also for single enhancer lines
[control_means control_sem] = errorprop('divide',control_means, control_sem,... 
    repmat(tm35_mean204, size(control_means)), repmat(tm35_sem204, size(control_means)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3: Calculate the fold change for each line (have to do this one at a time)

fold_change = ones(length(line_numbers), 4);
fold_error = ones(length(line_numbers), 4);
%[fold_change fold_error] = stripe_fold_change(stripe_mean, stripe_err, control_mean, control_err)

%line 240
fprintf('Analyzing line 240.\n');
[fold_change(i240,1:2) fold_error(i240, 1:2)] = stripe_fold_change(...
    stripe_means(i240,1:2),stripe_sem(i240,1:2),... 
    control_means(i453c, 1:2), control_sem(i453c, 1:2));

[fold_change(i240,3:4) fold_error(i240, 3:4)] = stripe_fold_change(...
    stripe_means(i240,3:4),stripe_sem(i240,3:4),... 
    stripe_means(i214, 3:4), stripe_sem(i214, 3:4));

%line 450
% fold_change(i450,:) = stripe_fold_change(i450,[i204c i455c], stripe_means, control_means);
fprintf('Analyzing line 450.\n');
[fold_change(i450,1:2) fold_error(i450, 1:2)] = stripe_fold_change(...
    stripe_means(i450,1:2),stripe_sem(i450,1:2),... 
    stripe_means(i204, 1:2), stripe_sem(i204, 1:2));

[fold_change(i450,3:4) fold_error(i450, 3:4)] = stripe_fold_change(...
    stripe_means(i450,3:4),stripe_sem(i450,3:4),... 
    control_means(i455c, 3:4), control_sem(i455c, 3:4));

%line 451
% fold_change(i451,:) = stripe_fold_change(i451,[i454c i214c], stripe_means, control_means);
fprintf('Analyzing line 451.\n');
[fold_change(i451,1:2) fold_error(i451, 1:2)] = stripe_fold_change(...
    stripe_means(i451,1:2),stripe_sem(i451,1:2),... 
    control_means(i454c, 1:2), control_sem(i454c, 1:2));

[fold_change(i451,3:4) fold_error(i451, 3:4)] = stripe_fold_change(...
    stripe_means(i451,3:4),stripe_sem(i451,3:4),... 
    stripe_means(i214, 3:4), stripe_sem(i214, 3:4));

%line 452
% fold_change(i452,:) = stripe_fold_change(i452,[i453c i328c], stripe_means, control_means);
fprintf('Analyzing line 452.\n');
[fold_change(i452,1:2) fold_error(i452, 1:2)] = stripe_fold_change(...
    stripe_means(i452,1:2),stripe_sem(i452,1:2),... 
    control_means(i453c, 1:2), control_sem(i453c, 1:2));

[fold_change(i452,3:4) fold_error(i452, 3:4)] = stripe_fold_change(...
    stripe_means(i452,3:4),stripe_sem(i452,3:4),... 
    control_means(i328c, 3:4), control_sem(i328c, 3:4));

%%
%line 356
% fold_change(i356,:) = stripe_fold_change(i356,[i453c i214c], stripe_means, control_means);
fprintf('Analyzing line 356.\n');
[fold_change(i356,1:2) fold_error(i356, 1:2)] = stripe_fold_change(...
    stripe_means(i356,1:2),stripe_sem(i356,1:2),... 
    control_means(i453c, 1:2), control_sem(i453c, 1:2));

[fold_change(i356,3:4) fold_error(i356, 3:4)] = stripe_fold_change(...
    stripe_means(i356,3:4),stripe_sem(i356,3:4),... 
    stripe_means(i214, 3:4), stripe_sem(i214, 3:4));

%line 322
% fold_change(i322,:) = stripe_fold_change(i322,[i204c i329c], stripe_means, control_means);
fprintf('Analyzing line 322.\n');
[fold_change(i322,1:2) fold_error(i322, 1:2)] = stripe_fold_change(...
    stripe_means(i322,1:2),stripe_sem(i322,1:2),... 
    stripe_means(i204, 1:2), stripe_sem(i204, 1:2));

[fold_change(i322,3:4) fold_error(i322, 3:4)] = stripe_fold_change(...
    stripe_means(i322,3:4),stripe_sem(i322,3:4),... 
    control_means(i329c, 3:4), control_sem(i329c, 3:4));

%line 323
% fold_change(i323,:) = stripe_fold_change(i323,[i454c i214c], stripe_means, control_means);
fprintf('Analyzing line 323.\n');
[fold_change(i323,1:2) fold_error(i323, 1:2)] = stripe_fold_change(...
    stripe_means(i323,1:2),stripe_sem(i323,1:2),... 
    control_means(i454c, 1:2), control_sem(i454c, 1:2));

[fold_change(i323,3:4) fold_error(i323, 3:4)] = stripe_fold_change(...
    stripe_means(i323,3:4),stripe_sem(i323,3:4),... 
    stripe_means(i214, 3:4), stripe_sem(i214, 3:4));

%line 324
% fold_change(i324,:) = stripe_fold_change(i324,[i453c i328c], stripe_means, control_means);
fprintf('Analyzing line 324.\n');
[fold_change(i324,1:2) fold_error(i324, 1:2)] = stripe_fold_change(...
    stripe_means(i324,1:2),stripe_sem(i324,1:2),... 
    control_means(i453c, 1:2), control_sem(i453c, 1:2));

[fold_change(i324,3:4) fold_error(i324, 3:4)] = stripe_fold_change(...
    stripe_means(i324,3:4),stripe_sem(i324,3:4),... 
    control_means(i328c, 3:4), control_sem(i328c, 3:4));

%% Adjust fold_error to plot 95% confidence interval
fold_error = 1.96*fold_error;


%%
% figure(1);
% bar(log2(fold_change([i240 i450 i451 i452], :)))
% ylim([-2 2])
% 
% 
% figure(2);
% bar(log2(fold_change([i356 i322 i323 i324], :)))
% ylim([-2 2])

close all;
myfiguredefaults;
xlimits = [0 6];
xticks = [1 2 4 5];
ylimits = [-2.8 2.8];
yticks = [-2 0 2];
xSize = 3;  ySize = 13; xLeft = 0.1; yBottom = 0.1; % Set width and height of figure in cm



%% Figure 3: 1000bp spacers

fh = figure(1);
fig_intialize(fh, xlimits, xticks, ylimits, yticks, [xLeft yBottom xSize ySize]);


subplot(412) % equivalent to line 207
barwitherr(fold_error(i240,:), xticks, fold_change(i240, :))
set(gca, 'XTickLabel',{'3','7','4','6' });
set(gca, 'Box', 'off')
ylabel('log2(fold-change)')

subplot(414) % equivalent to line 208
barwitherr(fold_error(i450,:),xticks, fold_change(i450, :))
set(gca, 'XTickLabel',{'3','7','4','6' });
set(gca, 'Box', 'off')
ylabel('log2(fold-change)')

subplot(411) % equivalent to line 209
barwitherr(fold_error(i451,:),xticks, fold_change(i451, :))
set(gca, 'XTickLabel',{'3','7','4','6' });
set(gca, 'Box', 'off')
ylabel('log2(fold-change)')

subplot(413) % equivalent to line 210
barwitherr(fold_error(i452,:),xticks, fold_change(i452, :))
set(gca, 'XTickLabel',{'3','7','4','6' });
set(gca, 'Box', 'off')
ylabel('log2(fold-change)')

colormap([0.7 0.7 0.7])

% Save figure
saveas(gcf, fullfile(save_dir,'Figure3_fold_change.eps'), 'epsc')

%% Figure 4: 200bp spacers

fh = figure(2);
fig_intialize(fh, xlimits, xticks, ylimits, yticks, [xLeft yBottom xSize ySize]);


subplot(412) % equivalent to line 207
barwitherr(fold_error(i356,:),xticks, fold_change(i356, :))
set(gca, 'XTickLabel',{'3','7','4','6' });
set(gca, 'Box', 'off')
ylabel('log2(fold-change)')

subplot(414) % equivalent to line 208
barwitherr(fold_error(i322,:),xticks, fold_change(i322, :))
set(gca, 'XTickLabel',{'3','7','4','6' });
set(gca, 'Box', 'off')
ylabel('log2(fold-change)')

subplot(411) % equivalent to line 209
barwitherr(fold_error(i323,:),xticks, fold_change(i323, :))
set(gca, 'XTickLabel',{'3','7','4','6' });
set(gca, 'Box', 'off')
ylabel('log2(fold-change)')

subplot(413) % equivalent to line 210
barwitherr(fold_error(i324,:),xticks, fold_change(i324, :))
set(gca, 'XTickLabel',{'3','7','4','6' });
set(gca, 'Box', 'off')
ylabel('log2(fold-change)')

colormap([0.7 0.7 0.7])

% Save figure
saveas(gcf, fullfile(save_dir,'Figure4_fold_change.eps'), 'epsc')


%% Supplementary figure

%Calculate stripe ratios using 95% values from stripe_95_controls
fold_diff = ones(length(line_numbers_controls), 1);
diff_error = ones(length(line_numbers_controls), 1);

%index line numbers for each enhancer
eve37i = [i204c i325 i326c i327c ci453c i454c];
eve46i = [i214c i328c i329c i330c i455c i456c];
enhancer_dist = [0 0 -500 -500 -1000 -1000];

%calculate stripe ratio within each enhancer
[fold_diff(eve37i) diff_error(eve37i)] = errorprop('divide',...
    stripe_95_controls(eve37i,2), stripe_95_sem_controls(eve37i,2), ...
    stripe_95_controls(eve37i,1), stripe_95_sem_controls(eve37i,1));

[fold_diff(eve46i) diff_error(eve46i)] = errorprop('divide',...
    stripe_95_controls(eve46i,4), stripe_95_sem_controls(eve46i,4), ...
    stripe_95_controls(eve46i,3), stripe_95_sem_controls(eve46i,3));

%use 95% confidence interval for error
diff_error = 1.96*diff_error;

%%
xlimits = [0.5 3.5];
xticks = [1 2 3];
ylimits = [0 2.5];
yticks = [0 1 2];
panel_size = [0.1 0.1 7 10]; %[xLeft yBottom xSize ySize] in cm

fh = figure(3);
fig_intialize(fh, xlimits, xticks, ylimits, yticks, panel_size);

%Eve 3/7
subplot(211) 
barwitherr([diff_error([i453c i326c i204c]), diff_error([i454c i327c i325c])],...
    [fold_diff([i453c i326c i204c]), fold_diff([i454c i327c i325c])])
set(gca, 'XTickLabel',{'-1000bp','-500bp','0bp'});
set(gca, 'Box', 'off')
ylabel('ratio stripe7:stripe3')
legend('forward', 'reverse')

%Eve 4/6
subplot(212) 
barwitherr([diff_error([i455c i329c i214c]), diff_error([i456c i330c i328c])],...
    [fold_diff([i455c i329c i214c]), fold_diff([i456c i330c i328c])])
set(gca, 'XTickLabel',{'-1000bp','-500bp','0bp'});
set(gca, 'Box', 'off')
ylabel('ratio stripe6:stripe4')

colormap([0.8 0.8 0.8; 0.2 0.2 0.2])

% Save figure
saveas(gcf, fullfile(save_dir,'SuppFig_stripe_ratios.eps'), 'epsc')


%%
% figure(3);
% ylimits = [-2 2.5];
% 
% subplot(412) % equivalent to line 207
% barwitherr([fold_error(i240,:); fold_error(i356,:)]',...
%             [fold_change(i240, :); fold_change(i356, :)]')
% ylim(ylimits)
% 
% subplot(414) % equivalent to line 208
% barwitherr([fold_error(i450,:); fold_error(i322,:)]',...
%             [fold_change(i450, :); fold_change(i322, :)]')
% ylim(ylimits)
% 
% subplot(411) % equivalent to line 209
% barwitherr([fold_error(i451,:); fold_error(i323,:)]',...
%             [fold_change(i451, :); fold_change(i323, :)]')
% ylim(ylimits)
% 
% subplot(413) % equivalent to line 210
% barwitherr([fold_error(i452,:); fold_error(i324,:)]',...
%             [fold_change(i452, :); fold_change(i324, :)]')
% ylim(ylimits)
% 
% % Save figure
% saveas(gcf, fullfile(save_dir,'Figure4_fold_change.eps'), 'epsc')


