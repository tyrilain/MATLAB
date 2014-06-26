%% Script to make PCA plots one cohort at a time

% From pcaForTara
% ZW 2012-03-06

%% Choose whether to run fusions or spacer lines:

fusions = 1;
lineNumbers = {'LacZ_0207', 'LacZ_0208', 'LacZ_0209', 'LacZ_0210'};

% fusions = 0;
% lineNumbers = {'LacZ_0356', 'LacZ_0322', 'LacZ_0323', 'LacZ_0324'};

% And a cohort
xgenicCohorts = [6];

atlasFile = '/Users/tmartin/Desktop/atlas_20130125_freeze/export/atlas.vpc';
save_dir = '/Users/tmartin/Dropbox/results/';
atlasData = readpointcloud(atlasFile);

for m=1:length(lineNumbers)
        eval(['a(m,:) = geneValues(atlasData, ''' lineNumbers{m} ''', xgenicCohorts);'])
end

[coeff,score,latent] = princomp(a, 'econ');
b = cumsum(latent)./sum(latent);
if fusions==1
    if xgenicCohorts ==5
        coeff(:,1) = -1*coeff(:,1);
    elseif xgenicCohorts ==6
        coeff(:,[1 3]) = -1*coeff(:,[1 3]);
    elseif xgenicCohorts==9
        coeff(:,2) = -1*coeff(:,2);
    end
    
end

%% Grab stripe boundary info for plotting
data_dir = '~/Desktop/pointclouds';

%% Generate summary.mat files (only once)
% SummaryStripes(fullfile(data_dir, 'line204', 'measure'), 'lacZ', 2);
% SummaryStripes(fullfile(data_dir, 'line214', 'measure'), 'lacZ', 2);

% Load summary.mat files (once)
line204 = load(fullfile(data_dir,'line204/measure/summary.mat'), '-mat');
line214 = load(fullfile(data_dir,'line214/measure/summary.mat'), '-mat');

%% Parameter settings
cohort = xgenicCohorts-3;

% DV angle index
phi = [-pi:(pi/8):pi]';
phi2 = repmat(phi, 1, 2); % use for 2 stripes
phi4 = repmat(phi, 1, 4); % use for 4 stripes

%% Pull out mean and standard error of mean (95% conf interval) for lines
line204.avg = squeeze(line204.boundary_mean(cohort,:,:));
line204.avg = [line204.avg; line204.avg(1,:)];

line214.avg = squeeze(line214.boundary_mean(cohort,:,:));
line214.avg = [line214.avg; line214.avg(1,:)];

%%
% make coordinate matrix for given cohort
c_str = num2str(xgenicCohorts);
eval(['coords = [atlasData.x__' c_str ' atlasData.y__' c_str ' atlasData.z__' c_str '];']);

%%

xSize = 10;  ySize = xSize;  % Set width and height of Figure 4 in cm
% xSize = 5;  ySize = 5;  % Set width and height of supp figures in cm
lims = [min(min(coeff)) max(max(coeff))];

for i = 1:3
    figure(i);
    pca1 = pointcloud(coords, coeff(:,i));
    pca1 = egglengthnormalize(pca1);
    pca1 = rotation(pca1,pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
    pca1 = align(pca1); %transform x,y,z coordinates to match a/p-phi adjustments above

    disp_small(pca1, 'print');  % edit the disp_small script at line 145 to make Fig 4 plots
    % p1 = gca;
    % lims = caxis
%     if fusions == 1
%         switch xgenicCohorts
%             case 5
%                 lims = [-0.0821 0.0670];
%             case 6
%                 lims = [-0.0613 0.0723];
%             case 7
%                 lims = [-0.0554 0.0897];
%             case 8
%                 lims = [-0.0613 0.0889];
%             case 9
%                 lims = [-0.0823 0.0853];
%         end
%     else
%         switch xgenicCohorts
%             case 5
%                 lims = [-0.0511 0.0879];
%             case 6
%                 lims = [-0.0489 0.0688];
%             case 7
%                 lims = [-0.0498 0.0687];
%             case 8
%                 lims = [-0.0594 0.0721];
%             case 9
%                 lims = [-0.0670 0.0738];
%         end
%     end
    set_PCA;
    
    hold on;
    plot(line204.avg, phi4, 'k');
    plot(line214.avg, phi4, 'k');
    if fusions ==1
        saveas(i, fullfile(save_dir,['pca' num2str(i) '.eps']), 'epsc');    % for Figure 4
%         saveas(i, fullfile(save_dir,['pca' num2str(i) '_' c_str '.eps']),'epsc'); % for Supp Figs
    else
        saveas(i, fullfile(save_dir,['pca' num2str(i) '_spacer_' c_str '.eps']), 'epsc')
    end
end

% %%
% figure(2);
% pca2 = pointcloud(coords, coeff(:,2));
% pca2 = egglengthnormalize(pca2);
% pca2 = rotation(pca2,pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
% pca2 = align(pca2); %transform x,y,z coordinates to match a/p-phi adjustments above
% disp(pca2, 'print')
% p2 = gca;
% % lims = caxis
% set_PCA;
% 
% hold on;
% plot(line204.avg, phi4, 'k');
% plot(line214.avg, phi4, 'k');
% 
% %%
% figure(3);
% pca3 = pointcloud(coords, coeff(:,3));
% pca3 = egglengthnormalize(pca3);
% pca3 = rotation(pca3,pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
% pca3 = align(pca3); %transform x,y,z coordinates to match a/p-phi adjustments above
% disp(pca3, 'print')
% p3 = gca;
% % lims = caxis
% set_PCA;
% 
% hold on;
% plot(line204.avg, phi4, 'k');
% plot(line214.avg, phi4, 'k');
% 
% %%
% saveas(1, fullfile(save_dir,['pca1_spacer_' c_str '.eps']), 'epsc')
% saveas(2, fullfile(save_dir,['pca2_spacer_' c_str '.eps']), 'epsc')
% saveas(3, fullfile(save_dir,['pca3_spacer_' c_str '.eps']), 'epsc')

%%
figure(4);
% plot(latent/sum(latent), 's-', 'LineWidth', 2)
% xlim([0.8 3.2]); ylim([0 1]);
bar(latent/sum(latent));
% set(gca, 'XTick', [1 2 3]);
set(gca, 'Box', 'off');
set(gca, 'FontSize', fontsize);
set(gca, 'FontName', fontname);
set(gca, 'LineWidth', 1);
%set(gca, 'YTick', [0 1]);
xlabel('PC');
ylabel('Fraction explained');


xSize = 3;  ySize = 3;  % Set width and height of figure in cm
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0.5 0.5 xSize ySize]);
set(gcf, 'Color', 'w');


if fusions==1
    saveas(gcf, fullfile(save_dir,['pca4_' c_str '.eps']), 'epsc')
else
    saveas(gcf, fullfile(save_dir,['pca4_spacer_' c_str '.eps']), 'epsc')
end

