% Testing whether positional information has changed using thresholded
% expression patterns and Fisher's Exact Test
%
% Tara Martin, 2014-02-09
% Edits:
% 2014-03-21: 


%%

myfiguredefaults; %set a variety of figure defaults

base_dir = '/Users/tmartin/'; % set manually because changes depending on computer
save_dir = [base_dir 'Dropbox/results/LM2014/position/'];

%tm-35
data_dir_controls = fullfile(base_dir, 'Documents/atlas_20140209_tm35/');
tm35 = load([data_dir_controls 'Data.mat']);
atlas_data_controls = tm35.Edata(2:13,6,1);  %remove ftz and other cohorts

%% Initialize variables

nEmbryos = zeros(12,1);

%% Trim embryo to only analyze 0.2-0.9 EL

tempPC = pointcloud(tm35.X{6}, ones(1,6078));   %make a dummy pc
tempPC = alignpce(tempPC);                      %normalize length and DV rotation
coords = tempPC.coords;                         %extract new coordinates
I_cells = find(coords(:,1)>0.2 & coords(:,1)<0.9); %get indices of kept cells
coords_trimmed = coords(I_cells,:);                     %hold onto coordinates
nCells = length(I_cells);                       %how many cells are left

%% Threshold atlas data
for gg = 1:12  %12 genes in tm35
    nEmbryos(gg) = size(atlas_data_controls{gg},1);      %number of embryos for this gene
    embryo_thresh = zeros(nEmbryos(gg),nCells); %holder variable
    
    % threshold each embryo individually
    for i = 1:nEmbryos(gg)
        embryo = atlas_data_controls{gg}(i,I_cells); %only use middle cells
        t = findThreshold(embryo);
        embryo_thresh(i,:) = embryo>t;  %save thresholded values
    end
    
    atlas_thresh{gg} = embryo_thresh;   %save thresholded values
end


%% Compare lines using fexact (Fisher's Exact Test)

% 3/7 enhancers: compare each line to each other
% lines 204, 325, 326, 327, 453, and 454
I37 = [1 3 4 5 9 10];
p37 = zeros(6,6,nCells);

for ii = 1:6
    for jj = 1:6
        x = [atlas_thresh{I37(ii)}; atlas_thresh{I37(jj)}];
        y = [zeros(nEmbryos(I37(ii)),1); ones(nEmbryos(I37(jj)),1)];
        
        %use 'perm' option to shuffle data 10 times to correct for multiple
        %hypothesis testing
        p37(ii,jj,:) = fexact(x,y, 'perm', 10); %calculates p-value for each cell
    end
end

%% Plot results to see where changes have occured

%plot parameters
cmap = cbrewer('seq', 'Blues', 9); %get a grey-scale color map
red = [1 0.1 0.1];
cmap = [cmap(1:6,:); red]; % drop darkest blue and add red at top end
s = 0.1; %size in pts of each cell plotted
ap_range = [0.2 0.9]; %this is range in which cells were tested for significance
panel_size = [0 0 3 3]; %dimensions in cm

forward_plot_order = [1 3 5];
reverse_plot_order = [2 4 6];
distances = {'0', '500', '1000'};

for i = 1:3
    %% Figure 1A (eve 3/7 in forward orientation)
    figure(i)
    set(gcf,'PaperPosition', panel_size); %figure dimensions in cm
    
    %assemble a pointcloud of average on/off values
    pce = pointcloud(coords_trimmed, mean(atlas_thresh{I37(forward_plot_order(i))}));
    %set values of p<0.05 cells to 2 so they are plotted in different color
    pce(squeeze(p37(1,forward_plot_order(i),:))<0.05) = 2;
    %remove cells with no expression to simplify plotting
    pce = pce(pce>0);
    %use custom plotting function
    plotpce(pce, ap_range, s);

    colormap(cmap);
    caxis([0 2]);
    set(gca,'FontSize', 6); % use 6pt font
    set(gca,'LineWidth', 0.5); % 0.5pt axis line width
    set(gca,'xtick', [], 'ytick',[]);
    
    % Save figure
    figname = sprintf('fig2A_%s.eps', distances{i});
    saveas(gcf, fullfile(save_dir,figname), 'epsc')
    
    %% Figure 1B (eve 3/7 in reverse orientation)
    figure(i+3)
    set(gcf,'PaperPosition', panel_size); %figure dimensions in cm
    
    %assemble a pointcloud of average on/off values
    pce = pointcloud(coords_trimmed, mean(atlas_thresh{I37(reverse_plot_order(i))}));
    %set values of p<0.05 cells to 2 so they are plotted in different color
    pce(squeeze(p37(1,reverse_plot_order(i),:))<0.05) = 2;
    %remove cells with no expression to simplify plotting
    pce = pce(pce>0);
    %use custom plotting function
    plotpce(pce, ap_range, s);

    colormap(cmap);
    caxis([0 2]);
    set(gca,'FontSize', 6); % use 6pt font
    set(gca,'LineWidth', 0.5); % 0.5pt axis line width
    set(gca,'xtick', [], 'ytick',[]);

    % Save figure
    figname = sprintf('fig2B_%s.eps', distances{i});
    saveas(gcf, fullfile(save_dir,figname), 'epsc')

end

% Save figure
% saveas(gcf, fullfile(save_dir,'stripes37_thresholded.eps'), 'epsc')


%%
% 4/6 enhancers: compare each line to each other
% lines 214, 328, 329, 330, 455, 456
I46 = [2 6 7 8 11 12];
p46 = zeros(6,6,nCells);

for ii = 1:6
    for jj = 1:6
        x = [atlas_thresh{I46(ii)}; atlas_thresh{I46(jj)}];
        y = [zeros(nEmbryos(I46(ii)),1); ones(nEmbryos(I46(jj)),1)];
        p46(ii,jj,:) = fexact(x,y, 'perm', 10);
    end
end

%% Plot results
%use same plot parameters as for eve 3/7 panels

for i = 1:3
    %% Figure 1A (eve 3/7 in forward orientation)
    figure(i)
    set(gcf,'PaperPosition', panel_size); %figure dimensions in cm
    
    %assemble a pointcloud of average on/off values
    pce = pointcloud(coords_trimmed, mean(atlas_thresh{I46(forward_plot_order(i))}));
    %set values of p<0.05 cells to 2 so they are plotted in different color
    pce(squeeze(p46(1,forward_plot_order(i),:))<0.05) = 2;
    %remove cells with no expression to simplify plotting
    pce = pce(pce>0);
    %use custom plotting function
    plotpce(pce, ap_range, s);

    colormap(cmap);
    caxis([0 2]);
    set(gca,'FontSize', 6); % use 8pt font
    set(gca,'LineWidth', 0.5); % 0.5pt axis line width
    set(gca,'xtick', [], 'ytick',[]);

    % Save figure
    figname = sprintf('fig2C_%s.eps', distances{i});
    saveas(gcf, fullfile(save_dir,figname), 'epsc')
    
    %% Figure 1B (eve 3/7 in reverse orientation)
    figure(i+3)
    set(gcf,'PaperPosition', panel_size); %figure dimensions in cm
    
    %assemble a pointcloud of average on/off values
    pce = pointcloud(coords_trimmed, mean(atlas_thresh{I46(reverse_plot_order(i))}));
    %set values of p<0.05 cells to 2 so they are plotted in different color
    pce(squeeze(p46(1,reverse_plot_order(i),:))<0.05) = 2;
    %remove cells with no expression to simplify plotting
    pce = pce(pce>0);
    %use custom plotting function
    plotpce(pce, ap_range, s);

    colormap(cmap);
    caxis([0 2]);
    set(gca,'FontSize', 6); % use 8pt font
    set(gca,'LineWidth', 0.5); % 0.5pt axis line width
    set(gca,'xtick', [], 'ytick',[]);

    % Save figure
    figname = sprintf('fig2D_%s.eps', distances{i});
    saveas(gcf, fullfile(save_dir,figname), 'epsc')

end


% Save figure
% saveas(gcf, fullfile(save_dir,'stripes46_thresholded.eps'), 'epsc')


%% Make panel for Figure 1B

figure()
set(gcf,'PaperPosition', [0 0 3.5 3.2]); %figure dimensions in cm
pce = pointcloud(coords, mean(atlas_data_controls{I46(1)}));
pce = pce(pce>0.1);
pce = pce([0.2 .85],:,:);
plotpce(pce, [0 1], 1);
colormap(cbrewer('seq', 'Reds', 6));

set(gca,'FontSize', 6); % use 8pt font
set(gca,'LineWidth', 0.5); % 0.5pt axis line width
set(gca,'xtick', [], 'yticklabel',{'dorsal','ventral','dorsal'});


% Save figure 
%(have nice version so making sure it doesn't get overwritten 2014-04-10)
% saveas(gcf, fullfile(save_dir,'fig1B.eps'), 'epsc')


%% Plot results to see where changes have occured
% 
% cmap2 = cmap(end:-1:1,:); %flip the map so it goes in the other direction
% 
% figure(3);
% % subplot(2,3,1);
% % pce = pointcloud(coords, mean(atlas_thresh{1}));
% % disp(pce, 'print','quick');
% 
% plot_order = [1 3 5 2 4 6];
% 
% for i = 2:6
%     %  subplot(6,6,6*(ii-1)+jj)
%     subplot(2,3,i);
% %     pce = pointcloud(coords, mean(atlas_thresh{I37(plot_order(i))}));
% %     disp(pce, 'print');
% %     subplot(2,6,i+6)
%     pce = pointcloud(coords, squeeze(p37(1,plot_order(i),:)));
%     pce(pce>0.05) = 1;
%     disp(pce, 'print','quick');
%     colormap(cmap2);
%     caxis([0 1]);
%     title(tm35.gene_names(I37(plot_order(i))+1))
%     % TODO: I can't figure out how to plot the p-values using a different
%     % colormap than the mean pces. colormap(gca, cmap) still changes all
%     % the colormaps in the figure.
% end
% 
% % Save figure
% saveas(gcf, fullfile(save_dir,'stripes37_pvalues.eps'), 'epsc')
% 
% %% Plot results to see where changes have occured
% 
% figure(4);
% 
% for i = 2:6
%     subplot(2,3,i);
%     pce = pointcloud(coords, squeeze(p46(1,plot_order(i),:)));
%     pce(pce>0.05) = 1;
%     disp(pce, 'print','quick');
%     colormap(cmap2);
%     caxis([0 1]);
%     title(tm35.gene_names(I46(plot_order(i))+1))
%     % TODO: I can't figure out how to plot the p-values using a different
%     % colormap than the mean pces. colormap(gca, cmap) still changes all
%     % the colormaps in the figure.
% end
% 
% % Save figure
% saveas(gcf, fullfile(save_dir,'stripes46_pvalues.eps'), 'epsc')

