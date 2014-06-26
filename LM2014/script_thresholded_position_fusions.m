% Testing whether positional information has changed using thresholded
% expression patterns and Fisher's Exact Test
%
% Tara Martin, 2014-02-09

%%
base_dir = '/Users/tmartin/'; % set manually because changes depending on computer
save_dir = [base_dir 'Dropbox/results/LM2014/position/'];

%tm-36
data_dir = fullfile(base_dir, 'Documents/atlas_20140128_tm36/');
tm36 = load([data_dir 'Data.mat']);

%tm-35
% data_dir_controls = fullfile(base_dir, 'Documents/atlas_20140209_tm35/');
% tm35 = load([data_dir_controls 'Data.mat']);

%%
atlas_data = tm36.Edata([2:21],6,1);  %remove ftz and other cohorts

% Comparing against additive model, so sum line 204 and 214
% control_data = tm36.Edata(2,6,1)+tm36.Edata(7,6,1);

%% Initialize variables

nEmbryos = zeros(20,1);
cmap = cbrewer('seq', 'Greys', 9); %get a grey-scale color map
cmap = cmap(end:-1:1,:); %flip the map so it goes in the other direction

%% Trim embryo to only analyze 0.2-0.9 EL

tempPC = pointcloud(tm36.X{6}, ones(1,6078));   %make a dummy pc
tempPC = alignpce(tempPC);            %normalize and rotate
coords = tempPC.coords;                         %extract new coordinates
I_cells = find(coords(:,1)>0.2 & coords(:,1)<0.9); %get indices of kept cells
coords = coords(I_cells,:);                     %hold onto coordinates

nCells = length(I_cells);                       %how many cells are left

%% Threshold atlas data

for gg = 1:20  %20 genes in tm36
    nEmbryos(gg) = size(atlas_data{gg},1);      %number of embryos for this gene
    embryo_thresh = zeros(nEmbryos(gg),nCells); %holder variable
    
    % threshold each embryo individually
    for i = 1:nEmbryos(gg)
        embryo = atlas_data{gg}(i,I_cells);
        t = findThreshold(embryo);
        embryo_thresh(i,:) = embryo>t;  %save thresholded values
    end
    
    atlas_thresh{gg} = embryo_thresh;   %save thresholded values
end


%% Compare lines using fexact (Fisher's Exact Test)

% TODO: make a model of additive 3/7 and 4/6.  I think this will involve
% hkb normalizing embryos, then randomly sampling and adding embryos from
% each line to generate an "additive" dataset.  For now I'm just going to
% compare to line 240 (gg=7) because it has all four stripes nicely
% defined and nEmbryos=17, which is pretty high.

control_line = 7;

p = zeros(20,nCells);

for ii = 1:20
        x = [atlas_thresh{control_line}; atlas_thresh{ii}];
        y = [zeros(nEmbryos(control_line),1); ones(nEmbryos(ii),1)];
        p(ii,:) = fexact(x,y, 'perm', 10);
end

%% Plot results to see where changes have occured

line_order = [2 13 7 18;
              3 10 14 20;
              4 11 15 17;
              5 12 16 19];

figure(1);

for i = 1:16
    subplot(4,4,i);
    pce = pointcloud(coords, mean(atlas_thresh{line_order(i)}));
    disp(pce, 'print','quick');
    colormap(cmap)
    caxis([0 1]);
end

figure(2);

for i= 1:16
    subplot(4,4,i)
    pce = pointcloud(coords, p(line_order(i),:));
    pce(pce>0.05) = 1;
    disp(pce, 'print','quick');
    colormap(cmap);
    caxis([0 1]);
    % TODO: I can't figure out how to plot the p-values using a different
    % colormap than the mean pces. colormap(gca, cmap) still changes all
    % the colormaps in the figure.
end


