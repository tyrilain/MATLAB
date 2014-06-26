%PCA_FOR_PCES    Evaluate a set of registered PCEs using PCA
% 
% Reads in a Data.mat file and uses Edata_scaled to run PCA on a set of
% individual pointclouds (after registering and scaling as for an atlas).
% Used to compare within line variation to across line variation.
% 
% Tara Martin
% Aug 21, 2012
% 

%% Load data, extract and reformat relevant channels

base_dir = '/Users/tmartin/'; % set manually because changes depending on computer
data_dir = [base_dir 'Dropbox/data/atlas_20120810_gemstat/'];
save_dir = [base_dir 'Dropbox/results/'];

load([data_dir 'Data.mat']);

% get list of included lines
channels = {'LacZ:0207', 'LacZ:0208', 'LacZ:0209', 'LacZ:0210'};

I = false(size(gene_names)); % Index of included lines
for i = 1:length(channels)
    I = I | strcmp(channels{i}, gene_names);
end

% extract only relevant time points and channels from scaled data structure
% (assumes only stage 5 embyros and no protein stains)
atlas_data = Edata_scaled(I, 4:9,1);

%% Pick a cohort, any cohort...
cohorts = [1:6];
a = [];

for i=1:length(channels)
    for c = 1:length(cohorts)
        a = [a; atlas_data{i, cohorts(c)}];
    end
end

%% Run PCA on dataset

[coeff,score,latent] = princomp(a, 'econ');
b = cumsum(latent)./sum(latent);

% plot variance explained by components
figure;
plot(b, '.-', 'LineWidth', 2)
ylim([0 1]);
xlabel('Component');
ylabel('Fraction variation explained');
b(1:10)
% set(gca, 'XTick', [1 2 3]);
% set_fig;
% saveas(gcf, fullfile(save_dir,'pca4.eps'), 'epsc')

%% Visualize components

figure;
disp(pointcloud(X{6}, coeff(:,1)), 'print')

figure;
disp(pointcloud(X{6}, coeff(:,2)), 'print')

figure;
disp(pointcloud(X{6}, coeff(:,3)), 'print')

% figure_params;
% xSize = 8;  ySize = 8;  % Set width and height of figure in cm
% 
% atlasFile = '/Users/tmartin/Dropbox/data/atlas20120222.vpc';
% atlasData = readpointcloud(atlasFile);
% figure;
% disp(pointcloud([atlasData.x__6(1:5000) atlasData.y__6(1:5000) atlasData.z__6(1:5000)], coeff(:,1)), 'print')
