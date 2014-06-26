 %FIGURE_PREP_TRACES    Used to read in atlas data and generate line traces
% 
% Reads in a Data.mat file, finds the line trace of each pointcloud
% included in atlas and then finds mean and standard deviation of the line
% traces (mimics non-3D data collection)
% 
% Tara Martin
% Aug 21, 2012
% 

%% Load data, extract and reformat relevant channels
load([data_dir 'Data.mat']);

% get list of genes (only LacZ channels)
channels = gene_names;
I = strncmp('LacZ', channels, 4); %index of LacZ channels
channels = channels(I); %remove non-LacZ channels from list

% extract only relevant time points and channels from scaled data structure
% (assumes only stage 5 embyros and no protein stains)
atlas_data = Edata_scaled(I, 4:9,1);

%% Choose cohort(s) and initialize some variables
cohorts = [3];
dv_slice = [7*pi/16 9*pi/16]; % set the DV slice along side
num_samples = 100; % bins along AP axis
sigma = 0.005; % smoothing parameter

% intializing some variables
LacZ_data = [];
LacZ_mean = [];
LacZ_std = [];


%% get line traces from individual embryos and calculate mean and std

for c=1:length(cohorts) %for each cohort
    cht = cohorts(c);
    for i = 1:length(channels)  % for each gene/channel
        for p = 1:size(atlas_data{i,cht},1) % for each pointcloud in that cohort and channel
            
            % make a pointcloud, orient it
            tempPC = pointcloud(X{cht+3}, atlas_data{i,cht}(p,:));  %X{c+4} contains x,y,z coordinates
            tempPC = egglengthnormalize(tempPC);
            tempPC = rotation(tempPC,pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
            tempPC = align(tempPC); %transform x,y,z coordinates to match a/p-phi adjustments above
            
            % extract line trace for single embryo
            trace = getapprojection(tempPC([0,1],dv_slice), num_samples, sigma, false);
            
            % store traces
            LacZ_data{i,c}(p,:) = trace;  % i=gene, c=cohort
        end
        
        % take mean and std dev of each gene per cohort
        LacZ_mean(i,:,c) = mean(LacZ_data{i,c});
        LacZ_std(i,:,c) = std(LacZ_data{i,c});
        
        %%%% 
        % insert other statistics here
        %%%%
    end
end









%% OLD CODE

% %% (old) normalize the std dev measuremtents by cohort
% 
% load([data_dir 'atlas_data.mat']);
% 
% % reformat data to the form: atlas_exp(nucleus, channel, cohort)
% atlas_exp = Emean_scaled(:,I,4:9,1); % extract relevant expression channels
% 
% maxgainfit = repmat(max(gainfit.*(ct>0),[],2),[1 size(Estd, 3) 1]);
% sgainfit = gainfit ./ (maxgainfit+(maxgainfit==0));
% Estd_scaled = zeros(size(Estd));
% 
% for gg = 1:size(Estd,2);  % each gene
%     for cht = 1:size(Estd, 3);  % each cohort
%         p = prctile(Emean(:,gg,cht,1),99); %99 prctile
%         Estd_scaled(:,gg,cht,1) = Estd(:,gg,cht,1)*sgainfit(gg,cht,1)/(p+(p==0)); % do the actual scaling
%     end
% end
% 
% atlas_std = Estd_scaled(:,I,4:9,1); % extract relevant std dev channels
% 
% LacZ_expr = [];
% LacZ_std = [];
% 
% %% Get expression traces
% for c=1:length(cohorts)
%     cht = cohorts(c);
%     % extract expression values for each channel
%     for i = 1:length(channels)
%         tempPC = pointcloud(X{cht+3}, atlas_exp(:,i,cht));  %X{c+4} contains x,y,z coordinates
%         tempPC = egglengthnormalize(tempPC);
%         tempPC = rotation(tempPC,pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
%         tempPC = align(tempPC); %transform x,y,z coordinates to match a/p-phi adjustments above
%         trace = getapprojection(tempPC([0,1],dv_slice), num_samples, sigma, false);
%         LacZ_expr(i,:,c) = trace;
%     end
% end
% 
% LacZ_expr = squeeze(LacZ_expr); % if only one cohort get rid of extra dimension
% 
% 
% %% Get std dev traces
% for c=1:length(cohorts)
%     cht = cohorts(c);
%     % extract expression values for each channel
%     for i = 1:length(channels)
%         tempPC = pointcloud(X{cht+3}, atlas_std(:,i,cht));  %X{c+4} contains x,y,z coordinates
%         tempPC = egglengthnormalize(tempPC);
%         tempPC = rotation(tempPC,pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
%         tempPC = align(tempPC); %transform x,y,z coordinates to match a/p-phi adjustments above
%         trace = getapprojection(tempPC([0,1],dv_slice), num_samples, sigma, false);
%         LacZ_std(i,:,c) = trace;
%     end
% end
% 
% LacZ_std = squeeze(LacZ_std); % if only one cohort get rid of extra dimension
% 

