% Load atlas data and normalize each embryo by mean hkb. This script is
% called once and then data can be used to make a variety of figures.

%% Load atlas data
base_dir = '/Users/tmartin/'; % set manually because changes depending on computer
save_dir = [base_dir 'Dropbox/results/LM2014/'];

%tm-36
data_dir = fullfile(base_dir, 'Documents/atlas_20140210_tm36/');
tm36 = load([data_dir 'Data.mat']);

%tm-35
data_dir_controls = fullfile(base_dir, 'Documents/atlas_20140209_tm35/');
tm35 = load([data_dir_controls 'Data.mat']);

cohorts = 3;
dv_slice = [7*pi/16 9*pi/16]; % take a slice along side

%% Extract the data for lines we care about

% Extract data we care about from atlas

% tm-36 data
line_numbers = {'LacZ:0204','LacZ:0207','LacZ:0208','LacZ:0209',...
    'LacZ:0210','LacZ:0214','LacZ:0240','LacZ:0320','LacZ:0321',...
    'LacZ:0322','LacZ:0323','LacZ:0324','LacZ:0356','LacZ:0450',...
    'LacZ:0451','LacZ:0452','LacZ:0457','LacZ:0458','LacZ:0459','LacZ:0460'};

% Extract only the data we care about
I = zeros(1,length(tm36.gene_names));
for i = 1:length(line_numbers)
    I = I+strcmp(line_numbers{i}, tm36.gene_names);
    I = logical(I);
end
line_numbers = tm36.gene_names(I); % changes order in user defined list to atlas order
atlas_data = tm36.Edata(I, cohorts+3,1); % trimmed atlas data
coords = tm36.X{cohorts+3};  % nuclei X, Y, Z coordinates

% tm-35 data
line_numbers_controls = {'LacZ:0204','LacZ:0214','LacZ:0325','LacZ:0326',...
    'LacZ:0327','LacZ:0328','LacZ:0329','LacZ:0330','LacZ:0453',...
    'LacZ:0454','LacZ:0455','LacZ:0456'};

% Extract only the data we care about
I = zeros(1,length(tm35.gene_names));
for i = 1:length(line_numbers_controls)
    I = I+strcmp(line_numbers_controls{i}, tm35.gene_names);
    I = logical(I);
end
line_numbers_controls = tm35.gene_names(I); % changes order in user defined list to atlas order
atlas_data_controls = tm35.Edata(I, cohorts+3,1); % trimmed atlas data


%% Initialize variables

hkb_mean = cell(length(line_numbers),1);
hkb_mean_controls = cell(length(line_numbers_controls),1);


%% For each line in tm-36, calculate hkb normalization

for gg = 1:length(line_numbers)
%     fprintf('Analyzing line %s\n', line_numbers{gg});

    %for each pce belonging to a line, read in and calculate normalization
    %factors
    [corrected_pces hkb middle] = NormalizeByHkb_atlas(atlas_data{gg}, coords);
%     [I_trimmed] = verifyNorm(dirname, pcestage(I), hkb, middle, []);
        
    % save chosen normalization factor
    % hkb(:,1) is geometric mean of anterior and posterior
    hkb_mean{gg} = hkb(:,1); % save mean hkb values
    
    % save mean lacZ level of line 204 for cross-hybe normalization
    if strcmp('LacZ:0204', line_numbers{gg})
        tm36_mean204 = mean(middle(:,2)./hkb_mean{gg});
        tm36_sem204 = std(middle(:,2)./hkb_mean{gg})/sqrt(length(hkb_mean{gg}));
    end
    
    % save atlas data with background subtracted
    atlas_data{gg} = corrected_pces;
end


%% Now process control lines from hybridization tm-35

for gg = 1:length(line_numbers_controls)
%     fprintf('Analyzing line %s\n', line_numbers_controls{gg});

    %for each pce belonging to a line, read in and calculate normalization
    %factors
    [corrected_pces hkb middle] = NormalizeByHkb_atlas(atlas_data_controls{gg}, coords);
%     [I_trimmed] = verifyNorm(dirname, pcestage(I), hkb, middle, []);
        
    % save chosen normalization factor
    % hkb(:,1) is geometric mean of anterior and posterior
    hkb_mean_controls{gg} = hkb(:,1); % save mean hkb values
    
    % save mean lacZ level of line 204 for cross-hybe normalization
    if strcmp('LacZ:0204', line_numbers_controls{gg})
        tm35_mean204 = mean(middle(:,2)./hkb_mean_controls{gg});
        tm35_sem204 = std(middle(:,2)./hkb_mean_controls{gg})/sqrt(length(hkb_mean_controls{gg})); 
    end
    
    % save atlas data with background subtracted
    atlas_data_controls{gg} = corrected_pces;
end

%% Calculate the cross-hybe normalization factor
% for now I'm using line 204 to normalize everything, but I think it might
% be better to use a geometric mean of 204 and 214.  I just don't know how
% to keep track of SEM with that, so will have to work through that first.
% I doubt it would change the results all that much, but would probably
% make me feel more comfortable.
[cross_hybe cross_hybe_sem] = errorprop('divide', tm36_mean204, tm36_sem204, tm35_mean204, tm35_sem204);

%%
% automatically figure out how to index the appropriate values

%fusions
i207 = find(strcmp('LacZ:0207', line_numbers));
i208 = find(strcmp('LacZ:0208', line_numbers));
i209 = find(strcmp('LacZ:0209', line_numbers));
i210 = find(strcmp('LacZ:0210', line_numbers));

%inversions
i320 = find(strcmp('LacZ:0320', line_numbers));
i321 = find(strcmp('LacZ:0321', line_numbers));

%200bp spacers
i356 = find(strcmp('LacZ:0356', line_numbers));
i322 = find(strcmp('LacZ:0322', line_numbers));
i323 = find(strcmp('LacZ:0323', line_numbers));
i324 = find(strcmp('LacZ:0324', line_numbers));

%1000bp spacers
i240 = find(strcmp('LacZ:0240', line_numbers));
i450 = find(strcmp('LacZ:0450', line_numbers));
i451 = find(strcmp('LacZ:0451', line_numbers));
i452 = find(strcmp('LacZ:0452', line_numbers));

%fusions 1000bp from promoter
i457 = find(strcmp('LacZ:0457', line_numbers));
i458 = find(strcmp('LacZ:0458', line_numbers));
i459 = find(strcmp('LacZ:0459', line_numbers));
i460 = find(strcmp('LacZ:0460', line_numbers));

%single enhancer controls
i204c = find(strcmp('LacZ:0204', line_numbers_controls));
i325c = find(strcmp('LacZ:0325', line_numbers_controls));
i326c = find(strcmp('LacZ:0326', line_numbers_controls));
i327c = find(strcmp('LacZ:0327', line_numbers_controls));
i453c = find(strcmp('LacZ:0453', line_numbers_controls));
i454c = find(strcmp('LacZ:0454', line_numbers_controls));

i214c = find(strcmp('LacZ:0214', line_numbers_controls));
i328c = find(strcmp('LacZ:0328', line_numbers_controls));
i329c = find(strcmp('LacZ:0329', line_numbers_controls));
i330c = find(strcmp('LacZ:0330', line_numbers_controls));
i455c = find(strcmp('LacZ:0455', line_numbers_controls));
i456c = find(strcmp('LacZ:0456', line_numbers_controls));

i204 = find(strcmp('LacZ:0204', line_numbers));
i214 = find(strcmp('LacZ:0214', line_numbers));


%% Trim embryo to only analyze 0.2-0.9 EL

tempPC = pointcloud(tm36.X{6}, ones(1,6078));   %make a dummy pc
tempPC = alignpce(tempPC);                      %normalize length and DV rotation
coords = tempPC.coords;                         %extract new nuclear coordinates
I_cells = find(coords(:,1)>0.2 & coords(:,1)<0.9); %get indices of cells in trunk
coords_trimmed = coords(I_cells,:);             %hold onto trunk coordinates
nCells = length(I_cells);                       %how many cells are left

