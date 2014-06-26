% NORMALIZEBYHKB Reads in a set of pointclouds and computes a
% normalization factor using the hkb co-stain. Intended for use with LacZ
% stains of eve stripe enhancers.
%
%   normalizeByHkb(dirname, pcelist) returns the normalization factor for a
%   set of pces.  dirname is the path to the pce directory and
%   pcelist is a list of pce filenames.
%

function normalizeByHkb(dirname, pcelist)

nPces = length(pcelist);


%% Initiating holder variables

hkb_means = zeros(1, nPces);
hkb_stds = zeros(1, nPces);
hkb_medians = zeros(1, nPces);
hkb_95 = zeros(1, nPces);
nHkbCells = zeros(1, nPces);
hkb_normal_threshold = zeros(1, nPces);
hkb_means_threshold = zeros(1, nPces);
hkb_stds_threshold = zeros(1, nPces);

middle_means = zeros(1, nPces);
middle_stds = zeros(1, nPces);
middle_medians = zeros(1, nPces);
middle_99 = zeros(1, nPces);
middle_1 = zeros(1, nPces);
nMiddleCells = zeros(1, nPces);

% substages = [pcedata.substage];
% names ={pcedata.name};
% dvrots = [pcedata.orientation];

%% find properties of the expression values on the posterior end and middle
%part of the embryo

for ii=1:nPces
    ii
            
    %filename = ['v5_' pcedata(ii).name];
    
    %read in and align PC
    pc = readpointcloud(fullfile(dirname, pcelist{ii}));
    pce = pointcloud(pc,'lacZ', 'correct'); %load a single channel of pointcloud
    pce = egglengthnormalize(pce);
    pce = align(pce);
    pce = rotation(pce,pc.metadata.DVrotation+pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
        
    %find mode and threshold of expression values
    t = findThreshold(double(pce));
    
    %find hkb expression values
    hkb_vals = double(pce([0.9,1], :));  % use the last 10% egg length to calculate hkb expression
    hkb_means(ii) = mean(hkb_vals);
    hkb_stds(ii) = std(hkb_vals);
    hkb_medians(ii) = median(hkb_vals);
    xx=sort(hkb_vals);
    hkb_95(ii) = xx(round(length(xx)*.95));
    nHkbCells(ii) = length(hkb_vals);
    
    %do normality test of "on" hkb cells (probably not normal)
    if length(hkb_vals(hkb_vals>t)) > 1
        hkb_normal_threshold(ii) = jbtest(hkb_vals(hkb_vals>t));
    else
        hkb_normal_threshold(ii) = 1;
    end
    
    %find thresholded hkb values
    if isempty(hkb_vals(hkb_vals>t))
        hkb_means_threshold(ii) = 0;
        hkb_stds_threshold(ii) = 0;
    else
        hkb_means_threshold(ii) = mean(hkb_vals(hkb_vals>t));
        hkb_stds_threshold(ii) = std(hkb_vals(hkb_vals>t));
    end
    
    
    middle_vals = double(pce([0.1,0.9], :)); %exclude ends of embryo that are expressing hkb
    middle_means(ii) = mean(middle_vals);
    middle_stds(ii) = std(middle_vals);
    middle_medians(ii) = median(middle_vals);
    xx=sort(middle_vals);
    middle_99(ii) = xx(round(length(xx)*.99));
    middle_1(ii) = xx(round(length(xx)*.01));
    nMiddleCells(ii)= length(middle_vals);
    
    
%     %find thresholded middle values values
%     if isempty(middle_vals(middle_vals>tMiddle))
%         middle_means_threshold(ii) = 0;
%         middle_stds_threshold(ii) = 0;
%         middle_99_threshold(ii) = 0;
%     else
%         
%         middle_means_threshold(ii) = mean(middle_vals(middle_vals>tMiddle));
%         middle_stds_threshold(ii) = std(middle_vals(middle_vals>tMiddle));
%         
%         xx=sort(middle_vals(middle_vals>tMiddle));
%         middle_99_threshold(ii) = xx(round(length(xx)*.99));
%     end
end

