%NORMALIZEBYHKB  Extract means of cells expressing hkb and other probe
%(assumed to be lacZ), as well as returning lateral line traces for each
%pce.  Subtracts background before calculating normalization factors and
%making traces. Uses thresholded cells to calculate statistics.
%
%   [VALUES HKB MIDDLE] = NormalizeByHkb(dirname, pcelist)
%   VALUES is a set of line traces (one per pce).
%   HKB is a matrix containing: 
%       [geometric_mean ant_mean post_mean nCells_ant nCells_post]
%   MIDDLE is a matrix containing:
%       [middle_mean middle_95 middle_99 middle_nCells]
%
%
% Tara Martin, 2013-01-19


function [values hkb middle] = NormalizeByHkb(dirname, pcelist)

dv_slice = [7*pi/16 9*pi/16]; % take a slice along side

%% Initiating holder variables
nPces = length(pcelist);
% embryo_modes = zeros(1, nPces);
pcestage = zeros(1, nPces); %embryo age

hkb_means = zeros(nPces, 3);
% hkb_stds = zeros(nPces, 2);
% hkb_medians = zeros(1, nPces);
% hkb_95 = zeros(1, nPces);
nHkbCells = zeros(nPces, 2);
hkb = zeros(nPces, 5);

middle_means = zeros(nPces, 1);
% middle_stds = zeros(1, nPces);
% middle_medians = zeros(1, nPces);
middle_95 = zeros(nPces, 1);
middle_99 = zeros(nPces, 1);
% middle_1 = zeros(1, nPces);
nMiddleCells = zeros(nPces, 1);
middle = zeros(nPces, 4);

values = zeros(nPces, 100);
% values_hkb = zeros(nPces, 100);
% values_slope = zeros(nPces, 100);
% 
% nTotal = zeros(1, length(line_numbers));
% nTrimmed = zeros(1, length(line_numbers));
% corrTotal = zeros(1, length(line_numbers));
% corrThresholdHkb = zeros(1, length(line_numbers));
% corr95Hkb = zeros(1, length(line_numbers));
% corrTrimmed = zeros(1, length(line_numbers));
% slopeTrimmed = zeros(1, length(line_numbers));
% slopeBoundTrimmed = zeros(1, length(line_numbers));
% intersectTrimmed = zeros(1, length(line_numbers));
% intersectBoundTrimmed = zeros(1, length(line_numbers));
% FvalTrimmed = zeros(1, length(line_numbers));
% nonoutliers = cell(1, length(line_numbers));

%% Calculate normalization

% I = [1 2];
I = 1:length(pcelist);

for ii = I
    %read in and align PC
    filename = fullfile(dirname, 'pointClouds', pcelist{ii});
    pc = readpointcloud(filename);
    pcestage(ii) = pc.metadata.stage_percent;
    
    pce = pointcloud(pc, 'lacZ');
    pce = egglengthnormalize(pce);
    pce = align(pce);
    pce = rotation(pce, pc.metadata.DVrotation+pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
    
    %find threshold of expression values
    t = findThreshold(double(pce));  %threshold=mode+std
    
    %subtract background
%     background = mean(pce(pce<t)); %define background as mean of "off" cells
    background = t-std(double(pce)); %define background as mode of fluorescence distribution
    pce = pce - background;
    t = t - background;
    
    %find hkb expression statistics
    pce_thresh = pce(pce>t);
    hkb_anterior_vals = double(pce_thresh([0 0.15], :));
    hkb_posterior_vals = double(pce_thresh([0.9 1], :));
    
    if isempty(hkb_anterior_vals)
        temp = double(pce([0 0.15], :));
        hkb_anterior_vals = temp(temp>0);
    end
    if isempty(hkb_posterior_vals)
        temp = double(pce([0.9 1], :));
        hkb_posterior_vals = temp(temp>0);
    end
    nHkbCells(ii, :) = [length(hkb_anterior_vals) length(hkb_posterior_vals)];

    %find thresholded hkb values
%     if isempty(hkb_anterior_vals) || isempty(hkb_posterior_vals)
%         hkb_means(ii, ) = [0 0 0];
%     else
    hkb_means(ii, 2) = mean(hkb_anterior_vals); %mean of the anterior

    hkb_means(ii, 3) = mean(hkb_posterior_vals); %mean of the posterior

    hkb_means(ii, 1) = geomean(hkb_means(ii,2:3)); %geometric mean of both
    
%     hkb_stds(ii,1) = std(hkb_anterior_vals);
%     hkb_stds(ii,2) = std(hkb_posterior_vals);
    
%     xx=sort(hkb_anterior_vals);
%     hkb_95(ii) = xx(round(length(xx)*.95));
%     
%     %do normality test of "on" hkb cells
%     if length(hkb_vals(hkb_vals>t)) > 1
%         hkb_normal_threshold(ii) = jbtest(hkb_vals(hkb_vals>t));
%     else
%         hkb_normal_threshold(ii) = 1;
%     end
    
%     %find thresholded hkb values
%     if isempty(hkb_anterior_vals(hkb_anterior_vals>t))
%         hkb_means_threshold(ii) = 0;
%         hkb_stds_threshold(ii) = 0;
%     else
%         hkb_means_threshold(ii) = mean(hkb_anterior_vals(hkb_anterior_vals>t));
%         hkb_stds_threshold(ii) = std(hkb_anterior_vals(hkb_anterior_vals>t));
%     end
    
    %find lacZ expression statistics
    middle_vals = double(pce_thresh([0.15,0.9], :));
    middle_means(ii) = mean(middle_vals);
%     middle_stds(ii) = std(middle_vals);
%     middle_medians(ii) = median(middle_vals);
    xx=sort(middle_vals);
    middle_95(ii) = xx(round(length(xx)*.95));
    middle_99(ii) = xx(round(length(xx)*.99));
%     middle_1(ii) = xx(round(length(xx)*.01));
    nMiddleCells(ii)= length(middle_vals);
    
    %% choose normalization statistic
%     hkb(ii) = hkb_means(ii,1);
    
    %extract pattern for normalization
%     temp = extractpattern(pce,100,0.02, false);
    values(ii,:) = getapprojection(pce([0,1],dv_slice), 100,0.005, false);
%     mean(temp([4 13],:)); %keep mean lateral line trace pattern
    %this version subtracts out background before normalization:
%     values(ii,:) = mean(temp([4 13],:))-min(mean(temp([4 13],:)));
%     values_hkb(ii,:) = values(ii,:)./hkb(ii); %simple divide by hkb
end

hkb = [hkb_means nHkbCells];
middle = [middle_means  middle_95 middle_99 nMiddleCells];
