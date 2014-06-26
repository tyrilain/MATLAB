%NORMALIZEBYHKB_ATLAS  Extract means of cells expressing hkb and other probe
%(assumed to be lacZ), as well as returning pces with background subtracted.
%Uses thresholded cells to calculate statistics.
%
%   [CORRECTED_PCES HKB MIDDLE] = NormalizeByHkb(dirname, pcelist)
%   CORRECTED_PCES contains atlas data with background subtracted.
%   HKB is a matrix containing: 
%       [geometric_mean ant_mean post_mean nCells_ant nCells_post]
%   MIDDLE is a matrix containing:
%       [middle_mean middle_95 middle_99 middle_nCells]
%
%
% Tara Martin, 2013-02-04
% adapted from NORMALIZEBYHKB


function [corrected_pces hkb middle] = NormalizeByHkb_atlas(pcedata,coords)

%% Initiating holder variables
nPces = size(pcedata,1);

hkb_means = zeros(nPces, 3);
nHkbCells = zeros(nPces, 2);
% hkb_stds = zeros(nPces, 2);
% hkb_medians = zeros(1, nPces);
% hkb_95 = zeros(1, nPces);

middle_means = zeros(nPces, 1);
middle_95 = zeros(nPces, 1);
middle_99 = zeros(nPces, 1);
nMiddleCells = zeros(nPces, 1);
% middle_stds = zeros(1, nPces);
% middle_medians = zeros(1, nPces);

corrected_pces = pcedata; % will subtract background from atlas data


%% Calculate normalization

% I = [1 2];
I = 1:nPces;

for ii = I
%     fprintf('pce number %d', ii);
    %find threshold of expression values
    t = findThreshold(pcedata(ii,:));  %threshold=mode+std
    
    %subtract background
%     background = mean(pce(pce<t)); %define background as mean of "off" cells
    background = t-std(pcedata(ii,:)); %define background as mode of fluorescence distribution
    t = t - background;
    pce = pointcloud(coords, (pcedata(ii,:) - background));
    pce = alignpce(pce);
    
    %find thresholded hkb expression statistics
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

    %save hkb means
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
    
    %% Save pces with background subtracted
    corrected_pces(ii,:) = double(pce);
    %getapprojection(pce([0,1],dv_slice), 100,0.005, false);

end
fprintf('%d embryos\n', ii)
hkb = [hkb_means nHkbCells];
middle = [middle_means  middle_95 middle_99 nMiddleCells];
