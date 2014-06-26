%FIGURE_PREP_TRACES2    Used to generate line traces and peak heights
% 
% Using a Data.mat file, finds the line trace of each pointcloud
% included in atlas for each gene in 'channels' and then finds mean and 
% standard deviation of the line traces (mimics non-3D data collection).  
% Also finds 'num_stripes' peak heights within 'search_range' region of 
% embryo. (See figure_controls for usage example).
% 
% Tara Martin
% Nov 09, 2012
% based on FIGURE_PREP_TRACES (Aug 21, 2012)
% 


%% Choose D/V slice and initialize some variables
dv_slice = [7*pi/16 9*pi/16]; % take a slice along side
num_samples = 100;
sigma = 0.005;

LacZ_data = [];
LacZ_peaks = [];
LacZ_mean = [];
LacZ_std = [];
LacZ_numPC = [];

%% get line traces from individual embryos and calculate mean and std

for c=1:length(cohorts) %for each cohort
    for i = 1:length(channels)  % for each gene/channel
        for p = 1:size(atlas_data{i,c},1) % for each pointcloud in that cohort and channel
            
            % make a pointcloud, orient it and extract line trace
            tempPC = pointcloud(X{cohorts(c)+3}, atlas_data{i,c}(p,:));  %X{c+4} contains x,y,z coordinates
            tempPC = egglengthnormalize(tempPC);
            tempPC = rotation(tempPC,pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
            tempPC = align(tempPC); %transform x,y,z coordinates to match a/p-phi adjustments above
            trace = getapprojection(tempPC([0,1],dv_slice), num_samples, sigma, false);
            % now find peaks
%             [trace_max ~] = peakdet(trace(search_range), 0.05);
%             [temp tempI] = sortrows(trace_max, 2); %find the top peaks to ignore background
%             tempI = flipud(tempI);
%             tempI = sort(tempI(1:num_stripes));
%             LacZ_peaks{i,c}(p,:) = [trace_max(tempI,2)'];

            % store traces
            LacZ_data{i,c}(p,:) = trace;
        end
        
        % take mean and std dev of each gene per cohort
        LacZ_mean(i,:,c) = mean(LacZ_data{i,c});
        LacZ_std(i,:,c) = std(LacZ_data{i,c});
        LacZ_numPC(i,c) = p;
    end
end


