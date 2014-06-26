% Tara 20120301
% Read in an atlas, extract channels to get an expression trace for each
% line.

% get list of genes (only LacZ channels)
channels = {atlas.metadata.column_info.name};
I = strncmp('LacZ', channels, 4); %index of LacZ channels
channels = channels(I); %remove non-LacZ channels from list

LacZ_expr = [];
    
for c=1:length(cohorts)

    % Select cohort to use
    atlas_cohort = num2str(cohorts(c)+3);  % 6 is 3rd cohort

    % make coordinates matrix for given cohort
    eval(['x = atlas.x__' atlas_cohort ';']);
    eval(['y = atlas.y__' atlas_cohort ';']);
    eval(['z = atlas.z__' atlas_cohort ';']);
    coords = [x y z];

    % extract expression values for each channel
    for i = 1:length(channels)
        eval(['temp = atlas.' channels{i} '__' atlas_cohort ';']);
        tempPC = pointcloud(coords, temp);
        tempPC = egglengthnormalize(tempPC);
        tempPC = rotation(tempPC,pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
        tempPC = align(tempPC); %transform x,y,z coordinates to match a/p-phi adjustments above
        dv_slice = [7*pi/16 9*pi/16]; % take a slice along side
        trace = getapprojection(tempPC([0,1],dv_slice), 100,0.005, false);
        LacZ_expr(i,:,c) = trace;
    end
end

squeeze(LacZ_expr); % if only one cohort get rid of extra dimension

% %% read in atlas
% atlas = readpointcloud(atlas_dir);
% 
% %% Select cohort to use
% atlas_cohort = num2str(cohort+3);  % 6 is 3rd cohort
% 
% % make coordinates matrix for given cohort
% eval(['x = atlas.x__' atlas_cohort ';']);
% eval(['y = atlas.y__' atlas_cohort ';']);
% eval(['z = atlas.z__' atlas_cohort ';']);
% coords = [x y z];
% 
% %% get list of genes (only LacZ channels)
% channels = {atlas.metadata.column_info.name};
% I = strncmp('LacZ', channels, 4); %index of LacZ channels
% channels = channels(I); %remove non-LacZ channels from list
% 
% % LacZ_expr = [];
% % for i = 1:length(channels)
% %     eval(['temp = atlas.' channels{i} '__' cohort ';']);
% %     LacZ_expr(:, i) = temp;
% % end
% 
% %% extract expression values for each channel
% LacZ_expr = [];
% for i = 1:length(channels)
%     eval(['temp = atlas.' channels{i} '__' atlas_cohort ';']);
%     tempPC = pointcloud(coords, temp);
%     tempPC = egglengthnormalize(tempPC);
%     tempPC = rotation(tempPC,pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
%     tempPC = align(tempPC); %transform x,y,z coordinates to match a/p-phi adjustments above
%     dv_slice = [7*pi/16 9*pi/16];
%     trace = getapprojection(tempPC([0,1],dv_slice), 100,0.005, false);
%     LacZ_expr(i,:) = trace;
% end

% %% Fancy plots
% xaxis = 30:100;
% figure;
% 
% for i = 1:length(channels)
%     subplot(length(channels),1,i)
%     plot(xaxis, LacZ_expr(i,xaxis))
%     title(channels{i})
% end

%% testing area


