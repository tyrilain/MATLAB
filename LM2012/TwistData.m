% TWISTDATA Reads in a set of pointclouds and extracts a set of twist
% expressing cells and lacZ lateral line trace.  It saves the expression
% values to a .mat file named for the line.
% Intended for use with LacZ stains of eve stripe enhancers.
%

function TwistData(dirname, linenames)

channel = 'LacZ';
dv_slice = [7*pi/16 9*pi/16]; % take a slice along side
num_samples = 100;
sigma = 0.005;

%% Read in list of PCEs
for d = 1:length(linenames)
    
    files = dir(fullfile(dirname, linenames{d})); %get a list of files in the directory
    
    ii = 1;
    cohort3.name = {};
    cohort3.twi = {};
    cohort3.lac = {};
    cohort4.name = {};
    cohort4.twi = {};
    cohort4.lac = {};
        
    %% Calculate values for each PCE
    for i=1:length(files)
        tic; %track how long each iteration takes

        %break down file name into parts
        [~, pcname, t3] = fileparts(files(i).name);
        
        %skip non-pce files
        if(~strcmp(t3, '.pce'))
            continue  %ie skip files that don't end in .pce
        end
        if(pcname(1) == '.')
            continue
        end
        
        %% Load the PCE file, assign cohort and check DV rotation
        
        %print out file names as they load
        fprintf('\nProcessing...\t%s\n', files(i).name);
        filename = fullfile(dirname, linenames{d}, files(i).name);
        pc = readpointcloud(filename);
        
        % define cohort PCE belongs to
        if (pc.metadata.stage_percent<4)
            cohort=1;
        elseif (pc.metadata.stage_percent<9)
            cohort=2;
        elseif (pc.metadata.stage_percent<26)
            cohort=3;
        elseif (pc.metadata.stage_percent<51)
            cohort=4;
        elseif (pc.metadata.stage_percent<76)
            cohort=5;
        else
            cohort=6;
        end
        
%         % only using the 3rd cohort
%         if ~(cohort==3)
%             continue;
%         end
        
        %check DV rotation
        if ~isfield(pc.metadata,'DVrotation')
            warning('Skipped:noDV', '%s file needs to be processed by PROCESS1POINTCLOUD first.', files(i).name)
            skipped = [skipped; {files(i).name}];
            continue
        end
        if isempty(pc.metadata.DVrotation)
            warning('Skipped:noDV', '%s does not have a D/V rotation.', files(i).name);
            skipped = [skipped; {files(i).name}];
            continue
        end
        
        %% Align and normalize the PointCloud
        pattern = pointcloud(pc,channel); %load a single channel of pointcloud
        pattern = egglengthnormalize(pattern); %normalize AP length (does not affect x,y,z)
        pattern = rotation(pattern,pc.metadata.DVrotation+pi/2); % Rotate so ventral side is at 0 degrees, rather than 90 degrees.
        pattern = align(pattern); %transform x,y,z coordinates to match a/p-phi adjustments above
        pattern = getneighbors(pattern,'quick','silent'); %fills in neighbors
        
        %% find a threshold to filter out "off" cells
        
        % New thresholding using known "off" region of embryo
        thresh_pattern = pattern([0.05 0.35], [pi/3, 5*pi/3]);
        [N X] = hist(double(thresh_pattern), 100);
        mode_expr = min(X(N==max(N)));
        thresh = mode_expr+2*std(double(thresh_pattern));
        
        % Old method using one std dev from mode of all cells
%         [N X] = hist(double(pattern), 100);
%         mode_expr(ii) = min(X(N==max(N)));
%         thresh(ii) = mode_expr(ii)+std(double(pattern));  

        twi_pattern = pattern([0.05 0.35],:);
        twi_values = double(twi_pattern(twi_pattern>thresh));
        
        lac_values = getapprojection(pattern([0,1],dv_slice), num_samples, sigma, false);
%         lac_pattern = pattern([0.35 0.9], [pi/3, 5*pi/3]);
%         lac_values = double(lac_pattern(lac_pattern>thresh));


        %%
        subplot(1,2,1)
        plot(1:100, lac_values, 1:100, thresh)
        subplot(1,2,2)
        disp(twi_pattern(twi_pattern>thresh), 'quick')
        
        storeEmbryo = input('Use this embryo? (Y/N) (default=Y)', 's');
        switch lower(storeEmbryo)
            case {'n', 'no'}
                continue;
            otherwise
                
                if cohort==3
                    cohort3.name = [cohort3.name files(i).name];
                    cohort3.twi = [cohort3.twi twi_values];
                    cohort3.lac = [cohort3.lac lac_values];
                elseif cohort==4
                    cohort4.name = [cohort4.name files(i).name];
                    cohort4.twi = [cohort4.twi twi_values];
                    cohort4.lac = [cohort4.lac lac_values];
                else
                    warning;
                end
                ii = ii+1; %increment counter
        end
    end
    
    result_file = fullfile(dirname, linenames{d}, linenames{d});
    save(result_file, 'cohort3', 'cohort4');
end
