% PROCESSSPEAKS   Computes stripe peak locations for each PCE and saves to 
% individual .peak files.
%
% Generally, called by a summary script.
%
% PROCESSSTRIPES(PCE_DIR, 'CHANNEL', NUMSTRIPES, MANUAL, INIT) reads in 
% all pointclouds from PCE_DIR and computes the peak locations of NUMSTRIPES 
% stripes in 'CHANNEL'. If MANUAL = 'on', uses INIT as initial guess of
% peaks.
%
% PROCESSSTRIPES(PCE_DIR, 'CHANNEL', NUMSTRIPES) reads in all pointclouds 
% from PCE_DIR and computes the peak locations of NUMSTRIPES stripes in 
% 'CHANNEL'. Tries to find peak initial guess by thresholding pce to find 
% stripe domains.
%
% Function returns a "results" structure saved in "PCENAME.peak".
% 
% Saved variables:
%  pcname
%  cohort
%  pos -> (16xNUMSTRIPES matrix) [stripe1_peak stripe2_peak... ]
%
% Tara Martin, May 2012
% Adapted from PROCESSSTRIPES which is adapted from DEMO_EVESTRIPES 
%  (Cris Luengo, Jan-Feb 2009)

function ProcessPeaks(dirname, channel, Nstripes, peak_guess)

% Initialize stuff
% init = nan(16,Nstripes);
Npixels = 400; %number of bins along AP-axis that we will sample
skipDB = fullfile(dirname, 'skippedembryos_peaks.txt');
skipped = {}; %list of skipped embryos
fig1 = 1;
% fig2 = 2;


%% Read in list of PCEs

files = dir(dirname); %get a list of files in the directory

%% Calculate values for each PCE
for i = 1:length(files);
    tic; %track how long each iteration takes
    
    %% break down file name into parts
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
    filename = fullfile(dirname, files(i).name);
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
    pattern = stretch(pattern,0.1,99.9); %rescales expression values from 0 to 100
    pattern = getneighbors(pattern,'quick','silent'); %fills in neighbors
    
    
    %% display the normalized pce
    figure(fig1);
    subplot(2,1,1);
    disp(pattern, 'quick');
    
    
    %% Find peaks
    done = 0;
    counter = 1;
    %%
    while done == 0
        %% Get the 16 a/p strips
        [onedim,~,phi] = extractpattern(pattern,Npixels);
        
        if counter==1
            % First try provided guess
            init = peak_guess(:,:,cohort);
        elseif counter==2
            % Try to automatically find the stripes (to a first approximation)
            thresholdedPC = segmentpairrule(pattern,Nstripes);
            
            if max(thresholdedPC)==Nstripes   % Did we get correct number stripes?
                stripes = getstrips(thresholdedPC,phi);
                for ii=1:16
                    s = measure(stripes{ii},[],'APextent');
                    for jj=s.id
                        init(ii,jj) = (s(jj).APextent(1)+s(jj).APextent(2))/2;
                    end
                end
                % Check for empty locations (ie fill in initial conditions for strips
                % that did not detect a stripe)
                for jj=1:Nstripes
                    I = isnan(init(:,jj));
                    if any(I)
                        warning('NaNs:found', 'Filling in %d missing values', sum(I))
                        if I(1)
                            % if the first point is NaN, copy the last found point into it
                            tmp = init(~I,jj);
                            init(1,jj) = tmp(end);
                            I(1) = false;
                        end
                        % Now propagate non-NaN values forward into NaN-valued locations
                        for ii=find(I(:)')
                            init(ii,jj) = init(ii-1,jj);
                        end
                    end
                end
                
            else
                warning('Skipped:noStripes','Found wrong number of stripes.\n')
                counter = counter+1;
                continue;
            end
        elseif counter==3
            % Try to manually find stripes, not yet implemented
            counter = counter+1;
            continue;
        else
            warning('Skipped:noStripes','Could not find peaks, %s discarded from dataset.', files(i).name)
            skipped = [skipped; {files(i).name}];
            done = 1;  % don't bother trying to find stripes anymore
            continue;
        end
        

        
        %% Find precise peak location
        onedim_gauss = double(gaussf(onedim,[Npixels/100,0]));
        [pos vals] = locateappeaks(onedim_gauss, init);
            
        %% Display calculated boundaries
        % first get indexing right
        phi = (phi(1:end-1)+phi(2:end))/2;  %finds middle phi value of each strip
        phi = [phi(9:end)-2*pi;phi(1:8)];   %re-order so ventral is middle
        pos = [pos(9:end,:);pos(1:8,:)];    %also re-order position matrix
        vals = [vals(9:end,:);vals(1:8,:)];    %also re-order values matrix
        
        % then plot
        figure(fig1);
        subplot(2,1,2)
        plot(pos,repmat(phi,1,Nstripes),'rx-','linewidth',1);
        xlim([0 1]);
        
%         toc %print how long this file took to process
%         beep;
        
        % Ask user to approve embryo and store results
        beep;
        storeEmbryo = input('Stripes look okay? (Y/N) (default=Y)', 's');
        
        switch lower(storeEmbryo)
            case {'n', 'no'}
                counter = counter+1; % Will try the next method
            otherwise
                % Collect results
                % results = {pcname,cohort,pos};
                result_file = fullfile(dirname, 'peaks', [pcname '.peak']);
                save(result_file, 'pcname', 'cohort', 'pos', 'vals');
                done = 1;
        end
    end
    toc
end

%save list of skipped embryos
fid = fopen(skipDB, 'w');
fprintf(fid, '%s\n', skipped{:});
fclose(fid);