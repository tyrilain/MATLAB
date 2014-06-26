% PROCESSSTRIPES   Computes stripe locations, intensities and number of 
% nuclei for each PCE and saves to individual .measure files.
%
% Generally, called by a summary script.
%
% PROCESSSTRIPES(PCE_DIR, 'CHANNEL', NUMSTRIPES, [APmin APmax]) reads in 
% all pointclouds from PCE_DIR and computes the locations of NUMSTRIPES 
% stripes in 'CHANNEL' between APmin and APmax (normalized to a/p length 1). 
% It returns a "results" structure saved in "PCENAME.mat". Note [APmin
% APmax] is an optional input, default is [0.2 .99].
% 
% Saved variables:
%  pcname
%  cohort
%  numNuc -> [total; stripe 1; stripe 2; ...]
%  pos -> (16xNboundaries matrix) [stripe1_anterior stripe1_posterior...
%                                  stripe2_anterior stripe2_posterior...]
%  intensities -> (Nstripesx3) [stripe 1 95%, mean, std; 
%                               stripe 2 95%, mean, std;...]
%
% Tara Martin, Feb 2012
% Adapted from KRSTRIPES which is adapted from DEMO_EVESTRIPES 
%  (Cris Luengo, Jan-Feb 2009)

function ProcessStripes(dirname, channel, Nstripes, APthresh)

if nargin<4
    APmin = 0.2;
    APmax = 0.99;
else
    APmin = APthresh(1);
    APmax = APthresh(2);
end

% Initialize stuff
Npixels = 400; %number of bins along AP-axis that we will sample
polarities = repmat([1,-1],1,Nstripes);

skipDB = fullfile(dirname, 'skippedembryos.txt');
skipped = {}; %list of skipped embryos
fig1 = 1;
fig2 = 2;

%% Read in list of PCEs

files = dir(dirname); %get a list of files in the directory

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

    % Get the 16 a/p strips
    [onedim,~,phi] = extractpattern(pattern,Npixels);
    
    %% display the normalized pce
    figure(fig1);
    subplot(2,1,1);
    disp(pattern, 'quick');

    %% Try to automatically find the stripes (to a first approximation)
    init = nan(16,Nstripes*2);
    thresholdedPC = segmentgap(pattern,Nstripes, [APmin APmax]);
    
    %% find initial borders
    if max(thresholdedPC)==Nstripes   % Did we get correct number stripes?
        stripes = getstrips(thresholdedPC,phi);
        for ii=1:16
            s = measure(stripes{ii},[],'APextent');
            for jj=s.id
                init(ii,jj*2-1) = s(jj).APextent(1);
                init(ii,jj*2)   = s(jj).APextent(2);
            end
        end
    end
    
    % If automatic domain detection failed, have the user do it
    if all(isnan(init(:,1)))
        fprintf('\nStarting manual stripe finding. \n');
        fprintf('Click on each edge sequentially.  Press backspace to go back or Q to quit. \n');
        storeEmbryo = input('Hit enter to continue on to manual stripe finding or "n" to skip.\n', 's');
        if strcmpi(storeEmbryo, 'n')
            warning('Skipped:noStripe', 'Stripe detection failed, skipping %s', files(i).name)
            skipped = [skipped; {files(i).name}];
            continue
        end
        
        figure(fig2);
        edges = findlines(pattern,polarities);
        if isempty(edges)
            warning('Skipped:noStripe', 'Stripe detection failed, skipping %s', files(i).name)
            skipped = [skipped; {files(i).name}];
            continue
        end
        % Get the initial positions
        for jj=1:Nstripes*2
            s = getstrips(pattern(edges{jj}),phi);
            for ii=1:16
                p = s{ii}.x;
                if ~isempty(p)
                    init(ii,jj) = mean(p);
                    % x-coordinates match a/p coordinates because of align & egglengthnormalize calls above.
                end
            end
        end
    end
    
%     if all(isnan(init(:,1)))
%         warning('Skipped:noStripe', 'Automatic stripe detection failed, skipping %s', files(i).name)
%         skipped = [skipped; {files(i).name}];
%         continue
%     end

    % Check for empty locations (ie fill in initial conditions for strips
    % that did not detect a stripe)
    for jj=1:Nstripes*2
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
    
    %% Get the accurate stripe boundaries
    pos = locateapboundaries(gaussf(onedim,[Npixels/100,0]),init,polarities);
    
    %% Accurately count nuclei and calculate intensities within each stripe
    thresholdedPC(:) = 0;  %clear holder
    numNuc= nan(Nstripes+1,1);
    numNuc(1) = length(pattern); %record total nuclei in embryo
    intensity = nan(Nstripes,3); %columns for max, mean, std
    
    % create holder (thresholdedPC) where nuclei values represent the
    % stripe number they are part of (eg all nuclei in stripe 2 have value
    % of 2, off cells have value of 0)
    for stripe = 1:Nstripes
       for strip = 1:16
           thresholdedPC([pos(strip, (2*stripe-1)), pos(strip, (2*stripe))],...
               [phi(strip), phi(strip+1)]) = stripe;
       end
       %count nuclei in each stripe
       numNuc(stripe+1) = length(thresholdedPC(thresholdedPC==stripe));
       %calculate gene expression intensity values within stripe
       intensity(stripe,1) = prctile(double(pattern(thresholdedPC==stripe)), 95);
       intensity(stripe,2) = mean(pattern(thresholdedPC==stripe));
       intensity(stripe,3) = std(pattern(thresholdedPC==stripe));

    end
    
    %display thresholded pce
%     subplot(3,1,2)
%     thresholdedPC

    %% Display calculated boundaries
    % first get indexing right
    phi = (phi(1:end-1)+phi(2:end))/2;  %finds middle phi value of each strip
    phi = [phi(9:end)-2*pi;phi(1:8)];   %re-order so ventral is middle
    pos = [pos(9:end,:);pos(1:8,:)];    %also re-order position matrix
    
    % then plot
    figure(fig1);
    subplot(2,1,2)
    plot(pos,repmat(phi,1,Nstripes*2),'rx-','linewidth',1);
    xlim([0 1]);
    
    toc %print how long this file took to process
    beep;
 
    % Ask user to approve embryo and store results
    storeEmbryo = input('Use this embryo? (Y/N) (default=Y)', 's');
    
    switch lower(storeEmbryo)
        case {'n', 'no'}
            warning('Skipped:byUser','%s discarded from dataset.', files(i).name)
            skipped = [skipped; {files(i).name}];
        otherwise
            % Collect results
%           results = {pcname,cohort,numNuc,pos,intensity};
            result_file = fullfile(dirname, 'measure', [pcname '.measure']);
            save(result_file, 'pcname', 'cohort', 'numNuc', 'pos', 'intensity');
    end
    
end

%save list of skipped embryos
fid = fopen(skipDB, 'w');
fprintf(fid, '%s\n', skipped{:});
fclose(fid);

