% This is a wrapper script for analyzing stripe boundaries of the control
% lines for stripe 3/7 and stripe 4/6.

% Calls ProcessStripes to find boundaries for all lines.
% Calls ProcessPeaks to find peaks.

%% Run this part once to generate initial peak guess
% ProcessPeaks_Atlas('/Users/tmartin/Desktop/atlas_20130125_freeze/export/atlas.vpc');

% I then manually renamed the files to use the same names as the folders
% containing pointclouds (eg 'line204')

%% User set values
basedir = '/Users/tmartin/Desktop/pointclouds/';
channel = 'lacZ'; % name of gene
NumStripes = 2; % expected number of stripes

% dirname = {'line204', 'line325', 'line326', 'line327'};  % stripe 3/7 lines
% APthresh = [0.3 0.95]; % portion of ap-axis to search for stripes if processing

% dirname = {'line214', 'line328', 'line329', 'line330'};  % stripe 4/6 lines
% APthresh = [0.5 0.85]; % portion of ap-axis to search for stripes if processing

%% Get stripe boundaries and peak info -> saved to .measure and .peak files
% 
% for i = 1:length(dirname)
%     peak_guess = [];
%     load(fullfile(basedir, [dirname{i} '.peaks']), '-mat');
%     peak_guess = peaks;
%     ProcessEmbryo(fullfile(basedir, dirname{i}), channel, NumStripes, peak_guess, APthresh);
% end

%% Now aggregate data into tables

% Basics
% dirname = {'line204', 'line325', 'line326', 'line327'};  % stripe 3/7 lines
% outputfile = '37measurements.txt';

dirname = {'line214', 'line328', 'line329', 'line330'};  % stripe 4/6 lines
outputfile = '46measurements.txt';

orientation_list = {'F', 'R', 'F', 'R'};    % Orientation code (F=forward, R=reverse)
distance_list = {'P', 'P', 'D', 'D'};       % Distance code (P=proximal, D=distal)

%% Standard characteristics
Nstrips  = 16; %number of DV strips
Nlines   = NumStripes*2; %number of boundaries
filelist = {};


%% Read in list of PCE measurement files from dir_list

for i = 1:length(dirname)
    files = [];
    %get a list of files in a line directory
    files = dir(fullfile(basedir, dirname{i}, 'measure'));
    
    %make a list of files to load
    for j=1:length(files)
        %break down file name into parts
        [~, pcname, t3] = fileparts(files(j).name);
        
        %skip non-measurement files
        if(~strcmp(t3, '.measure'))
            continue  %ie skip files that don't end in .measure
        end
        if(pcname(1) == '.')
            continue
        end
        
        %list of measurement files and PC info
        filelist = [filelist; {dirname{i}, pcname,...
            distance_list{i}, orientation_list{i}}];
    end
end

%create holders for results
summary = cell(length(filelist)*16,9); %16 slices per embryo by 8 columns of info

% DV slice index (folded so that dorsal = 1, ventral = 9, and left and
% right sides are symmetric
% dv_index = [1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2];
% replace 'j' below with 'dv_index(j)' to use folded indices

%% load .measure and .peak files into holder variable
for i=1:length(filelist)
    m = load(fullfile(basedir, filelist{i,1}, 'measure', [filelist{i,2} '.measure']), '-mat');
    m2 = load(fullfile(basedir, filelist{i,1}, 'peaks', [filelist{i,2} '.peak']), '-mat');
%    m.pcname;
%    m.cohort;
%    m.numNuc';
%    m.intensity;
%    m.pos;

%    m2.pcname;
%    m2.cohort;
%    m2.pos; <- peak locations
%    m2.vals; <- peak intensities

    % first aggregate boundary info into summary cell array
    for j = 1:Nstrips
        summary((i-1)*16+j,:) = {filelist{i,2},...  % PC Name
            filelist{i,3},...                       % Distance
            filelist{i,4},...                       % Orientation
            m.cohort,...                            % Cohort
            j, ...                                  % DV Slice
            m.pos(j,:), ...                         % Boundaries
            m2.pos(j,:), ...                        % Peaks
            m2.vals(j,:), ...                       % Peak intensities
            (m2.vals(j,2)/m2.vals(j,1))};           % Relative peak intensity
        
    end
    
end

%% Save measurements

% First save boundary info
fid = fopen(fullfile(basedir,outputfile), 'w');
% Write header line
fprintf(fid, 'Name\tDistance\tOrientation\tCohort\tDVslice\t');
fprintf(fid, 'Anterior1\tPosterior1\tAnterior2\tPosterior2\tPeak1\tPeak2\t');
fprintf(fid, 'PeakIntensity1\tPeakIntensity2\tRelPeakIntensity\n');
% Write data
for i = 1:size(summary,1)
    fprintf(fid, '%s\t', summary{i,1}); % Name
    fprintf(fid, '%s\t', summary{i,2}); % Distance
    fprintf(fid, '%s\t', summary{i,3}); % Orientation
    fprintf(fid, '%d\t', summary{i,4}); % Cohort
    fprintf(fid, '%d\t', summary{i,5}); % DV Slice
    fprintf(fid, '%10.4f\t%10.4f\t%10.4f\t%10.4f\t', [summary{i,6}]); % Boundaries
    fprintf(fid, '%10.4f\t%10.4f\t', [summary{i,7}]); % Peaks
    fprintf(fid, '%10.4f\t%10.4f\t', [summary{i,8}]); % Peak intensities
    fprintf(fid, '%10.4f\n', [summary{i,9}]); % Relative peak intensity
end
fclose(fid);

