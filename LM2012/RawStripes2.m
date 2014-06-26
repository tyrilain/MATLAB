% This script is used to aggregate the raw measurements from a bunch of
% pce.measure files.  Outputs a table to use in an external statistics
% program (JMP).
%
% Tara Martin, April 2012

clear all;

% Basics
% dir_list = {'line204', 'line325', 'line326', 'line327'}; %lines to analyze
% outputfile1 = '37boundaries.txt';
% outputfile2 = '37nuclei.txt';

dir_list = {'line214', 'line328', 'line329', 'line330'};
outputfile1 = '46boundaries.txt';
outputfile2 = '46nuclei.txt';

orientation_list = {'F', 'R', 'F', 'R'};    % Orientation code (F=forward, R=reverse)
distance_list = {'P', 'P', 'D', 'D'};       % Distance code (P=proximal, D=distal)
NumStripes = 2;

%% Standard characteristics
Nstrips  = 16; %number of DV strips
Nlines   = NumStripes*2; %number of boundaries
filelist = {};

%% Read in list of PCE measurement files from dir_list

for dir_i = 1:length(dir_list)
    files = [];
    files = dir(dir_list{dir_i}); %get a list of files from the directory
    
    %make a list of files to load
    for i=1:length(files)
        %break down file name into parts
        [~, pcname, t3] = fileparts(files(i).name);
        
        %skip non-measurement files
        if(~strcmp(t3, '.measure'))
            continue  %ie skip files that don't end in .measure
        end
        if(pcname(1) == '.')
            continue
        end
        
        %list of measurement files and PC info
        filelist = [filelist; {dir_list{dir_i}, pcname,...
            distance_list{dir_i}, orientation_list{dir_i}}];
    end
end

%create holders for results
summary = cell(length(filelist)*16,6); %16 slices per embryo by 6 columns of info
summary2 = cell(length(filelist), 7);  % number embryos by 6 columns of info

%load measurement files into holder variable
for i=1:length(filelist)
    m = load(fullfile(filelist{i,1}, [filelist{i,2} '.measure']), '-mat');
%    m.pcname;
%    m.cohort;
%    m.numNuc';
%    m.intensity;
%    m.pos;

    % first aggregate boundary info into summary cell array
    for j = 1:Nstrips
        summary((i-1)*16+j,:) = {filelist{i,2},...  % PC Name
            filelist{i,3},...                       % Distance
            filelist{i,4},...                       % Orientation
            m.cohort,...                            % Cohort
            j,...                                   % DV Slice
            m.pos(j,:)};                            % Boundaries
    end
    
    % then aggregate nuclei and intensity info into cell array
    summary2(i,:) = {filelist{i,2},...              % PC Name
            filelist{i,3},...                       % Distance
            filelist{i,4},...                       % Orientation
            m.cohort,...                            % Cohort
            m.numNuc(2:3),...                       % Number nuclei
            m.intensity(:,1),...                    % Stripe intensity (95%)
            m.intensity(1,1)/m.intensity(2,1)};
end



%% Save measurements

% First save boundary info
fid = fopen(outputfile1, 'w');
fprintf(fid, 'Name\tDistance\tOrientation\tCohort\tSlice\tAnterior1\tPosterior1\tAnterior2\tPosterior2\n');
for i = 1:length(filelist)*16
    fprintf(fid, '%s\t', summary{i,1}); % Name
    fprintf(fid, '%s\t', summary{i,2}); % Distance
    fprintf(fid, '%s\t', summary{i,3}); % Orientation
    fprintf(fid, '%d\t', summary{i,4}); % Cohort
    fprintf(fid, '%d\t', summary{i,5}); % DV Slice
    fprintf(fid, '%10.4f\t%10.4f\t%10.4f\t%10.4f\n', [summary{i,6}]); % Boundaries
end
fclose(fid);

% Next save nuclei and intensity info

fid = fopen(outputfile2, 'w');
fprintf(fid, 'Name\tDistance\tOrientation\tCohort\tS1_nuc\tS2_nuc\tS1_intensity\tS2_intensity\tRelIntensity\n');
for i = 1:length(filelist)
    fprintf(fid, '%s\t', summary2{i,1}); % Name
    fprintf(fid, '%s\t', summary2{i,2}); % Distance
    fprintf(fid, '%s\t', summary2{i,3}); % Orientation
    fprintf(fid, '%d\t', summary2{i,4}); % Cohort
    fprintf(fid, '%d\t%d\t', [summary2{i,5}]); % Number nuclei each stripe
    fprintf(fid, '%10.4f\t%10.4f\t', [summary2{i,6}]);   % Intensity (95%) each stripe
    fprintf(fid, '%10.4f\n', summary2{i,7});   % Relative stripe intensity
end
fclose(fid);
