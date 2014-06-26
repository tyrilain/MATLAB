% This script specifies some variables common to all figures

fontsize = 12;
fontname = 'Helvetica';

% commonly used colors
wh = [1 1 1];
bl = [0 0 0];
dark_grey = [0.4 0.4 0.4];
light_grey = [0.7 0.7 0.7];
steel = [70 130 180]/255;

paper = 'off';

%color order from color brewer (10 colors)
brewer10 = [
    31,120,180
    51,160,44
    227,26,28
    255,127,0
    106,61,154
    166,206,227
    178,223,138
    251,154,153
    253,191,111
    202,178,214
    ]./256;

% Zeba's color order:
% colororder = [202 206 194; 233 109 31; 113 111 179; 105 106 108; 181 18 27; 55 3 86]./256;


% Usage example:
%
% figure;
% set(gcf,'color',bl, 'InvertHardCopy', paper); %sets background to black (stored in bl variable), InvertHardCopy, if paper = on (black background), if paper = off (white background)
% hold all %vital for using 'ColorOrder'
% 
% %sets axis parameters
% set(gca, 'Box', 'Off', 'FontSize', fontsize, 'FontName', fontname, 'LineWidth', 1, 'ColorOrder', colororder, 'color', bl, 'XColor', wh, 'YColor', wh, 'ZColor', wh)
% 
% %sets figure size, xSize and ySize are in centimeters
% xSize = 5.5; ySize = 6.3; xLeft = 0.5; yTop = 0.5;
% set(gcf,'PaperUnits','centimeters')
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
