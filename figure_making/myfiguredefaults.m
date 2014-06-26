%MYFIGUREDEFAULTS Default settings for my figures, axes and lines.
%
% Tara Martin, 2014-03-10


function myfiguredefaults()

set(0,'DefaultFigureColor',[1 1 1]); % set background of figures to white
set(0,'DefaultFigureInvertHardCopy', 'off'); % force matlab to print figures as they look on screen
set(0,'DefaultFigurePaperUnits','centimeters'); % use cm to measure figure size


set(0,'DefaultAxesBox', 'off'); % do not put box around plot
set(0,'DefaultAxesColor', 'none'); % just use figure color
set(0,'DefaultAxesFontName', 'Helvetica'); % use Helvetica font
set(0,'DefaultAxesFontSize', 8); % use 10pt font
set(0,'DefaultAxesLineWidth', 0.5); % 1pt axis line width

set(0,'DefaultLineLineWidth', 1); % 1.5pt plot lines


% set(0,'DefaultFigurePaperPosition',[xLeft yTop xSize ySize])
