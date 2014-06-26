dirname = '/Users/tmartin/Dropbox/data/twist/';
linenames = {'204', '325', '326', '327', '214', '328', '329', '330'};
G = {};
lineNum = [];
ii = 1;
xx = 40:90;
numPC = [];

% Run TwistData to generate data files
% TwistData(dirname, linenames);

% Load data for each line and extract relevant info
for d=1:length(linenames)
    load(fullfile(dirname, linenames{d}, [linenames{d} '.mat']));
    
    %%
    numPC(d) = length(cohort3.name); % want to keep track of sample size per line
    for i = 1:numPC(d)
        
        % use this code to check expression histograms for consistency
%         figure(1)
%         subplot(numPC,2,2*i-1)
%         hist(cohort3.twi{i},100)
%         title('twi')
%         subplot(numPC,2,2*i)
%         hist(cohort3.lac{i},100)
%         title('lacZ')
        
        % calculate mean twist
        twi.mean(ii) = mean(cohort3.twi{i});
%         twi.prct(ii,:) = prctile(cohort3.twi{i}, [05 50 95]);
        
        % normalized traces
        lac.trace(ii,:) = cohort3.lac{i}/twi.mean(ii);
        
        % find stripe peaks
        [temp ~] = peakdet(lac.trace(ii,xx), 0.1*max(lac.trace(ii,xx)));
        lac.peaksI(ii,1:2) = temp(:,1)+min(xx-1); % store peak indices
        lac.peaks(ii,1:2) = temp(:,2); % store peak values
        
        G{ii} = [linenames{d}]; % keeps track of line associated with index for making box plots
        lineNum(ii) = d;  % create bins for plotting individual embryo data

%         plot(1:100, lac.trace(ii,:), lac.peaksI(ii,:), lac.peaks(ii,:),'sr');
%         storeEmbryo = input('Hit enter to continue', 's');
        ii = ii+1;
    end

    I = (lineNum==d); % get indices for given line
    meantrace(d,:) = mean(lac.trace(I,:),1);
    peaks(d,:) = mean(lac.peaks(I,:),1);
    peakstd(d,:) = std(lac.peaks(I,:),1);
    peakI(d,:) = mean(lac.peaksI(I,:),1);

%     plot(1:100, meantrace(d,:), peakI(d,:), peaks(d,:),'sr')   
%     storeEmbryo = input('Hit enter to continue to next line...', 's');  
end

%% NOTES
% Plotting the twi mean vs. median or 95th percentile give pretty straight
% lines, with a nearly 1:1 ratio of mean to median, and 1:1.5 of mean to
% 95th percentile.  LacZ gives similar results, except with a 1:2 ratio of
% mean to 95th.  This says to me that these summary statistics are fairly
% well correlated, so it doesn't really matter which one I use (based on
% examination of line 204), but the ratio of lacZ to twi I get will vary 
% depending on summary statistic .  I should probably check that this 
% relationship holds for all of the lines.


%% Make figure of twist normalized peak heights for controls
% After running TwistNorm

figure(1)
xSize = 8;  ySize = 10;  % Set width and height of figure in cm

% stripes 3/7
temp = 1:4;
subplot(2,1,1);
barwitherr(peakstd(temp,:)./sqrt(repmat(numPC(temp)', 1,2)), peaks(temp,:));
set_fig;
% errorbar(temp, peaks(temp,1), peakstd(temp,1)'./sqrt(numPC(temp)), '.');
% hold on;
% errorbar(temp, peaks(temp,2), peakstd(temp,2)'./sqrt(numPC(temp)), '.r');
% hold off;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel', ['F-P'; 'R-P'; 'F-D'; 'R-D']);
title('Stripe 3/7 (standard error of means)');
legend('stripe 3', 'stripe 7');
legend BOXOFF;
ylim([0 3.5]);
xlabel('F=forward, R=reverse, P=proximal, D=distal');
ylabel('Peak fluorescence/mean twi');

% stripes 4/6
temp2 = 5:8;
subplot(2,1,2);
barwitherr(peakstd(temp2,:)./sqrt(repmat(numPC(temp2)', 1,2)), peaks(temp2,:));
set_fig;
% errorbar(temp2, peaks(temp2,1), peakstd(temp2,1)'./sqrt(numPC(temp2)), '.')
% hold on;
% errorbar(temp2, peaks(temp2,2), peakstd(temp2,2)'./sqrt(numPC(temp2)), '.r')
% hold off;
set(gca, 'XTick', temp, 'XTickLabel', ['F-P'; 'R-P'; 'F-D'; 'R-D']);
title('Stripe 4/6 (standard error of means)');
legend('stripe 4', 'stripe 6', 'Location','Northwest');
legend BOXOFF;
ylim([0 3.5]);
xlabel('F=forward, R=reverse, P=proximal, D=distal');
ylabel('Peak fluorescence/mean twi');

% Save
base_dir = '/Users/tmartin/'; % set manually because changes depending on computer
save_dir = [base_dir 'Dropbox/results/'];
saveas(gcf, fullfile(save_dir,'twi_controls.eps'), 'epsc')

%%
% figure(2)
% errorbar(1:length(linenames), peaks(:,1), peakstd(:,1), '.')
% hold on;
% errorbar(1:length(linenames), peaks(:,2), peakstd(:,2), '.r')
% hold off;
% title('standard deviation');
% 
% %%
% figure(3)
% temp = [[lac.mean./twi.mean]' lac.mean'./twi.prct(:,3) lac.prct(:,3)./twi.prct(:,3)];
% boxplot(temp, {'mean', 'mean/95th', '95th'}, 'notch', 'on');

