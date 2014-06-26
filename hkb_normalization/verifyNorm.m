%VERIFYNORM  Makes plots to visually inspect quality of hkb normalization
%
%   [VALUES HKB MIDDLE] = NormalizeByHkb(dirname, pcelist)
%   VALUES is a set of line traces (one per pce).
%   HKB is a matrix containing: 
%       [geometric_mean ant_mean post_mean nCells_ant nCells_post]
%   MIDDLE is a matrix containing:
%       [middle_mean middle_95 middle_99 middle_nCells]
%
%
% Tara Martin, 2013-01-19


function [I_trimmed] = verifyNorm(dirname, pcestage, hkb, middle, ff)


%% check linear regression of hkb and lacZ expression

nTotal = length(pcestage); %number of embryos per genotype
I = 1:nTotal;

%calculate correlations
% corrTotal = corr(hkb(:,1)', middle(:,1)');
% corrThresholdHkb= corr(hkb_means_threshold(I)', middle_99(I)');
% corr95Hkb = corr(hkb_95(I)', middle_99(I)');

%% choose measurements to use for normalization
yy = middle(:,2); %middle columns are [mean 95% 99% nCells]
xx = hkb(:,1);  %hkb columns are [geomean ant_mean post_mean nCells]

%check for outlier embryos
stats = regstats(yy, xx, 'linear');
I_trimmed = I(stats.cookd<(4/(length(I)-2)));
 
%% remove outliers from calculations
yyTrimmed = yy(stats.cookd<(4/(length(I)-2)));
xxTrimmed = xx(stats.cookd<(4/(length(I)-2)));

nTrimmed = length(I_trimmed); %number of embryos kept
corrTrimmed = corr(xxTrimmed, yyTrimmed);
fprintf('%d embryos removed from set, %d kept.\nCorrelation = %1.2f. \n',...
    [nTotal-nTrimmed nTrimmed corrTrimmed]);

%% calculate linear regression
[aa, bb]=regress(yyTrimmed, [ones(size(I_trimmed))' xxTrimmed]);
slopeTrimmed = aa(2);
% slopeBoundTrimmed = aa(2) - bb(2, 1);
intersectTrimmed = aa(1);
% intersectBoundTrimmed = aa(1) - bb(1, 1);

fprintf('Slope = %1.2f.\n\n', slopeTrimmed);

stats2 = regstats(yyTrimmed, xxTrimmed, 'linear');
% FvalTrimmed = stats2.fstat.pval;

%% Plot trimmed data
if ~isempty(ff)
    figure(ff)
    %%
    subplot(311)
    hold on;
    scatter(xxTrimmed, yyTrimmed, 5, pcestage(I_trimmed)', 'filled');
    % plot(xxTrimmed, yyTrimmed, '.', 'MarkerSize', 10, 'Color', brewer10(mod(gg-1,10)+1, :))
    xlabel('geomean of hkb domains')
    ylabel('95% quantile of middle region')
    title('Linear regression')
%     caxis([0 25]);
%     colorbar;
    
    xlimits = get(gca, 'xlim');
%     ylimits = get(gca, 'ylim');
    plot(xlimits, xlimits.*slopeTrimmed+intersectTrimmed, 'LineWidth', 0.5, 'Color', 'k')
    
    %plot residuals
    subplot(312)
    hold on;
    scatter(xxTrimmed, stats2.standres, 5, pcestage(I_trimmed)', 'filled')
    xlimits = get(gca, 'xlim');
    plot(xlimits, [0 0], 'Color', 'k', 'LineWidth', 0.5);
    xlim(xlimits)
    ylim([-4 4])
    xlabel('mean of anterior hkb region')
    ylabel('standardized residual')
    title('Regression residuals')
%     caxis([0 25]);
%     colorbar;
    
    %plot normalized traces in third subplot
end
%     clear I xx yy xxTrimmed yyTrimmed

%% proof of principle analysis using mean normalization
%two different ways of normalizing with hkb
%values -- unnormlized line trace
%values_slope -- fancy normalized line trace (uses slope)
%values_hkb -- simple normalized line trace (just divide by hkb,
%already calculated previously)
% 
% for ii = I
%     %average of two methods: divide by hkb and using slope
%     temp = (hkb(ii)+middle_99(ii)/slopeTrimmed)/2;
%     values_slope(ii,:) = values(ii,:)./temp;
% end
