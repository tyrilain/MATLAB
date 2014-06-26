%%
% Calculating stripe ratios for fusion lines based on twi normalized
% controls.
%
% Important math for ratios (Z=X/Y):
% The mean ratio is just the ratio of two means... so far easy.
%   E(Z) = E(X)/E(Y)
% The standard deviation (SD) or standard error of the mean (SEM) of a 
% ratio is a lot harder to propagate.  I am assuming that I can treat them 
% the same way.
% This leads to the formula (assuming no covariance):
%  SEM(Z) = E(Z)*sqrt[ (SEM(X)/E(X))^2 + (SEM(Y)/E(Y))^2 ]

%% First get null hypothesis ratios from twi co-stained controls
% Run TwistNorm to generate data structures.
% TwistNorm;

% First calculate the SEM for all lines
peakSEM = [peakstd(:,:)'./[sqrt(numPC); sqrt(numPC)]]';

% Fusion A = lines 327+214
a = 4;  % line 327
b = 5;  % line 214
fusion_twi(1,:) = peaks(a,:)./peaks(b,:); % stripe3/4 and stripe7/6
fusionSEM_twi(1,:) = fusion_twi(1,:).* sqrt( (peakSEM(a,:)./peaks(a,:)).^2 +...
    (peakSEM(b,:)./peaks(b,:)).^2) ;

% Fusion B = lines 326+214
a = 3;  % line 326
b = 5;  % line 214
fusion_twi(2,:) = peaks(a,:)./peaks(b,:); % stripe3/4 and stripe7/6
fusionSEM_twi(2,:) = fusion_twi(2,:).* sqrt( (peakSEM(a,:)./peaks(a,:)).^2 +...
    (peakSEM(b,:)./peaks(b,:)).^2) ;

% Fusion C = lines 326+328
a = 3;  % line 326
b = 6;  % line 328
fusion_twi(3,:) = peaks(a,:)./peaks(b,:); % stripe3/4 and stripe7/6
fusionSEM_twi(3,:) = fusion_twi(3,:).* sqrt( (peakSEM(a,:)./peaks(a,:)).^2 +...
    (peakSEM(b,:)./peaks(b,:)).^2) ;

% Fusion D = lines 204+329
b = 1;  % line 204
a = 7;  % line 329
fusion_twi(4,:) = peaks(a,:)./peaks(b,:); % stripe3/4 and stripe7/6
fusionSEM_twi(4,:) = fusion_twi(4,:).* sqrt( (peakSEM(a,:)./peaks(a,:)).^2 +...
    (peakSEM(b,:)./peaks(b,:)).^2) ;

figure(2)
errorbar(1:4, fusion_twi(:,1), fusionSEM_twi(:,1), '.')
hold on;
errorbar(1:4, fusion_twi(:,2), fusionSEM_twi(:,2), '.r')
hold off;
title('SEM');

% for later convenience, re-assign numPC to new name
numPC_twi = numPC;


%% Now extract stripe ratios from actual fusion traces
% First run figure_fusions to generate data structures
% figure_fusions;

% Now extract peak positions for each stripe in each fusion, using control
% peak positions to estimate location.

% Get peaks for all lines, but only going to use controls
for i = 1:size(LacZ_mean,1);
    [temp ~] = peakdet(LacZ_mean(i,45:90), 0.1);
    temp(:,1) = temp(:,1)+44;
    LacZ_peak{i} = temp;
end

% Get index number for relevant lines
line204 = strcmp('LacZ:0204', channels);
line207 = strcmp('LacZ:0207', channels);
line208 = strcmp('LacZ:0208', channels);
line209 = strcmp('LacZ:0209', channels);
line210 = strcmp('LacZ:0210', channels);
line214 = strcmp('LacZ:0214', channels);
line326 = strcmp('LacZ:0326', channels);
line327 = strcmp('LacZ:0327', channels);
line328 = strcmp('LacZ:0328', channels);
line329 = strcmp('LacZ:0329', channels);

%%
w = 2; % window around control peak to search for fusion peak

% find Fusion A peaks
i = line209; % Fusion 1, line 209
a = line327;  % line 327
b = line214;  % line 214

% find peak within a certain window of control peak
peak(1,2) = max(LacZ_mean(i,(LacZ_peak{a}(1,1)-w:LacZ_peak{a}(1,1)+w)));
peak(2,2) = max(LacZ_mean(i,(LacZ_peak{b}(1,1)-1:LacZ_peak{b}(1,1)+w)));
peak(3,2) = max(LacZ_mean(i,(LacZ_peak{b}(2,1)-w:LacZ_peak{b}(2,1)+w)));
peak(4,2) = max(LacZ_mean(i,(LacZ_peak{a}(2,1)-w:LacZ_peak{a}(2,1)+w)));

% find index of peak
for j = 1:4
    peak(j,1) = find(LacZ_mean(i,:)==peak(j,2));
end

% update peak info
LacZ_peak{i} = peak;
LacZ_peakSEM{i} = LacZ_std(i,peak(:,1))./sqrt(LacZ_numPC(i));

% plot to check everything looks in order
figure(3);
plot(1:100, LacZ_mean(i,:), LacZ_peak{i}(:,1), LacZ_peak{i}(:,2), 'sr')
hold on;
errorbar(LacZ_peak{i}(:,1), LacZ_peak{i}(:,2), LacZ_peakSEM{i}, 'sr')
hold off;

fusion(1,:) = LacZ_peak{i}([1 4],2)./LacZ_peak{i}([2,3],2); % stripe3/4 and stripe7/6
fusionSEM(1,:) = fusion(1,:).* sqrt( (LacZ_peakSEM{i}([1 4])./LacZ_peak{i}([1 4],2)').^2 ...
    + (LacZ_peakSEM{i}([2 3])./LacZ_peak{i}([2 3],2)').^2) ;
numPC_LacZ(1) = LacZ_numPC(i);

%% find Fusion B peaks
i = line207; % Fusion 2, line 207
a = line326;  % line 326
b = line214;  % line 214

% find peak within a certain window of control peak
peak(1,2) = max(LacZ_mean(i,(LacZ_peak{a}(1,1)-w:LacZ_peak{a}(1,1)+w)));
peak(2,2) = max(LacZ_mean(i,(LacZ_peak{b}(1,1)-1:LacZ_peak{b}(1,1)+w)));
peak(3,2) = max(LacZ_mean(i,(LacZ_peak{b}(2,1)-2*w:LacZ_peak{b}(2,1))));
peak(4,2) = max(LacZ_mean(i,(LacZ_peak{a}(2,1)-2*w:LacZ_peak{a}(2,1))));

% find index of peak
for j = 1:4
    peak(j,1) = find(LacZ_mean(i,:)==peak(j,2));
end

% update peak info
LacZ_peak{i} = peak;
LacZ_peakSEM{i} = LacZ_std(i,peak(:,1))./sqrt(LacZ_numPC(i));

% plot to check everything looks in order
figure(3);
plot(1:100, LacZ_mean(i,:), LacZ_peak{i}(:,1), LacZ_peak{i}(:,2), 'sr')
hold on;
errorbar(LacZ_peak{i}(:,1), LacZ_peak{i}(:,2), LacZ_peakSEM{i}, 'sr')
hold off;

fusion(2,:) = LacZ_peak{i}([1 4],2)./LacZ_peak{i}([2,3],2); % stripe3/4 and stripe7/6
fusionSEM(2,:) = fusion(2,:).* sqrt( (LacZ_peakSEM{i}([1 4])./LacZ_peak{i}([1 4],2)').^2 ...
    + (LacZ_peakSEM{i}([2 3])./LacZ_peak{i}([2 3],2)').^2) ;
numPC_LacZ(2) = LacZ_numPC(i);

%% find Fusion C peaks
i = line210; % Fusion 3, line 210

% note for this line peakdet above finds the four peaks appropriately
LacZ_peakSEM{i} = LacZ_std(i,peak(:,1))./sqrt(LacZ_numPC(i));

% plot to check everything looks in order
figure(3);
plot(1:100, LacZ_mean(i,:), LacZ_peak{i}(:,1), LacZ_peak{i}(:,2), 'sr')
hold on;
errorbar(LacZ_peak{i}(:,1), LacZ_peak{i}(:,2), LacZ_peakSEM{i}, 'sr')
hold off;

fusion(3,:) = LacZ_peak{i}([1 4],2)./LacZ_peak{i}([2,3],2); % stripe3/4 and stripe7/6
fusionSEM(3,:) = fusion(3,:).* sqrt( (LacZ_peakSEM{i}([1 4])./LacZ_peak{i}([1 4],2)').^2 ...
    + (LacZ_peakSEM{i}([2 3])./LacZ_peak{i}([2 3],2)').^2) ;
numPC_LacZ(3) = LacZ_numPC(i);

%% find Fusion D peaks
i = line208; % Fusion 3, line 208
a = line204;  % line 204
b = line329;  % line 329

% find peak within a certain window of control peak
peak(1,2) = max(LacZ_mean(i,(LacZ_peak{a}(1,1)-w:LacZ_peak{a}(1,1)+w)));
peak(2,2) = max(LacZ_mean(i,(LacZ_peak{b}(1,1)-w:LacZ_peak{b}(1,1)+w)));
peak(3,2) = max(LacZ_mean(i,(LacZ_peak{b}(2,1)-w:LacZ_peak{b}(2,1)+w)));
peak(4,2) = max(LacZ_mean(i,(LacZ_peak{a}(2,1):LacZ_peak{a}(2,1)+w)));

% find index of peak
for j = 1:4
    peak(j,1) = find(LacZ_mean(i,:)==peak(j,2));
end

% update peak info
LacZ_peak{i} = peak;
LacZ_peakSEM{i} = LacZ_std(i,peak(:,1))./sqrt(LacZ_numPC(i));

% plot to check everything looks in order
figure(3);
plot(1:100, LacZ_mean(i,:), LacZ_peak{i}(:,1), LacZ_peak{i}(:,2), 'sr')
hold on;
errorbar(LacZ_peak{i}(:,1), LacZ_peak{i}(:,2), LacZ_peakSEM{i}, 'sr')
hold off;

fusion(4,:) = LacZ_peak{i}([2 3],2)./LacZ_peak{i}([1 4],2); % stripe3/4 and stripe7/6
fusionSEM(4,:) = fusion(4,:).* sqrt( (LacZ_peakSEM{i}([2 3])./LacZ_peak{i}([2 3],2)').^2 ...
    + (LacZ_peakSEM{i}([1 4])./LacZ_peak{i}([1 4],2)').^2) ;
numPC_LacZ(4) = LacZ_numPC(i);

%% Plot

% figure(3);
% errorbar(1:4, fusion_twi(:,1), 2*fusionSEM_twi(:,1), '.b')
% hold on;
% errorbar(1:4, fusion_twi(:,2), 2*fusionSEM_twi(:,2), '.c')
% errorbar(1:4, fusion(:,1), 2*fusionSEM(:,1), '.r')
% errorbar(1:4, fusion(:,2), 2*fusionSEM(:,2), '.m')
% hold off;
% legend('twi3/4', 'twi7/6', 'fusion3/4', 'fusion7/6')
% title('SEM');

figure(2);
xSize = 10;  ySize = 6;  % Set width and height of figure in cm

barwitherr([fusionSEM_twi fusionSEM],[fusion_twi fusion]);
set_fig;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel', ['Fusion A'; 'Fusion B'; 'Fusion C'; 'Fusion D']);
ylabel('Fold difference')
legend('twi3/4', 'twi7/6', 'fusion3/4', 'fusion7/6')
legend BOXOFF;
title('SEM');

% save figure
saveas(gcf, fullfile(save_dir, 'twi_fusions.eps'), 'epsc')

%% Tests of statistical significance
% Since we don't know if variances are equal, want to use a specific type
% of t-test (known as Welch t-test):
% t = (E(v1)-E(v2))/sqrt(s1^2/n1 + s2^2/n2)
%   where v1 and v2 are variables, s1 and s2 are respective standard
%   deviation, n1 and n2 are respective sample sizes
% v = (s1^2/n1 + s2^2/n2)^2 /( s1^4/(n1^2*(n1-1)) + s2^4/(n2^2*(n2-1)) )

% translating these to use the SEM instead of standard deviation gives:
% t = (E(v1)-E(v2))/sqrt(sem1^2 +sem2^2)
% v = (sem1^2 + sem2^2)^2 /( sem1^4/(n1-1) + sem2^4/(n2-1) )

% so I've got E(v1) = fusion_twi, sem1 = fusionSEM_twi and n1 = numPC_twi
%   E(v2) = fusion, sem2 = fusionSEM, and n2 = numPC

% for i = 1:4
%     t(i,:) = (fusion_twi(i,:) - fusion(i,:))./sqrt(fusionSEM_twi(i,:).^2 ...
%         + fusionSEM(i,:).^2);
% 
% end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    