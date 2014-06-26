% Author: Tara Martin
% Date: 20110415
%
% Purpose: extract input and output data from atlas.vpc, perform linear
% regression on data, store parameters, calculate predicted expression
% based on the regression parameters and store for comparison to actual
% expression.
%


%% Import atlas and initialize parameters

%uncomment to import atlas data:
% atlas =  readpointcloud('atlas.vpc');
% dmel_atlas = readpointcloud('D_mel_wt__atlas__08_02_26.vpc');

% PCs to extract
regList = {'cad' 'gt' 'kni' 'Kr' 'slp1'};  %mRNA inputs to regression
protRegList = {'bcdP' 'hbP'};  %protein inputs to regression
coeff_key = ['const' regList protRegList];
%geneList = {'LacZ_0204' 'eve'};  %outputs to fit
geneList = {'eve' 'LacZ_0204' 'LacZ_0207' 'LacZ_0208' 'LacZ_0209' 'LacZ_0210' 'LacZ_0211' 'LacZ_0214'};

cohort = 5:9;  %timepoints to include, normally 4:9
numNuc = 6078;  %number of nuclei in atlas
regulators = zeros(numNuc*length(cohort), length(regList));  %initialize data structure


%% Extract PCs of regulators from atlases

% first get regulators with mRNA patterns
for j=1:length(regList)
    for k=1:length(cohort)
        fieldString = [regList{j} '__' num2str(cohort(k))];
        if (~isfield(atlas, fieldString))  %check if PC is in atlas
            warning('pcRegression:ExtractData', 'Atlas does not contain %s', fieldString);
        else
            eval(['reg = atlas.' fieldString ';']);  %creates data holder 'reg'
            regulators(((k-1)*numNuc+1):(k*numNuc),j) = reg; %concatenates all cohorts
        end
    end
    warning('pcRegression:testing', 'Extracted %s from atlas.', regList{j});
end

% now get regulators with protein patterns
for j=1:length(protRegList)
    for k=1:length(cohort)
        fieldString = [protRegList{j} '__' num2str(cohort(k))];
        if (~isfield(dmel_atlas, fieldString))  %check if PC is in atlas
            warning('pcRegression:ExtractData', 'Dmel atlas does not contain %s', fieldString);
        else
            eval(['reg = dmel_atlas.' fieldString ';']);  %creates data holder 'reg'
            regulators(((k-1)*numNuc+1):(k*numNuc),(j+length(regList))) = reg; %concatenates all cohorts
        end
    end
    warning('pcRegression:testing', 'Extracted %s from dmel atlas.', protRegList{j});
end


%% Extract vector of output values from atlas, fit each

for j=1:length(geneList)
    pc =zeros(numNuc*length(cohort),1); %creates data holder 'pc'
    for k=1:length(cohort)
        
        fieldString = [geneList{j} '__' num2str(cohort(k))];
        if (~isfield(atlas, fieldString))
            warning('pcRegression:ExtractData', 'Atlas does not contain %s', fieldString);
        else
            eval(['pc(((k-1)*numNuc+1):(k*numNuc)) = atlas.' fieldString ';']); %stores cohort data in 'pc'
        end
    end
    
    %threshold output data -- should do this in more coherent way ->
    %mode+stdev
    pc = gt(pc, 0.35);
    
%     coeffs = glmfit(regulators, pc, 'binomial');  %logistic regression
%     reg_fit = glmval(coeffs, regulators, 'logit'); %predicted expression given fit parameters
    
    coeffs = glmfit(regulators, pc);  %linear regression
    reg_fit = glmval(coeffs, regulators, 'identity'); %predicted expression given fit parameters

    %store PCs, coefficients of regression and predicted expression (fit)
    eval(['pc_' geneList{j} '=pc;']);
    eval(['coeffs_' geneList{j} ' = coeffs;']);
    eval(['fit_' geneList{j} '= reg_fit;']);
    
end


%% Visualize results

%see pc_plots.m



