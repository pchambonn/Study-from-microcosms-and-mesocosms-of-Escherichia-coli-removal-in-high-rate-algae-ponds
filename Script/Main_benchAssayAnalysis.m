%% Main bench assay analysis
% This script is used to develop the statistical analysis of bench scale
% results from the publication "Study from microcosms and mesocosms reveals
% Escherichia coli removal in High Rate Alage Ponds is primarily caused by
% dark decay". This study is particularly presented in Supplementary
% Information S13.

optionLoadWorkspace = 1; 
if optionLoadWorkspace == 1
    load('../Datasets/bench_data_import_full_workspace.mat')
    optionBuildData = 0;
    optionImportData = 0;
end
optionFigure = 1;
optionStatAnalysis = 1;
optionPca = 1;
optionPcaSun = 1;
optionPcaDark = 1;

fs = 16;

if optionBuildData == 1
    % Here, data classifier are created for later hypothesis testing (also
    % used for data formatting)
    
    % Name of the experiment (for xls reading)
    expNames = {'0211 1A','0211 1B','0211 2A','0211 2B','0211 3A','0211 3B',...
        '1411 1A','1411 1B','1411 2A','1411 2B','1411 3A','1411 3B',...
        '1611 1A','1611 1B','1611 2A','1611 2B','1611 3A','1611 3B',...
        '2311 1A','2311 1B','2311 2A','2311 2B','2311 3A','2311 3B'};
    % Number of experiments included (i.e. reactors individually followed)
    nExp = 24;

    % Identifier of high pH (1) or low pH (0) reactors
    highPH = [ 1 , 1 , 1 , 1 , 1 , 1,...
        0 , 0 , 0 , 0 , 0 , 0 ,...
        1 , 0 , 1 , 0 , 1 , 0 ,...
        1 , 0 , 1 , 0 , 1 , 0];
    
    % Identifier of light attenuation coefficient
    sigmaType = [1 , 1 , 1 , 1 , 1 , 1 ,...
        2 , 2 , 2 , 2 , 2 , 2 ,...
        3 , 3 , 3 , 3 , 3 , 3 ,...
        4 , 4 , 4 , 4 , 4 , 4];
    
    % Identifier of light (1) or dark (0) reactors
    sunType = [1 , 1 , 1 , 1 , 0 , 0,...
        1 , 1 , 1 , 1 , 0 , 0,...
        1 , 1 , 1 , 1 , 0 , 0,...
        1 , 1 , 1 , 1 , 0 , 0];

    %% Data extraction
    % This paragraph serves the purpose of extracting the experimental to
    % matlab worskpace. Cell vector of data (one cell per reactors
    % followed)
    if optionImportData == 1
    
        timeDataFile = cell(nExp,1); 
        sunDataFile = cell(nExp,1); 
        pHDataFile = cell(nExp,1); 
        tempDataFile = cell(nExp,1);  
        DODataFile = cell(nExp,1);
        coliDataFile = cell(nExp,1);
        decayRateDataFile = cell(nExp,1);
        nDataFile = cell(nExp,1); 

        for i = 1:nExp
            A = xlsread('../../../Data sets/data_bench_assays.xlsx',expNames{i},'A:J');
            
            timeDataFile{i} = A(:,1);
            sunDataFile{i} = A(:,2);
            pHDataFile{i} = A(:,4);
            tempDataFile{i} = A(:,5);
            DODataFile{i} = A(:,6);
            coliDataFile{i} = A(:,8);
            decayRateDataFile{i} = A(:,9);
            nDataFile{i} = length(timeDataFile{i});
        end
    sigma0 = [55, 67, 66, 70]; % light attenuation coefficient for each bench experiment (m-1)
    d0 = 0.25; % bench reactor depth
    end
    
    %% Gathering of the E. coli and corresponding environmental parameters (averaged, min, max) measurements in single vectors 
    
    nData = 54; % Total number of data for E; coli cell count
    linearRegMeanVector = NaN(nData,6); % Vector which consist in averaged pH, DO, Temp, Sun, the decay rate and the log10 removal (LRV) calculated for each measurement interval.
    linearRegMaxVector = NaN(nData,4); % Vector which consist in max pH, DO, Temp, Sun.
    linearRegMinVector = NaN(nData,4); % Vector which consist in min pH, DO, Temp, Sun.

    iData = 1; % index of the data entered in the final tables.

    for iExp = 1:nExp

        indexColi = find(~isnan(coliDataFile{iExp}));
        if length(indexColi) >= 2
            for k_exp = 1:length(indexColi)-1
                linearRegMeanVector(iData + k_exp - 1,:) = [mean(pHDataFile{iExp}(indexColi(k_exp:k_exp + 1))),...
                    mean(DODataFile{iExp}(indexColi(k_exp:k_exp + 1))),...
                    mean(tempDataFile{iExp}(indexColi(k_exp:k_exp + 1))),...
                    mean(sunDataFile{iExp}(indexColi(k_exp:k_exp + 1))),...
                    log(coliDataFile{iExp}(indexColi(k_exp))/coliDataFile{iExp}(indexColi(k_exp+1)))/(timeDataFile{iExp}(indexColi(k_exp+1)) - timeDataFile{iExp}(indexColi(k_exp))),...
                    log10(coliDataFile{iExp}(indexColi(k_exp))/coliDataFile{iExp}(indexColi(k_exp+1)))];
    %                 decayrate_data_file{i_exp}(index_coli(k_exp + 1))];

                linearRegMaxVector(iData + k_exp - 1,:) = [max(pHDataFile{iExp}(indexColi(k_exp:k_exp + 1))),...
                    max(DODataFile{iExp}(indexColi(k_exp:k_exp + 1))),...
                    max(tempDataFile{iExp}(indexColi(k_exp:k_exp + 1))),...
                    max(sunDataFile{iExp}(indexColi(k_exp:k_exp + 1)))];

                linearRegMinVector(iData + k_exp - 1,:) = [min(pHDataFile{iExp}(indexColi(k_exp:k_exp + 1))),...
                    min(DODataFile{iExp}(indexColi(k_exp:k_exp + 1))),...
                    min(tempDataFile{iExp}(indexColi(k_exp:k_exp + 1))),...
                    min(sunDataFile{iExp}(indexColi(k_exp:k_exp + 1)))];    
            end
            iData = iData + length(indexColi)-1;
        end             
    end
end

%% Index of data subsamples creation

% Clustering of environmental parameters per class of magnitude.

categoryPH = zeros(nData,1);
categoryDO = zeros(nData,1);
categoryTemperature = zeros(nData,1);
categorySun = zeros(nData,1);

for i = 1:nData
    if linearRegMeanVector(i,1) < 9.3
        categoryPH(i) = 0;
    else
        categoryPH(i) = 1;
    end

    if linearRegMeanVector(i,2) < 8
        categoryDO(i) = 0;
    else
        categoryDO(i) = 1;
    end

    if linearRegMeanVector(i,3) <= median(linearRegMeanVector(:,3))
        categoryTemperature(i) = 0;
    else
        categoryTemperature(i) = 1;
    end

    if linearRegMeanVector(i,4) == 0
        categorySun(i) = 0;
    else
        categorySun(i) = 1;
    end
end

indexDark = find(linearRegMeanVector(:,4) == 0);
indexSun = find(linearRegMeanVector(:,4) > 0);

indexHighPH = find(categoryPH == 1);
indexLowPH = find(categoryPH == 0);

indexHighDO = find(categoryDO == 1);
indexLowDO = find(categoryDO == 0);

%% Plotting SI figures
if optionFigure == 1
    % Figure S13-1: Decay rates during bench in darkness
    figure(1), clf, hold on
    plot(linearRegMeanVector(~categorySun & categoryDO,1),linearRegMeanVector(~categorySun & categoryDO,5),'ok','MarkerFaceColor','k')
    plot(linearRegMeanVector(~categorySun & ~categoryDO,1),linearRegMeanVector(~categorySun & ~categoryDO,5),'+k')
    txt = text(linearRegMeanVector(~categorySun,1) - 0.1,linearRegMeanVector(~categorySun,5),num2str(linearRegMeanVector(~categorySun,3),3),'HorizontalAlignment','right','FontSize',fs-2);
    ax = gca;
    ax.FontSize = fs-2;
    xlabel('pH','FontSize',fs,'FontName','Arial')
    ylabel('{\itE. coli} decay rate (d^{-1})','FontSize',fs,'FontName','Arial')
    xlim([5 11])
    ylim([0 200])
    ax.XTick = 5:2:11;

    % DO FINAL EDITIING OF THIS FIGURE AND PRINT IT OUT: GOOD TO GO OTHERWISE

    % Figure S13-2: E. coli decay rates under sunlight measured during
    % bench assays grouped for pH < 8.4 and pH > 9.3
    figure(2), clf, hold on
    plot(1*ones(sum(categorySun & ~categoryPH),1),linearRegMeanVector(categorySun & ~categoryPH,5),'+k')
    plot(2*ones(sum(categorySun & categoryPH),1),linearRegMeanVector(categorySun & categoryPH,5),'ok','MarkerFaceColor','k')
    xlim([0.5 2.5])
    ax = gca;
    ax.FontSize = fs-2;
    ylabel('{\itE. coli} decay rate (d^{-1})','FontSize',fs,'FontName','Arial')
    ax.XTick = [1,2];
    ax.XTickLabel = {'pH < 8.4' , 'pH > 9.3'}
    ax.FontName = 'Arial'

    % Figure S13-3:  E. coli decay rates under sunlight measured during
    % bench assays grouped for DO < 2 mg/L and DO > 10 g/L
    figure(3), clf, hold on
    plot(1*ones(sum(categorySun & ~categoryDO),1),linearRegMeanVector(categorySun & ~categoryDO,5),'+k')
    plot(2*ones(sum(categorySun & categoryDO),1),linearRegMeanVector(categorySun & categoryDO,5),'ok')
    xlim([0.5 2.5])
    ax = gca;
    ax.FontSize = fs-2;
    ylabel('{\itE. coli} decay rate (d^{-1})','FontSize',fs,'FontName','Arial')
    ax.XTick = [1,2];
    ax.XTickLabel = {'DO < 2 mg{\cdot}L^{-1}' , 'DO > 10 mg{\cdot}L^{-1}'}
    ax.FontName = 'Arial'

end
%% STATISTICAL ANALYSIS

if optionStatAnalysis == 1
    
%% Tests for homogeneity of parameters on compared subsets
    
% In the following, parameters accross different tested subpopulation are
% tested for homogeneity based on Kruskal-Wallis hypothesis testing.
    % The naming below goes as, in camel writing, pKw'parameter impact for hypothesis
    % testing'_'subdataset tested'_'parameter compared across cluster'
    % For instance: 
        % pKwPH_All_Sun is the p-value testing the hypothesis that sunlight
        % has the same distribution for pH clusters when testing on the
        % whole dataset.
    
        
    % Tests on whole data set
    
    [pKwPH_All_Sun, tblKwPH_All_Sun, statsKwPH_All_Sun] = kruskalwallis([linearRegMeanVector(categoryPH == 1,4);...
        linearRegMeanVector(categoryPH == 0,4)],[zeros(sum(categoryPH == 1),1);ones(sum(categoryPH == 0),1)]);
    [pKwPH_All_DO, tblKwPH_All_DO, statsKwPH_All_DO] = kruskalwallis([linearRegMeanVector(categoryPH == 1,2);...
        linearRegMeanVector(categoryPH == 0,2)],[zeros(sum(categoryPH == 1),1);ones(sum(categoryPH == 0),1)]);
    [pKwPH_All_Temperature, tblKwPH_All_Temperature, statsKwPH_All_Temperature] = kruskalwallis([linearRegMeanVector(categoryPH == 1,3);...
        linearRegMeanVector(categoryPH == 0,3)],[zeros(sum(categoryPH == 1),1);ones(sum(categoryPH == 0),1)]);

    [pKwDO_All_Sun, tblKwDO_All_Sun, statsKwDO_All_Sun] = kruskalwallis([linearRegMeanVector(categoryDO == 1,4);...
        linearRegMeanVector(categoryDO == 0,4)],[zeros(sum(categoryDO == 1),1);ones(sum(categoryDO == 0),1)]);
    [pKwDO_All_PH, tblKwDO_All_PH, statsKwDO_All_PH] = kruskalwallis([linearRegMeanVector(categoryDO == 1,1);...
        linearRegMeanVector(categoryDO == 0,1)],[zeros(sum(categoryDO == 1),1);ones(sum(categoryDO == 0),1)]);
    [pKwDO_All_Temperature, tblKwDO_All_Temperature, statsKwDO_All_Temperature] = kruskalwallis([linearRegMeanVector(categoryDO == 1,3);...
        linearRegMeanVector(categoryDO == 0,3)],[zeros(sum(categoryDO == 1),1);ones(sum(categoryDO == 0),1)]);

    [pKwSun_All_DO, tblKwSun_All_DO, statsKwSun_All_DO] = kruskalwallis([linearRegMeanVector(categorySun == 1,2);...
        linearRegMeanVector(categorySun == 0,2)],[zeros(sum(categorySun == 1),1);ones(sum(categorySun == 0),1)]);
    [pKwSun_All_PH, tblKwSun_All_PH, statsKwSun_All_PH] = kruskalwallis([linearRegMeanVector(categorySun == 1,1);...
        linearRegMeanVector(categorySun == 0,1)],[zeros(sum(categorySun == 1),1);ones(sum(categorySun == 0),1)]);
    [pKwSun_All_Temperature, tblKwSun_All_Temperature, statsKwSun_All_Temperature] = kruskalwallis([linearRegMeanVector(categorySun == 1,3);...
        linearRegMeanVector(categorySun == 0,3)],[zeros(sum(categorySun == 1),1);ones(sum(categorySun == 0),1)]);
    
    % Tests on sunlight only data set
    [pKwPH_Sun_Sun, tblKwPH_Sun_Sun, statsKwPH_Sun_Sun] = kruskalwallis([linearRegMeanVector(categorySun == 1 & categoryPH == 1,4);...
        linearRegMeanVector(categorySun == 1 & categoryPH == 0,4)],[zeros(sum(categorySun == 1 & categoryPH == 1),1);ones(sum(categorySun == 1 & categoryPH == 0),1)]);
    [pKwPH_Sun_DO, tblKwPH_Sun_DO, statsKwPH_Sun_DO] = kruskalwallis([linearRegMeanVector(categorySun == 1 & categoryPH == 1,2);...
        linearRegMeanVector(categorySun == 1 & categoryPH == 0,2)],[zeros(sum(categorySun == 1 & categoryPH == 1),1);ones(sum(categorySun == 1 & categoryPH == 0),1)]);
    [pKwPH_Sun_Temperature, tblKwPH_Sun_Temperature, statsKwSun_Sun_Temperature] = kruskalwallis([linearRegMeanVector(categorySun == 1 & categoryPH == 1,3);...
        linearRegMeanVector(categorySun == 1 & categoryPH == 0,3)],[zeros(sum(categorySun == 1 & categoryPH == 1),1);ones(sum(categorySun == 1 & categoryPH == 0),1)]);

    [pKwDO_Sun_Sun, tblKwDO_Sun_Sun, statsKwDO_Sun_Sun] = kruskalwallis([linearRegMeanVector(categorySun == 1 & categoryDO == 1,4);...
        linearRegMeanVector(categorySun == 1 & categoryDO == 0,4)],[zeros(sum(categorySun == 1 & categoryDO == 1),1);ones(sum(categorySun == 1 & categoryDO == 0),1)]);
    [pKwDO_Sun_PH, tblKwDO_Sun_PH, statsKwDO_Sun_PH] = kruskalwallis([linearRegMeanVector(categorySun == 1 & categoryDO == 1,1);...
        linearRegMeanVector(categorySun == 1 & categoryDO == 0,1)],[zeros(sum(categorySun == 1 & categoryDO == 1),1);ones(sum(categorySun == 1 & categoryDO == 0),1)]);
    [pKwDO_Sun_Temperature, tblKwDO_Sun_Temperature, statsKwPH_Sun_Temperature] = kruskalwallis([linearRegMeanVector(categorySun == 1 & categoryDO == 1,3);...
        linearRegMeanVector(categorySun == 1 & categoryDO == 0,3)],[zeros(sum(categorySun == 1 & categoryDO == 1),1);ones(sum(categorySun == 1 & categoryDO == 0),1)]);
        
    
    %% ttest2 testing hypothesis of higher decay rate in groups: hypothesis of unequal variance equivalent to Welch anova
    
    % 2 samples ttest testing if the data grouped by (tested) parameters
    % categories (e.g. high pH vs low pH) have the same mean assuming
    % unequal variance.
    % The test are performed as one tailed t-test assuming:
        % high pH has higher impact than low pH (i.e. decay rate is higher
            % in the subgroup of high pH than in the subgroup of low pH).
        % high DO has higher impact than low DO
        % sunlight has higher impact than dark
        
    % On all data:
    [hPH_All,pPH_All] = ttest2(linearRegMeanVector(categoryPH == 1,5),linearRegMeanVector(categoryPH == 0,5),'tail','right','Vartype','unequal');
    [hDO_All,pDO_All] = ttest2(linearRegMeanVector(categoryDO == 1,5),linearRegMeanVector(categoryDO == 0,5),'tail','right','Vartype','unequal');
    [hSun_All,pSun_All] = ttest2(linearRegMeanVector(categorySun == 1,5),linearRegMeanVector(categorySun == 0,5),'tail','right','Vartype','unequal');

    % On sun data only
    [hPH_Sun,pPH_Sun] = ttest2(linearRegMeanVector(categoryPH == 1 & categorySun == 1,5),...
        linearRegMeanVector(categoryPH == 0 & categorySun == 1,5),'tail','right','Vartype','unequal');
    [hDO_Sun,pDO_Sun] = ttest2(linearRegMeanVector(categoryDO == 1 & categorySun == 1,5),...
        linearRegMeanVector(categoryDO == 0 & categorySun == 1,5),'tail','right','Vartype','unequal');

    %On dark data only
    [hPH_Dark,pPH_Dark] = ttest2(linearRegMeanVector(categoryPH == 1 & categorySun == 0,5),...
        linearRegMeanVector(categoryPH == 0 & categorySun == 0,5),'tail','right','Vartype','unequal');
    [hDO_Dark,pDO_Dark] = ttest2(linearRegMeanVector(categoryDO == 1 & categorySun == 0,5),...
        linearRegMeanVector(categoryDO == 0 & categorySun == 0,5),'tail','right','Vartype','unequal');
    
end


