%% Main lab assay analysis

% This script aims at producing all figures and statistical analysis
% obtained for lab assay experiments for the publication "Study from
% microcosms and mesocosms reveals Escherichia coli removal in High Rate
% Alage Ponds is primarily caused by dark decay" 

fs = 16; % Base fontsize used for figures
option_import = 0; % Option to import (=1) or not (=0) the laboratory experiments results when initializing the analysis

% List of option to target the generation of results on a per mechanism
% basis.
option_starvation = 0;
option_toxicity = 0;
option_pH = 1;
option_ammonia = 0;
option_sun_direct_decay = 0;
option_sun_pH = 0;
option_exogenousPO = 0;
option_exogenousPO_pH10 = 0;

%% Structure of the file
% The data from lab experiment is first loaded and stored into tables (both
% numeric and text).
% A section is created on a per mechanism invetigated basis.
% In each section, a boolean list is created which enables to select all
% the lines of interest. Having selected such data, the numerical
% statistical and figures needed are created.

%% Import results
if option_import == 1
    [table_results,table_results_txt] = xlsread('../../../Data sets/data_lab assays.xlsx');
    table_results_txt = cell2table(table_results_txt);
    table_results_txt = table_results_txt(2:end,:);
end

%% Starvation and inactivation

if option_starvation == 1
    % Creation of boolean list for reactors targeting starvation.
    list_dark = table_results(:,1) == 0;
    list_RO_water = strcmp(table_results_txt{:,4} , 'RO water');
    list_addition = strcmp(table_results_txt{:,5} , '') | strcmp(table_results_txt{:,5} , 'pH buffer');
    list_neutral_pH = not(table_results(:,6) > 7.5 | table_results(:,6) < 6.5 | table_results(:,5) > 7.5 | table_results(:,5) < 6.5);

    list_starvation = list_dark & list_RO_water & list_neutral_pH & list_addition;

    % Figure S7-1: Distribution of first order decay rate calculated from
    % laboratory experiments performed in dark microcosms filled with RO
    % water or neutral pH buffer
    figure(1), clf, hold on
    histogram(table_results(list_starvation,12), 25,'FaceColor',0.25*[1,1,1])
    axis = gca;
    axis.FontSize = fs - 3;
    xlabel('First order decay rate (d^{-1})','FontSize',fs)
    ylabel('Frequency','FontSize',fs)
    xlim([-10 10])

    %Figrue S7-2: E. coli decay rates in dark controls during temperature
    %controlled experiments according to the incubation temperature.
    
    % Addition of jitter: a jitter was added to separate visually the
    % different data points
    l = length(table_results(list_starvation,12));
    jitter = -0.5 + 1.5*rand(l,1);
    
    figure(2), clf, hold on
    scatter(table_results(list_starvation,9) + jitter,table_results(list_starvation,12),'dk','MarkerFaceColor' ,'k')
    xlim([0 40])
    x = linspace(0,40);
    y = linspace(0,0);
    plot(x,y,'-k')
    errorbar(table_results(list_starvation,9)+ jitter,table_results(list_starvation,12),table_results(list_starvation,15),'LineStyle','none','Color','k','LineWidth',1)
    axis = gca;
    axis.FontSize = fs - 3;
    xlabel('Temperature (°C)','FontSize',fs)
    ylabel('First order decay rate (d^{-1})','FontSize',fs)
    
    % Statistical analysis
    mean(table_results(list_starvation,12))
    std(table_results(list_starvation,12))
    numel(table_results(list_starvation,12))
    [h_starvation,p_starvation,ci_starvation,stats_starvation] = ttest(table_results(list_starvation,12),0,'Alpha',0.05,'Tail','both')

    % Visual check of experiments included duration

    table_duration_check = [table_results_txt(list_starvation,13),table(table_results(list_starvation,12))];
end
%% Toxicity in filtrates

if option_toxicity == 1
    % Reduction to data of interest
    list_dark = table_results(:,1) == 0;
    list_matrix = strcmp(table_results_txt{:,4} , 'Filtrated wastewater') | strcmp(table_results_txt{:,4} , 'Filtrated wastewater, autoclaved')...
        | strcmp(table_results_txt{:,4} , 'Filtrated HRAP broth');
    list_addition = strcmp(table_results_txt{:,5} , '') | strcmp(table_results_txt{:,5} , 'pH buffer');
    list_neutral_pH = not(table_results(:,6) > 8.5 | table_results(:,6) < 6.5 | table_results(:,5) > 8.5 | table_results(:,5) < 6.5);

    list_filtrates = list_dark & list_matrix & list_neutral_pH & list_addition;

    % Visual check of the data
    figure(1), clf,
    plot(table_results(list_filtrates,12),'x')
end

%% pH toxicity

if option_pH == 1
    % Reduction to the data of interest
    list_dark = strcmp(table_results_txt{:,2} , 'Indoor');
    list_matrix = strcmp(table_results_txt{:,4} , 'RO water');
    list_addition = strcmp(table_results_txt{:,5} , '') | strcmp(table_results_txt{:,5} , 'pH buffer');
    list_temperature = strcmp(table_results_txt{:,10} , 'Y');

    list_pH = list_dark & list_matrix & list_addition & list_temperature;

    list_pH_below_10_5 = table_results(:,7) < 10.6;
    list_pH_final = list_pH & list_pH_below_10_5;

    % Figure S9-1: INfluence of pH on E. coli decay rate at 30°C in
    % laboratory microcosms
    list_T_30 = table_results(:,9) == 30;
    list_pH_T_30 = list_pH_final & list_T_30;

    figure(1), clf, hold on
    plot(table_results(list_pH_T_30,7),table_results(list_pH_T_30,12),'dk','MarkerFaceColor' ,'k')
    % Error bars
    err = table_results(list_pH_T_30,15);
    errorbar(table_results(list_pH_T_30,7),table_results(list_pH_T_30,12),err,'LineStyle','None','Color','k','LineWidth',1)
    % Fit
    x = linspace(6,11);
    x0 = table_results(list_pH_T_30,7);
    y0 = table_results(list_pH_T_30,12);
    fun = @(a,x)(a*10.^(x-14));
    [fit_30, gof_30,~] = fit(x0,y0,fun);
    y = fit_30.a*10.^(x-14);
    plot(x,y,'-k')
    xlim([6 11])
    ylim([-25 350])
    axis = gca;
    axis.FontSize = fs - 3;
    xlabel('pH','FontSize',fs)
    ylabel('First order decay rate (d^{-1})','FontSize',fs)
    
    % In the following,values for a(T) at different temperatures are
    % obtained by fitting the pH dependency law for a given temperature.
    % Plots are also produced for visual checks
    figure(2), clf, hold on
        % 5°C
    list_T_5 = table_results(:,9) == 5;
    list_pH_T_5 = list_pH_final & list_T_5;
    x0 = table_results(list_pH_T_5,7);
    y0 = table_results(list_pH_T_5,12);
    fun = @(a,x)(a*10.^(x-14));
    [fit_5, gof_5,~] = fit(x0,y0,fun)
    subplot(7,1,1)
    plot(fit_5,x0,y0)

        % 10°C
    list_T_10 = table_results(:,9) == 10;
    list_pH_T_10 = list_pH_final & list_T_10;
    x0 = table_results(list_pH_T_10,7);
    y0 = table_results(list_pH_T_10,12);
    fun = @(a,x)(a*10.^(x-14));
    [fit_10, gof_10,~] = fit(x0,y0,fun)
    subplot(7,1,2)
    plot(fit_10,x0,y0)

        % 15°C
    list_T_15 = table_results(:,9) == 15;
    list_pH_T_15 = list_pH_final & list_T_15;
    x0 = table_results(list_pH_T_15,7);
    y0 = table_results(list_pH_T_15,12);
    fun = @(a,x)(a*10.^(x-14));
    [fit_15, gof_15,~] = fit(x0,y0,fun) 
    subplot(7,1,3)
    plot(fit_15,x0,y0)

        % 20°C
    list_T_20 = table_results(:,9) == 20;
    list_pH_T_20 = list_pH_final & list_T_20;
    x0 = table_results(list_pH_T_20,7);
    y0 = table_results(list_pH_T_20,12);
    fun = @(a,x)(a*10.^(x-14));
    [fit_20, gof_20,~] = fit(x0,y0,fun)
    subplot(7,1,4)
    plot(fit_20,x0,y0)

        % 25°C
    list_T_25 = table_results(:,9) == 25;
    list_pH_T_25 = list_pH_final & list_T_25;
    x0 = table_results(list_pH_T_25,7);
    y0 = table_results(list_pH_T_25,12);
    fun = @(a,x)(a*10.^(x-14));
    [fit_25, gof_25,~] = fit(x0,y0,fun)
    subplot(7,1,5)
    plot(fit_25,x0,y0)

        % 30°C
    list_T_30 = table_results(:,9) == 30;
    list_pH_T_30 = list_pH_final & list_T_30;
    x0 = table_results(list_pH_T_30,7);
    y0 = table_results(list_pH_T_30,12);
    fun = @(a,x)(a*10.^(x-14));
    [fit_30, gof_30,~] = fit(x0,y0,fun)
    subplot(7,1,6)
    plot(fit_30,x0,y0)

        % 35°C
    list_T_35 = table_results(:,9) == 35;
    list_pH_T_35 = list_pH_final & list_T_35;
    x0 = table_results(list_pH_T_35,7);
    y0 = table_results(list_pH_T_35,12);
    fun = @(a,x)(a*10.^(x-14));
    [fit_35, gof_35,~] = fit(x0,y0,fun)
    subplot(7,1,7)
    plot(fit_35,x0,y0)
    
    % Figure S9-2: Influence of temperature on E. coli decay rate at pH 10
    % in laboratory microcosms

    list_pH_10 = table_results(:,7) > 9.8 & table_results(:,7) < 10.2;
    list_pH_pH_10 = list_pH_final & list_pH_10;
    
    l = length(table_results(list_pH_pH_10,9));
    jitter = -0.5 + rand(l,1);
    figure(3), clf, hold on
    plot(table_results(list_pH_pH_10,9) + jitter ,table_results(list_pH_pH_10,12),'dk','MarkerFaceColor' ,'k')
    % Error bars
    err = table_results(list_pH_pH_10,15);
    errorbar(table_results(list_pH_pH_10,9) + jitter,table_results(list_pH_pH_10,12),err,'LineStyle','None','Color','k','LineWidth',1)
    % Fitting
    x = linspace(0,40);
    x0 = table_results(list_pH_pH_10,9);
    y0 = table_results(list_pH_pH_10,12);
    fun = @(a,b,x)(a*exp(b*x));
    [fit_pH10, gof_pH10,~] = fit(x0,y0,fun,'StartPoint',[1,0.1]);
    y = fit_pH10.a*exp(fit_pH10.b*x);
    plot(x,y,'-k')
    xlim([0 40])
    ylim([0 200])
    axis = gca;
    axis.FontSize = fs - 3;
    xlabel('Temperature (°C)','FontSize',fs)
    ylabel('First order decay rate (d^{-1})','FontSize',fs)

    % Figure S9-3: Variation of the log transformed E. coli decay at pH 14
    % (a(T)) according to temperature during laboratory assays
    T = [5, 10, 20 , 25 , 30 , 35];
    a = [fit_5.a , fit_10.a , fit_20.a , fit_25.a , fit_30.a , fit_35.a ];
    ci = [confint(fit_5) , confint(fit_10) , confint(fit_20) , confint(fit_25) , confint(fit_30) , confint(fit_35)];
    err_neg_a = a - ci(1,:);
    err_pos_a = ci(2,:) - a;
    err_neg = log(a) - log(a - err_neg_a);
    err_pos = log(a + err_neg_a) - log(a);

    figure(4), clf, hold on
    plot(T,log(a),'ok','MarkerSize',10,'LineWidth',2)
    xlim([0 40])
    ylim([0 15])
    axis = gca;
    axis.FontSize = fs - 3;
    xlabel('Temperature (°C)','FontSize',fs)
    ylabel('ln(a)','FontSize',fs)
    mdl = fitlm(T,log(a));
    x = linspace(0,40);
    y = mdl.Coefficients{1,1} + mdl.Coefficients{2,1}*x;
    plot(x,y,'-k')
    errorbar(T,log(a),err_neg,err_pos,'LineStyle','None','Color','k','LineWidth',1)
    txt = {char(strcat('y =',{' '},num2str(mdl.Coefficients{2,1},3),'\cdot','x + ',{' '},num2str(mdl.Coefficients{1,1},3))),...
            char(strcat('R^2 =',{' '},num2str(mdl.Rsquared.Ordinary,3)))};
    text(1,7,txt,'FontSize',fs - 2)
    
    % Statistical analysis
    % pH 10
    list_pH_10 = table_results(:,5) == 10;
    list_T_30 = table_results(:,9) == 30;
    list_pH_10_T_30 = list_pH_10 & list_T_30;
    mean(table_results(list_pH_10_T_30,12));
    std(table_results(list_pH_10_T_30,12));
    
    % Model fit:
    teta = exp(mdl.Coefficients{2,1});
    k_pH_20 = exp(mdl.Coefficients{1,1} + 20*mdl.Coefficients{2,1});

    % Figure 1: Measured vs fitted E. coli decay rates for all tested pH
    % and temperature from laboratory assays
    figure(5), clf, hold on
    plot(k_pH_20*teta.^(table_results(list_pH_final,9)-20).*(10.^(table_results(list_pH_final,7)-14)),table_results(list_pH_final,12),'dk','MarkerFaceColor' ,'k')
    err = table_results(list_pH_final,15);
    errorbar(k_pH_20*teta.^(table_results(list_pH_final,9)-20).*(10.^(table_results(list_pH_final,7)-14)),table_results(list_pH_final,12),err,'LineStyle','None','Color','k','LineWidth',1)
    x = linspace(-10,300);
    y = x;
    plot(x,y,'-k')
    xlim([-10 300])
    ylim([-10 300])
    axis = gca;
    axis.FontSize = fs - 3;
    xlabel('Predicted {\itE. coli} decay rate (d^{-1})','FontSize',fs,'FontName','Arial','FontWeight','bold')
    ylabel('Measured {\itE. coli} decay rate (d^{-1})','FontSize',fs,'FontName','Arial','FontWeight','bold')

    SS_tot = sum((table_results(list_pH_final,12)-mean(table_results(list_pH_final,12))).^2);
    SS_res = sum((table_results(list_pH_final,12) - k_pH_20*teta.^(table_results(list_pH_final,9)-20).*(10.^(table_results(list_pH_final,7)-14))).^2);
    det_coeff = 1 - SS_res/SS_tot;
    n_data = length(table_results(list_pH_final,12));
    SE = sqrt(1/(n_data - 2)*SS_res/(SS_tot));


    % This was added to have a look at the lab data at very high pH: 
    list_pH_above_10_5 = table_results(:,7) >= 10.6;
    list_pH_high = list_pH & list_pH_above_10_5;
    table_results(list_pH_high,:)
end



%% Ammonia toxicity

if option_ammonia == 1
    % Reduction to data of interest:
    date_1 = '18/01/2017';
    date_2 = '19/01/2017';
    date_3 = '25/01/2017';
    date_4 = '11/08/2017';
    date_5 = '25/08/2017';
    % Experiment dates are used to select data

    list_dark = table_results(:,1) == 0;
    list_matrix = strcmp(table_results_txt{:,4} , 'RO water');
    list_addition = strcmp(table_results_txt{:,1} ,date_1) | strcmp(table_results_txt{:,1} ,date_2) | strcmp(table_results_txt{:,1} ,date_3) |...
        strcmp(table_results_txt{:,1} ,date_4) ;

    list_ammonia_test = list_dark & list_matrix & list_addition;

    % Results are categorized into no ammonia, HRAP level, and wastewater
    % level of ammonia:
        % This is done by classifying ammonia test levels in three levels
        % 1,2,3 by order of magnitude
    % A table is created to gather all the needed info
    table_NH3_experiment = table_results(list_ammonia_test,:);
    table_NH3_experiment_txt = table_results_txt(list_ammonia_test,:);
    n =length(table_NH3_experiment_txt{:,5});
    NH3_levels = NaN(n,1);
    for i = 1:n
        if contains(table_NH3_experiment_txt{i,5},'1 mg/L N-NH3') | contains(table_NH3_experiment_txt{i,5},'0.5 mg N-NH3/L')
            NH3_levels(i) = 2;
        else
            if contains(table_NH3_experiment_txt{i,5},'30 mg/L N-NH3') | contains(table_NH3_experiment_txt{i,5},'50 mg N-NH3/L')
                NH3_levels(i) = 3;
            else
                NH3_levels(i) = 1;
            end
        end
    end
    
%     Figure S8-1: Effect of NH3 salt addition on E. coli removal
%     performance at different pH normalized for 20°C
% THe data for this figure is first normalized, error bars are calculated
% based on the standard deviation for the available data set.
    decay_NH3_T_normalized = table_NH3_experiment(:,12)./(teta.^(table_NH3_experiment(:,9)-20));

    list_no_NH3 = NH3_levels == 1;
    list_low_NH3 = NH3_levels == 2;
    list_high_NH3 = NH3_levels == 3;

    decay_norm_no_NH3_pH7 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 7 & NH3_levels == 1);
    decay_norm_no_NH3_pH8 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 8 & NH3_levels == 1);
    decay_norm_no_NH3_pH9 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 9 & NH3_levels == 1);
    decay_norm_no_NH3_pH10 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 10 & NH3_levels == 1);

    decay_norm_low_NH3_pH7 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 7 & NH3_levels == 2);
    decay_norm_low_NH3_pH8 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 8 & NH3_levels == 2);
    decay_norm_low_NH3_pH9 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 9 & NH3_levels == 2);
    decay_norm_low_NH3_pH10 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 10 & NH3_levels == 2);

    decay_norm_hi_NH3_pH7 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 7 & NH3_levels == 3);
    decay_norm_hi_NH3_pH8 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 8 & NH3_levels == 3);
    decay_norm_hi_NH3_pH9 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 9 & NH3_levels == 3);
    decay_norm_hi_NH3_pH10 = decay_NH3_T_normalized(table_NH3_experiment(:,5) == 10 & NH3_levels == 3);

    bar_matrix = [mean(decay_norm_no_NH3_pH7) , mean(decay_norm_low_NH3_pH7) , mean(decay_norm_hi_NH3_pH7) ;...
                    mean(decay_norm_no_NH3_pH8) , mean(decay_norm_low_NH3_pH8) , mean(decay_norm_hi_NH3_pH8) ;...
                    mean(decay_norm_no_NH3_pH9) , mean(decay_norm_low_NH3_pH9) , mean(decay_norm_hi_NH3_pH9) ;...
                    mean(decay_norm_no_NH3_pH10) , mean(decay_norm_low_NH3_pH10) , mean(decay_norm_hi_NH3_pH10)];
    error_bar_matrix = [std(decay_norm_no_NH3_pH7) , std(decay_norm_low_NH3_pH7) , std(decay_norm_hi_NH3_pH7) ;...
                    std(decay_norm_no_NH3_pH8) , std(decay_norm_low_NH3_pH8) , std(decay_norm_hi_NH3_pH8) ;...
                    std(decay_norm_no_NH3_pH9) , std(decay_norm_low_NH3_pH9) , std(decay_norm_hi_NH3_pH9) ;...
                    std(decay_norm_no_NH3_pH10) , std(decay_norm_low_NH3_pH10) , std(decay_norm_hi_NH3_pH10)];

    figure(1), clf, hold on
    a = bar(bar_matrix);
    set(gca,'xtick', [1,2,3,4],'xticklabel',{'7','8','9','10'},'FontSize',fs);
    a(1).FaceColor = [0 0 0];
    a(2).FaceColor = [0.5 0.5 0.5];
    a(3).FaceColor = [1 1 1];


    xlabel('pH','FontSize',fs + 2)
    ylabel('{\itE. coli} decay rate (d^{-1})','FontSize',fs + 2)


    for i = 1:4
        errorbar(i - 0.225,bar_matrix(i,1),error_bar_matrix(i,1),error_bar_matrix(i,1),'k');
        errorbar(i,bar_matrix(i,2),error_bar_matrix(i,2),error_bar_matrix(i,2),'k');
        errorbar(i + 0.225,bar_matrix(i,3),error_bar_matrix(i,3),error_bar_matrix(i,3),'k');
    end

    legend1 = legend('No NH_3 addition','NH_3 at HRAP level','NH_3 at wastewater level','Location','northwest','Fontsize',fs);


    % Figure S8-2: Effect of NH3 salt addition on E; coli decay at pH 10
    % for different temperatures
    % Gathering of the data needed and formatted into a matrix adapated to
    % Matlab bar function
    table_NH3_experiment = [table_NH3_experiment , NH3_levels];
    table_NH3_experiment_pH_10 = table_NH3_experiment(table_NH3_experiment(:,5) == 10,:);

    decay_norm_no_NH3_T10 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 10 & table_NH3_experiment_pH_10(:,end) == 1,12);
    decay_norm_no_NH3_T25 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 25 & table_NH3_experiment_pH_10(:,end) == 1,12);
    decay_norm_no_NH3_T30 = (table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 30 & table_NH3_experiment_pH_10(:,end) == 1,12));
    decay_norm_no_NH3_T35 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 35 & table_NH3_experiment_pH_10(:,end) == 1,12);

    decay_norm_low_NH3_T10 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 10 & table_NH3_experiment_pH_10(:,end) == 2,12);
    decay_norm_low_NH3_T25 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 25 & table_NH3_experiment_pH_10(:,end) == 2,12);
    decay_norm_low_NH3_T30 = (table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 30 & table_NH3_experiment_pH_10(:,end) == 2,12));
    decay_norm_low_NH3_T35 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 35 & table_NH3_experiment_pH_10(:,end) == 2,12);

    decay_norm_hi_NH3_T10 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 10 & table_NH3_experiment_pH_10(:,end) == 3,12);
    decay_norm_hi_NH3_T25 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 25 & table_NH3_experiment_pH_10(:,end) == 3,12);
    decay_norm_hi_NH3_T30 = (table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 30 & table_NH3_experiment_pH_10(:,end) == 3,12));
    decay_norm_hi_NH3_T35 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 35 & table_NH3_experiment_pH_10(:,end) == 3,12);

    bar_matrix = [decay_norm_no_NH3_T10 , decay_norm_low_NH3_T10 , decay_norm_hi_NH3_T10 ;...
                    decay_norm_no_NH3_T25 , decay_norm_low_NH3_T25 , decay_norm_hi_NH3_T25 ;...
                    decay_norm_no_NH3_T30 , decay_norm_low_NH3_T30 , decay_norm_hi_NH3_T30 ;...
                    decay_norm_no_NH3_T35 , NaN , decay_norm_hi_NH3_T35];

    % Obtaining error bars from measurment uncertainty:
    err_decay_norm_no_NH3_T10 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 10 & table_NH3_experiment_pH_10(:,end) == 1,15);
    err_decay_norm_no_NH3_T25 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 25 & table_NH3_experiment_pH_10(:,end) == 1,15);
    err_decay_norm_no_NH3_T30 = (table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 30 & table_NH3_experiment_pH_10(:,end) == 1,15));
    err_decay_norm_no_NH3_T35 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 35 & table_NH3_experiment_pH_10(:,end) == 1,15);

    err_decay_norm_low_NH3_T10 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 10 & table_NH3_experiment_pH_10(:,end) == 2,15);
    err_decay_norm_low_NH3_T25 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 25 & table_NH3_experiment_pH_10(:,end) == 2,15);
    err_decay_norm_low_NH3_T30 = (table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 30 & table_NH3_experiment_pH_10(:,end) == 2,15));
    err_decay_norm_low_NH3_T35 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 35 & table_NH3_experiment_pH_10(:,end) == 2,15);

    err_decay_norm_hi_NH3_T10 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 10 & table_NH3_experiment_pH_10(:,end) == 3,15);
    err_decay_norm_hi_NH3_T25 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 25 & table_NH3_experiment_pH_10(:,end) == 3,15);
    err_decay_norm_hi_NH3_T30 = (table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 30 & table_NH3_experiment_pH_10(:,end) == 3,15));
    err_decay_norm_hi_NH3_T35 = table_NH3_experiment_pH_10(table_NH3_experiment_pH_10(:,9) == 35 & table_NH3_experiment_pH_10(:,end) == 3,15);

    error_bar_matrix = [err_decay_norm_no_NH3_T10 , err_decay_norm_low_NH3_T10 , err_decay_norm_hi_NH3_T10 ;...
                    err_decay_norm_no_NH3_T25 , err_decay_norm_low_NH3_T25 , err_decay_norm_hi_NH3_T25 ;...
                    err_decay_norm_no_NH3_T30 , err_decay_norm_low_NH3_T30 , err_decay_norm_hi_NH3_T30 ;...
                    err_decay_norm_no_NH3_T35 , NaN , err_decay_norm_hi_NH3_T35];

    figure(2), clf, hold on
    a = bar(bar_matrix);
    set(gca,'xtick', [1,2,3,4],'xticklabel',{'10','25','30','35'},'FontSize',fs);
    a(1).FaceColor = [0 0 0];
    a(2).FaceColor = [0.5 0.5 0.5];
    a(3).FaceColor = [1 1 1];


    xlabel('Temperature (°C)','FontSize',fs + 2)
    ylabel('{\itE. coli} decay rate (d^{-1})','FontSize',fs + 2)
    % Manually inserting error bars
    for i = 1:4
        errorbar(i - 0.225,bar_matrix(i,1),error_bar_matrix(i,1),error_bar_matrix(i,1),'k');
        errorbar(i,bar_matrix(i,2),error_bar_matrix(i,2),error_bar_matrix(i,2),'k');
        errorbar(i + 0.225,bar_matrix(i,3),error_bar_matrix(i,3),error_bar_matrix(i,3),'k');
    end

    legend1 = legend('No NH_3 addition','NH_3 at HRAP level','NH_3 at wastewater level','Location','northwest','Fontsize',fs);
end

%% Sunlight direct decay

if option_sun_direct_decay == 1
    % Reduction to data of interest:    
    list_outdoor = strcmp(table_results_txt{:,2} , 'Outdoor');
    list_matrix = strcmp(table_results_txt{:,4} , 'RO water');
    list_addition = strcmp(table_results_txt{:,5} , '') | strcmp(table_results_txt{:,5} , 'pH buffer');
    list_pH = not(table_results(:,7) > 8);
    list_sun = table_results(:,1) > 0;
    list_sun_direct_decay = list_outdoor & list_matrix & list_addition & list_pH & list_sun;


    light_irradiance = table_results(list_sun_direct_decay,1); % MJ/m2
    duration = table_results(list_sun_direct_decay,11); % d-1;
    light_intensity = (table_results(list_sun_direct_decay,1)*10^6)./(table_results(list_sun_direct_decay,11)*3600*24); % W/m2
    decay_rate = table_results(list_sun_direct_decay,12); % d-1
    log_decay = decay_rate.*duration/log(10);
    error_bars = table_results(list_sun_direct_decay,15);
    error_bars_LRV = error_bars.*table_results(list_sun_direct_decay,11);
    
    % Outlier removal: 
    lin_reg = fitlm(light_irradiance,log_decay,'Intercept',false);
    % fitlm(light_intensity,decay_rate,'Intercept',false)
    residuals = log_decay - lin_reg.Coefficients{1,1}*light_irradiance;
    
    
    [~,outliers_index] = rmoutliers(residuals,'grubbs');
    while max(outliers_index) == 1

        light_irradiance = light_irradiance(not(outliers_index));
        log_decay = log_decay(not(outliers_index));
        light_intensity = light_intensity(not(outliers_index));
        decay_rate = decay_rate(not(outliers_index));
        error_bars = error_bars(not(outliers_index));
        error_bars_LRV = error_bars_LRV(not(outliers_index));

        lin_reg = fitlm(light_irradiance,log_decay,'Intercept',false);
        % fitlm(light_intensity,decay_rate,'Intercept',false)
        residuals = log_decay - lin_reg.Coefficients{1,1}*light_irradiance;
        [~,outliers_index] = rmoutliers(residuals,'grubbs');
    end

    % Figure 2: Effect of sunlight dose on E. coli LRV at neutral pH during
    % laboratory assays
    figure(1), clf, hold on
    plot(light_irradiance,log_decay,'ok','MarkerSize',8,'LineWidth',2)
    x = linspace(0,8);
    y = lin_reg.Coefficients{1,1}*x;
    plot(x,y,'-k','LineWidth',2)
    % xlim([-10 250])
    % ylim([-10 250])
    errorbar(light_irradiance,log_decay,error_bars_LRV,'LineStyle','None','Color','k','Linewidth',1.2)
    axis = gca;
    axis.FontSize = fs - 3;
    xlabel('Sunlight radiation (MJ{\cdot}m^{-2})','FontSize',fs,'FontName','Arial','FontWeight','bold')
    ylabel('{\itE. coli} LRV','FontSize',fs,'FontName','Arial','FontWeight','bold')

    R_square = 1 - (sum((log_decay-light_irradiance*lin_reg.Coefficients{1,1}).^2)/sum(log_decay.^2))
    
    lin_reg_rate = fitlm(light_intensity,decay_rate,'Intercept',false);
    R_square_rate = 1 - (sum((decay_rate-light_intensity*lin_reg_rate.Coefficients{1,1}).^2)/sum(decay_rate.^2))
    
end

%% Sun pH combined effect
if option_sun_pH == 1

    % Reduction to data of interest:
    date_1 = '24/02/2017';
    date_2 = '28/02/2017';

    list_outdoor = strcmp(table_results_txt{:,2} , 'Outdoor');
    list_matrix = strcmp(table_results_txt{:,4} , 'RO water');
    list_addition = strcmp(table_results_txt{:,5} , '') | strcmp(table_results_txt{:,5} , 'pH buffer');
    list_sun = table_results(:,1) > 0;

    list_exp = strcmp(table_results_txt{:,1} ,date_1) | strcmp(table_results_txt{:,1} ,date_2);

    list_sun_pH_decay = list_outdoor & list_matrix & list_addition & list_sun;

    sunlight_spec_LRV = table_results(list_sun_pH_decay,12)/log(10).*table_results(list_sun_pH_decay,11)./table_results(list_sun_pH_decay,1);
    
    error_bars_LRV = table_results(list_sun_pH_decay,15).*table_results(list_sun_pH_decay,11);
    
    
    % Figure S11-1: Influence of pH on E. coli LRV normalized for sunlight
    % exposure during laboratory assays
    l = length(sunlight_spec_LRV);
    % A jitter is added to visually separate data measured at equal pH.
    position_jitter = randn(l);
    position_jitter = 2*position_jitter/max(position_jitter);
    
    figure(1), clf, hold on
    plot(table_results(list_sun_pH_decay,5) + position_jitter,sunlight_spec_LRV,'dk','MarkerSize',8,'MarkerFaceColor','k')
    errorbar(table_results(list_sun_pH_decay,5) + position_jitter,sunlight_spec_LRV,error_bars_LRV,'LineStyle','None','Color','k','Linewidth',1.2)
    xlim([6 11])
    ylim([-1 3])
    yticks(-1:1:3)
    axis = gca;
    axis.FontSize = fs - 3;
    xlabel('pH','FontSize',fs,'FontName','Arial','FontWeight','bold')
    ylabel({'{\itE. coli} sunlight specific LRV' ; '(LRV{\cdot}m^{2}{\cdot}MJ^{-1})'},'FontSize',fs,'FontName','Arial','FontWeight','bold')

    
    % Numerical investigation of potential synergetic effect of pH with
    % sunlight
    list_sun_pH_decay_spec = list_sun_pH_decay & list_exp;
    table_sun_pH_decay = table_results(list_sun_pH_decay_spec,:);
    % In this table, the first and second line added correspond to the
    % total of the 4th line (second data missed hence the decay was considered
    % on a two-sampling period. For the sake of analysis, both lines are put
    % together in a single experiment.

    table_sun_pH_decay = [table_sun_pH_decay(1,1) + table_sun_pH_decay(2,1) , table_sun_pH_decay(1,2:10) , ...
                            table_sun_pH_decay(1,11) + table_sun_pH_decay(2,11) ,...
                            (table_sun_pH_decay(1,12)*table_sun_pH_decay(1,11) + table_sun_pH_decay(2,12)*table_sun_pH_decay(2,11))/( table_sun_pH_decay(1,11) + table_sun_pH_decay(2,11)),...
                            table_sun_pH_decay(1,13:end);...
                            table_sun_pH_decay(3:end,:)];

    table_sun_pH_decay = sortrows(table_sun_pH_decay);

     % Ordering per comparable experiment: in order of data: sunlight
     % irradiance, high pH measured, temeprature high pH measured, decay neutral pH, decay
     % high pH, decay theoretical at high pH (equation 2)
     table_synergy = NaN(5,6);
     aux = table_sun_pH_decay(1,1); 
     k = 1;
     i = 1;
    while i <= size(table_sun_pH_decay,1)
        aux =  table_sun_pH_decay(i,1);
        if table_sun_pH_decay(i,5) == 7
            table_synergy(k,4) = table_sun_pH_decay(i,12);
        end
        if table_sun_pH_decay(i,5) == 10
            table_synergy(k,1:3) = [table_sun_pH_decay(i,1),table_sun_pH_decay(i,6),table_sun_pH_decay(i,10)];
            table_synergy(k,5) = table_sun_pH_decay(i,12);
            table_synergy(k,6) = 149000*1.14^(table_synergy(k,3)-20)*10^(table_synergy(k,2) - 14);
        end
        i = i+1;
        if isnan(table_synergy(k,:)) == [0 0 0 0 0 0]
            k = k + 1;
        end
    end

    synergy_effect = (table_synergy(:,5) - table_synergy(:,6) - table_synergy(:,4))./table_synergy(:,4);

    
    
    % Analysis in presence of photosensitizers
    % Only experiment from 28/02 yielded exploitable resutls so data is
    % reduced to this experiment
    list_exp = strcmp(table_results_txt{:,1} ,date_2);

    list_matrix = strcmp(table_results_txt{:,4} , 'Filtrated wastewater') | ...
                    strcmp(table_results_txt{:,4} , 'Filtrated HRAP broth');
    list_high_pH = table_results(:,5) == 10;


    list_sun_high_pH_PS_decay = list_matrix & list_outdoor & list_high_pH & list_addition & list_sun & list_exp ;
    list_sun_high_pH_noPS_decay = not(list_matrix) & list_outdoor & list_high_pH & list_addition & list_sun & list_exp ;

    table_sun_high_pH_PS = table_results(list_sun_high_pH_PS_decay,:);
    table_sun_high_pH_noPS = table_results(list_sun_high_pH_noPS_decay,:);
end

%% Exogenous photo-oxidation

if option_exogenousPO == 1
    % Reduction to data of interest:
    list_outdoor = strcmp(table_results_txt{:,2} , 'Outdoor');
    list_addition = strcmp(table_results_txt{:,5} , '') | strcmp(table_results_txt{:,5} , 'pH buffer');
    list_sun = table_results(:,1) > 0;
    list_pH = table_results(:,5) <= 8 | isnan(table_results(:,5)); 

    list_sun_PS = list_outdoor & list_addition & list_sun & list_pH;

    table_PS = table_results(list_sun_PS,:);
    table_PS_txt = table_results_txt(list_sun_PS,:);

    % A table is created in which each line gives experimental metadata and
    % the decay rate in RO water and a filtrate for simultaneous measurements

        % Col 1: type of filtrate: 1 for HRAP, 2 for WW
        % Col 2: 2nd type of filtrate if applicable
        % Col 3: sunlight irradiance
        % Col 4: duration
        % Col 5: decay rate RO water
        % Col 6: 2nd decay rate RO water if applicable
        % Col 7: decay rate filtrate
        % Col 8: 2nd decay rate filtrate if appicable
        % Col 9: error for first decay rate RO water if applicable
        % Col 10: error for second decay rate RO water if applicable
        % Col 11: error for first decay rate filtrate if applicable
        % Col 12: error for second decay rate filtrate if applicable
        

    irradiance_list = unique(table_PS(:,1)); 
    histogram_PS = NaN(length(irradiance_list),12);

    for i = 1:size(irradiance_list,1)
        list = table_PS(:,1) == irradiance_list(i);
        table = table_PS(list,:);
        table_txt = table_PS_txt(list,:);
        if size(table,1) >= 2
            histogram_PS(i,3:4) = [table(1,1) , table(1,11)];
            for k = 1 : size(table,1)
                % Sun intensity, duration
                % Matrix and decay rate
                if strcmp(table_txt{k,4},'Filtrated wastewater') == 1
                    if isnan(histogram_PS(i,1)) == 1
                        histogram_PS(i,1) = 2;
                        histogram_PS(i,7) = table(k,12);
                        histogram_PS(i,11) = table(k,15);
                    else
                        histogram_PS(i,2) = 2;
                        histogram_PS(i,8) = table(k,12);
                        histogram_PS(i,12) = table(k,15);
                    end
                elseif strcmp(table_txt{k,4},'Filtrated HRAP broth') == 1
                    if isnan(histogram_PS(i,1)) == 1
                        histogram_PS(i,1) = 1;
                        histogram_PS(i,7) = table(k,12);
                        histogram_PS(i,11) = table(k,15);
                    else
                        histogram_PS(i,2) = 1;
                        histogram_PS(i,8) = table(k,12);
                        histogram_PS(i,12) = table(k,15);
                    end
                else
                    if strcmp(table_txt{k,4},'RO water') == 1
                        if isnan(histogram_PS(i,5)) == 1
                            histogram_PS(i,5) = table(k,12);
                            histogram_PS(i,9) = table(k,15);
                        else
                            histogram_PS(i,6) = table(k,12);
                            histogram_PS(i,10) = table(k,15);
                        end
                    end
                end
            end     
        end
    end
        
     %  All data with no PS reactor is removed from the table

    list_keep = not(isnan(histogram_PS(:,1)) & isnan(histogram_PS(:,2))); 
    histogram_PS = histogram_PS(list_keep,:);
    [h_PS,p_PS,ci_PS,stats_PS] = ttest(nanmean(histogram_PS(:,5:6),2),nanmean(histogram_PS(:,7:8),2),'tail','left')% tail left so that the alternative hypothesis is that decay is lower without PS


    light_intensity = histogram_PS(:,3)*10^6./(histogram_PS(:,4)*24*3600);
    [light_intensity,index] = sortrows(light_intensity);
    histogram_PS = histogram_PS(index,:);

    labels = strings(1,length(light_intensity));
    for i = 1:length(light_intensity)
        labels(i) = strcat(num2str(light_intensity(i),3),' W{\cdot}m^{-2}');
    end
    % Figure S12-2: E. coli decay rate under natural sunlight in probable
    % presence of photosensitizers or absence at neutral pH.
    figure(1), clf, hold on
    b = bar(histogram_PS(:,5:8))
    b(1).FaceColor = [1 1 1];
    b(2).FaceColor = [1 1 1];
    b(3).FaceColor = [0.6 0.6 0.6];
    b(4).FaceColor = [0.6 0.6 0.6];

    a = gca;
    set(a,'xtick',1:9,'xticklabel',labels,'FontSize',fs,'XTickLabelRotation',60,'FontName','Arial','FontWeight','bold');
    ylabel({'{\itE. coli} decay rate' ; '(d^-^1)'},'FontSize',fs,'FontName','Arial','FontWeight','bold')

    % Addition of * on wastewater filtrates
    eps = 0.09;
    X = [1 + eps ,  1 + 3*eps ,2 + eps ,  2 + 3*eps , 6 + eps ,  6 + 3*eps ];
    Y  = [histogram_PS(1,7) , histogram_PS(1,8) , histogram_PS(2,7) ,  histogram_PS(2,8) , histogram_PS(6,7) ,  histogram_PS(6,8)];
    plot(X,Y,'*k','MarkerSize',8,'LineWidth',1)
    % Addition of error bars
    errorbar(b(1).XData - 0.27,b(1).YData , histogram_PS(:,9),'LineStyle','None','Color','k','Linewidth',1.2)
    errorbar(b(2).XData - 0.08,b(2).YData , histogram_PS(:,10),'LineStyle','None','Color','k','Linewidth',1.2)
    errorbar(b(3).XData + 0.08,b(3).YData , histogram_PS(:,11),'LineStyle','None','Color','k','Linewidth',1.2)
    errorbar(b(4).XData + 0.27,b(4).YData , histogram_PS(:,12),'LineStyle','None','Color','k','Linewidth',1.2)    
    
    legend([b(1) b(3)],{'Absence of photosensitizers','Probable presence of photosensitizers'},'FontSize',fs-2,'FontName','Arial','Location','northwest','FontWeight','normal')

    %% PS pH interaction
    
    if option_exogenousPO_pH10 == 1
        % Reduction to data of interest:
        list_outdoor = strcmp(table_results_txt{:,2} , 'Outdoor');
        list_matrix_RO = strcmp(table_results_txt{:,4} , 'RO water');
        list_matrix_PS = strcmp(table_results_txt{:,4} , 'Filtrated HRAP broth');
        list_addition = strcmp(table_results_txt{:,5} , 'pH buffer');
        list_sun = table_results(:,1) > 0;
        list_pH = table_results(:,5) >= 8; 

        list_sun_PS_pH10_RO = list_outdoor & list_addition & list_sun & list_pH & list_matrix_RO;
        list_sun_PS_pH10_PS = list_outdoor & list_addition & list_sun & list_pH & list_matrix_PS;
    
    
        table_PS_pH10_RO = table_results(list_sun_PS_pH10_RO,:);
        table_PS_pH10_RO = [table_PS_pH10_RO(3:end,:),table_PS_pH10_RO(3:end,1)*10^6./table_PS_pH10_RO(3:end,11)/24/3600];
        table_PS_pH10_PS = [table_results(list_sun_PS_pH10_PS,:),table_results(list_sun_PS_pH10_PS,1)*10^6./table_results(list_sun_PS_pH10_PS,11)/24/3600];

        table_PS_pH10_RO = sortrows(table_PS_pH10_RO,size(table_PS_pH10_RO,2));   
        table_PS_pH10_PS = sortrows(table_PS_pH10_PS,size(table_PS_pH10_PS,2));

        % Figure S11-2: E. colid decay rate in RO water or filtrates from
        % HRAP both under natural sunlight at pH 10       
        figure(1), clf, hold on
        b = bar([table_PS_pH10_RO(:,12),table_PS_pH10_PS(:,12)])
        b(1).FaceColor = [1 1 1];
        b(2).FaceColor = [0.6 0.6 0.6];

        a = gca;
        labels = {'808 W{\cdot}m^{-2}','872 W{\cdot}m^{-2}','952 W{\cdot}m^{-2}'};
        set(a,'xtick',1:9,'xticklabel',labels,'FontSize',fs,'XTickLabelRotation',60,'FontName','Arial','FontWeight','bold');
        ylabel({'{\itE. coli} decay rate' ; '(d^-^1)'},'FontSize',fs,'FontName','Arial','FontWeight','bold')
        errorbar(b(1).XData - 0.15,b(1).YData , table_PS_pH10_RO(:,15),'LineStyle','None','Color','k','Linewidth',1.2)
        errorbar(b(2).XData + 0.15,b(2).YData , table_PS_pH10_RO(:,15),'LineStyle','None','Color','k','Linewidth',1.2)
    
        legend([b(1) b(2)],{'Absence of photosensitizers','Probable presence of photosensitizers'},'FontSize',fs-2,'FontName','Arial','Location','northeast','FontWeight','normal')
        
    end
    
end