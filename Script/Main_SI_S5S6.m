%% Main_SI_S5S6

% This script aims at generating the matlab figures for SI-5 and SI 6 as
% well as providing the analysis performed to determine approximate
% distributio of the Monte Carlo generated data sets for fitted parameters.


fs = 18;
option_MCFit = 1;
option_find_dsitribution = 1;
option_figure_S6_2 = 1;

load('../Datasets/bench_data_fitting_MC_analysis.mat')

if option_MCFit == 1
        prctile(k_nat_MC,50), mean(k_nat_MC) , prctile(k_nat_MC,5) , prctile(k_nat_MC,95),
        prctile(teta_nat_MC,50) , mean(teta_nat_MC), prctile(teta_nat_MC,5), prctile(teta_nat_MC,95),
        prctile(k_pH_MC,50) , mean(k_pH_MC) , prctile(k_pH_MC,5) , prctile(k_pH_MC,95),
        prctile(teta_pH_MC,50) , mean(teta_pH_MC) , prctile(teta_pH_MC,5) , prctile(teta_pH_MC,95),
        prctile(alpha_sun_MC,50) , mean(alpha_sun_MC) , prctile(alpha_sun_MC,5), prctile(alpha_sun_MC,95)
                                  
    
        figure(100), clf, hold on
            subplot(2,3,1)
            h1 = boxplot(k_nat_MC','colors','k','width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3);
            set(gca,'xticklabel',{[]},'FontSize',fs-4)
            ylabel('k_{20}^{dark} (d^-^1)','FontSize',fs)
            set(gca,'LineWidth',2)
            set(h1,'LineWidth',1)

            subplot(2,3,4)
            h2 = boxplot(teta_nat_MC','colors','k','width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3)
            set(gca,'xticklabel',{[]},'FontSize',fs-4)
            ylabel('\theta^{dark}','FontSize',fs)
            set(gca,'LineWidth',2)
            set(h2,'LineWidth',1)

            subplot(2,3,2)
            h3 = boxplot(k_pH_MC','colors','k','width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3)
            set(gca,'xticklabel',{[]},'FontSize',fs-4)
            ylabel('k_{20}^{pH} (d^-^1)','FontSize',fs)
            set(gca,'LineWidth',2)
            set(h3,'LineWidth',1)

            subplot(2,3,5)
            h4 = boxplot(teta_pH_MC,'colors','k','width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3)
            set(gca,'xticklabel',{[]},'FontSize',fs-4)
            ylabel('\theta^{pH}','FontSize',fs)
            set(gca,'LineWidth',2)
            set(h4,'LineWidth',1)

            subplot(2,3,3)
            h5 = boxplot(alpha_sun_MC,'colors','k','width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3)
            set(gca,'xticklabel',{[]},'FontSize',fs-4)
            ylabel('\alpha (m^2.W^-^1.d^-^1)','FontSize',fs)
            set(gca,'LineWidth',2)
            set(h4,'LineWidth',1)
end

if option_find_dsitribution == 1

    % Dark decay

    nbins = 100;

    h = hist(k_nat_MC,nbins);
    x0 = min(k_nat_MC); 
    x1 = max(k_nat_MC);
    x = x0:(x1-x0)/(nbins-1):x1;

    coeff_ev  = evfit(k_nat_MC);
    p = evpdf(x,coeff_ev(1),coeff_ev(2));

    figure(1), clf, hold on
    bar(x,h/length(k_nat_MC))
    line(x,p,'color','r')


    %% pH decay

    % An extreme outlier was found in the distribution of teta_pH which
    % compromised the possibility of building an extreme value distribution fit
    % satisfyingly matching the observed distribution. THis single extreme
    % value was disregarded in the building of a test distribution.

    teta_test_dist = teta_pH_MC;
    teta_test_dist(teta_test_dist == max(teta_test_dist)) = [];
    x0_fit = min(teta_test_dist);
    x1_fit = max(teta_test_dist);

    x_fit = x0_fit:(x1_fit-x0_fit)/(nbins-1):x1_fit;

    coeff_ev  = evfit(teta_test_dist);
    p = evpdf(x_fit,coeff_ev(1),coeff_ev(2));

    x0 = min(teta_pH_MC);
    x1 = max(teta_pH_MC);
    x = x0:(x1-x0)/(nbins-1):x1;
    h = hist(teta_pH_MC,nbins);

    figure(2), clf, hold on
    bar(x,h/100)
    line(x_fit,p,'color','r')

    %% Sun decay

    % another extreme value fit on alpha sun MC when disregarding the 0

    nbins = 100;
    alpha_fit = alpha_sun_MC(alpha_sun_MC > 0);
    x0 = min(alpha_fit);
    x1 = max(alpha_fit);
    x = x0:(x1-x0)/(nbins-1):x1;

    h = hist(alpha_fit,nbins);

    coeff_ev  = fitdist(alpha_fit','exponential');
    p = exppdf(x,coeff_ev.mu);

    figure(3), clf, hold on
    bar(x,h/length(h))
    line(x,p/10,'color','r')
end

if option_figure_S6_2 == 1
    clear all
    load('../Datasets/MC_analysis_contribution_uncertainty.mat')

    k_pH_MC_out = k_pH_MC;
    teta_pH_MC_out = teta_pH_MC;
    k_dark_MC_out = k_dark_MC;
    teta_dark_MC_out = teta_dark_MC;
    alpha_sun_MC_out = alpha_sun_MC;
    
    
    load('../Datasets/bench_data_fitting_MC_analysis.mat')
    
    index_1 = find(k_nat_MC > 9.98 & k_nat_MC < 58.4);
    index_2 = find(teta_nat_MC >= 1 & teta_nat_MC < 1.10);
    index_3 = find(k_pH_MC > 2220 & k_pH_MC < 25800);
    index_4 = find(teta_pH_MC > 1.19 & teta_pH_MC < 1.49);
    index_5 = find(alpha_sun_MC >= 0 & alpha_sun_MC < 0.491);

    k_nat_MC = k_nat_MC(index_1);
    teta_nat_MC = teta_nat_MC(index_2);
    k_pH_MC = k_pH_MC(index_3);
    teta_pH_MC = teta_pH_MC(index_4);
    alpha_sun_MC = alpha_sun_MC(index_5);

    k_pH_MC_in = k_pH_MC;
    teta_pH_MC_in = teta_pH_MC;
    k_dark_MC_in = k_nat_MC;
    teta_dark_MC_in = teta_nat_MC;
    alpha_sun_MC_in = alpha_sun_MC;
    
    %% Plot
    fs = 14;
    figure(1), clf, hold on
    %%
    subplot(2,5,1)
    histogram(k_dark_MC_in,100,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
    title({'k_{20}^{dark}'; '(d^{-1})'},'FontSize',fs)
    ylabel({'Monte Carlo outputs';'from uncertainty analysis';'Frequency'},'FontWeight','bold','FontSize',fs)
    ylim([0 0.025])
    
    subplot(2,5,6)
    histogram(k_dark_MC_out,100,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
    ylabel({'Monte Carlo inputs for' ; 'mechanism significance analysis';'Frequency'},'FontWeight','bold','FontSize',fs-2)
    ylim([0 0.025])
    %%
    subplot(2,5,2)
    histogram(teta_dark_MC_in,100,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
    title({'\theta^{dark}';''},'FontSize',fs)
	ylim([0 0.6])

    subplot(2,5,7)
    histogram(teta_dark_MC_out,100,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
	ylim([0 0.6])

    %%
    subplot(2,5,3)
    histogram(k_pH_MC_in,100,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
    title({'k_{20}^{pH}';'(d^{-1})'},'FontSize',fs)
	ylim([0 0.25])

    subplot(2,5,8)
    histogram(k_pH_MC_out,50,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
    ylim([0 0.25])

    %%
    subplot(2,5,4)
    histogram(teta_pH_MC_in,100,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
    title({'\theta^{pH}';''},'FontSize',fs)
	ylim([0 0.06])

    subplot(2,5,9)
    histogram(teta_pH_MC_out,50,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
	ylim([0 0.06])

    %%
    subplot(2,5,5)
    histogram(alpha_sun_MC_in,100,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
    title({'\alpha' ; '(m^2.W^{-1}.d^{-1})'},'FontSize',fs)
	ylim([0 0.6])

    subplot(2,5,10)
    histogram(alpha_sun_MC_out,100,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
	ylim([0 0.6])

end




