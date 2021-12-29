%% Main_Figure3

% This script aims at generating Figure 3 of the manuscript "Study from microcosms and mesocosms reveals
% Escherichia coli removal in High Rate Alage Ponds is primarily caused by
% dark decay", titled "Measured versus computed E. coli cell couns during
% bench assays using best fit model parameters".

fs = 18;

load('../Datasets/bench_data_fitting_best_fit.mat');
display(R2_MC)
display(alpha_sun_MC)
display(k_nat_MC)
display(teta_nat_MC)
display(k_pH_MC)
display(teta_pH_MC)
display(mean(abs(coli_MC_logcounts - coli_MC_logmodel)./coli_MC_logcounts))

figure(1), clf, hold on
plot(coli_MC_logmodel,coli_MC_logcounts,'dk','MarkerSize',8,'LineWidth',2)
xlim([0 8])
ylim([0 8])
x = linspace(0,8);
y = linspace(0,8);
plot(x,y,'-k','LineWidth',3)
errorbar(coli_MC_logmodel,coli_MC_logcounts,0.2*coli_MC_logcounts,0.2*coli_MC_logcounts,'Color','k','LineStyle','None','LineWidth',1)
ax = gca;
ax.FontSize = fs - 2
xlabel({'Predicted {\itE. coli} cell count' ; '(log MPN{\cdot}100 mL^{-1})'},...
    'FontSize',fs,'FontName','Arial','FontWeight','bold')
ylabel({'Measured {\itE. coli} cell count' ;  '(log MPN{\cdot}100 mL^{-1})'},...
    'FontSize',fs,'FontName','Arial','FontWeight','bold')