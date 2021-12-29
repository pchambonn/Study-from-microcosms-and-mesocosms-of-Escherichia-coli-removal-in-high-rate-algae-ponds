%% Main_Figure4

% This script aims at generating Figure 4 of the manuscript "Study from
% microcosms and mesocosms reveals Escherichia coli removal in High Rate
% Alage Ponds is primarily caused by dark decay" titled "Relative
% contribution of removal mechanisms contributing to E. coli decay during
% bench assays".

load('../Datasets/MC_analysis_contribution_uncertainty.mat')

% The parameters outside of the 5-95 percentiles are filtered out (see
% Table 3 of main manuscript)
index_1 = find(k_dark_MC > 9.98 & k_dark_MC < 58.4);
index_2 = find(teta_dark_MC >= 1 & teta_dark_MC < 1.10);
index_3 = find(k_pH_MC > 2220 & k_pH_MC < 25800);
index_4 = find(teta_pH_MC > 1.19 & teta_pH_MC < 1.49);
index_5 = find(alpha_sun_MC >= 0 & alpha_sun_MC < 0.491);

index_f = intersect(index_1,intersect(index_2,intersect(index_3,intersect(index_4,index_5))));

k_pH_MC = k_pH_MC(index_f);
teta_pH_MC = teta_pH_MC(index_f);
k_dark_MC = k_dark_MC(index_f);
teta_dark_MC = teta_dark_MC(index_f);
alpha_sun_MC = alpha_sun_MC(index_f);

Rel_contr_sun = Rel_contr_sun(index_f);
Rel_contr_dark = Rel_contr_dark(index_f);
Rel_contr_pH = Rel_contr_pH(index_f);

X = [Rel_contr_dark;Rel_contr_pH;Rel_contr_sun];

fs = 20;

figure(1), clf, hold on
h1 = boxplot(X','colors','k','width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3,'whisker',0.7193);
set(gca,'xticklabel',[],'FontSize',fs-2)
ylim([0 100])
% ylabel({'Mechanism relative contribution to overall' ; '{\itE. coli} decay in HRAP microcosms (%)'},'FontSize',fs,'FontWeight','bold')
ylabel({'Mechanism relative contribution' ;  'to total {\itE. coli} decay (%)'},'FontSize',fs + 2,'FontWeight','bold', 'FontName', 'Arial')

xlabel({'Uncharacterized     pH-induced   Sunlight-mediated' ;...
        'dark decay             toxicity              damage   '}, 'FontName', 'Arial')

set(gca,'LineWidth',3)
set(h1,'LineWidth',2)
yticks(0:20:100)

