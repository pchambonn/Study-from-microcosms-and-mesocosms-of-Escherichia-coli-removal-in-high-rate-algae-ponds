%% Main_monteCarlo_uncertaintyFittedParameters

% This script aims at performing the uncertainty analysis of the fitted
% parameters calculated for the model developed in the publication "Study
% from microcosms and mesocosms reveals Escherichia coli removal in High
% Rate Alage Ponds is primarily caused by dark decay".

% The full uncertainty analysis is based on a Monte Carlo method by
% generating random data for the uncertain parameter or measurement used
% for the model development and recompute fitting parameters following the
% base case scenatio routine.

% The uncertainty computed takes into account the uncertainty on data
% obtained for sunlight radiation, pH, temperature, E. coli cell counts,
% light extinction coefficient, reactor depth, and fitted parameters
% uncertainty.

%% Monte Carlo Parameters
N_MC = 2000; 

err_pH = 0.01; % pH unit
err_temp = 0.1; % °C
err_sun = 0.1; % relative, *100%
err_coli = 0.2; % relative, *100%
err_sigma = 10; % m-1
err_depth = 0.03; % m


% Distribution parameters for uncertainty of fitting parameters: see
% supplementary information S6 for details
err_k_pH = 0.109;
mu_teta_pH = 1.43;
sigma_teta_pH = 6.56*10^(-2);

mu_k_dark = 44.5;
sigma_k_dark = 12.4;
err_teta_dark = 0.0202;

mu_alpha_sun = 0.2388;

%% Initialization

load('../Datasets/bench_data_import_full_workspace.mat');

% Creation of final tables used for the storage of data generated for Monte
% Carlo Analysis

Rel_contr_sun = NaN(1,N_MC); % Relative contribution of sunlight direct decay to overall E. coli decay (%)
Rel_contr_pH = NaN(1,N_MC); % Relative contribution of alkaline pH induced decay to overall E. coli decay (%)
Rel_contr_dark = NaN(1,N_MC); % Relative contribution of dark decay to overall E. coli decay (%)

N_decay_dark_MC = NaN(1,N_MC); % Total number of E. coli cell decayed through dark decay (MPN)
N_decay_pH_MC = NaN(1,N_MC); % Total number of E. coli cell decayed through alkaline pH induced decay (MPN)
N_decay_sun_MC = NaN(1,N_MC); % Total number of E. coli cell decayed through sunlight direct decay (MPN)

% Creation of tables used for storage of the inputs for Monte Carlo
% simulations

coli_MC_logcounts = zeros(nData,N_MC); % Table of the E. coli cell count log transformed measured (log MPN/100mL) 
time_data_MC = cell(nExp,N_MC);  % Table of time data (d)
sun_data_MC = cell(nExp,N_MC); % Table of sunlight direct intensity (W/m2)
DO_data_MC = cell(nExp,N_MC); % Table of dissolved oxygen (mg/L)
pH_data_MC = cell(nExp,N_MC); % Table of pH
temp_data_MC = cell(nExp,N_MC); % Table of temperature (°C)
coli_data_MC = cell(nExp,N_MC); % Table of E. coli cell count measured (MPN/100mL)
coli_model_MC = cell(nExp,N_MC); % Table of E. coli cell count modelled (MPN/100mL)
sigma_MC = NaN(N_MC,4); % Light attenuation coefficient (m-1)
d_MC = NaN(N_MC,1); % Reactor depth (m)

% Creation of tables used for storage of the fitted parameters necessary
% for Monte Carlo analysis

k_pH_MC = NaN(N_MC,1); % coefficient of alkaline pH induced decay at 20°C and pH 14 (d-1)
teta_pH_MC = NaN(N_MC,1); % Arrhenius coefficient of alkaline pH induced decay (-)
k_dark_MC = NaN(N_MC,1); % coefficient of dark decay at 20°C (d-1)
teta_dark_MC = NaN(N_MC,1); % Arrhenius coefficient of dark decay (-)
alpha_sun_MC = NaN(N_MC,1); % Coefficient of direct sunlight decay (m2/W/d)



    
%% Monte Carlo loop

for i_MC = 1:N_MC
    
    % Initialization The total number of E. coli cell decayed are set to 0
    N_decay_dark_MC(i_MC) = 0;
    N_decay_pH_MC(i_MC) = 0;
    N_decay_sun_MC(i_MC) = 0;
    
    % Input data are given values from experimental results and randomized
    % if needed
    for i = 1:nExp
        time_data_MC{i,i_MC} = timeDataFile{i};
        sun_data_MC{i,i_MC} = sunDataFile{i}.*(1+err_sun/2*randn(nDataFile{i},1));
        pH_data_MC{i,i_MC} = pHDataFile{i}+err_pH/2*randn(nDataFile{i},1);
        temp_data_MC{i,i_MC} = tempDataFile{i}+err_temp/2*randn(nDataFile{i},1);
        coli_data_MC{i,i_MC} = 10.^(log10(coliDataFile{i}).*(1+err_coli/2*randn(nDataFile{i},1)));
        coli_model_MC{i,i_MC} = NaN(nDataFile{i},1); % creation of table of data output
    end    
    
    sigma_MC(i_MC,:) = max(0.01,sigma0 + err_sigma/2*randn);
    d_MC(i_MC) = d0 + err_depth/2*randn;
    
    % Generation of randomized values for the fitted parameters:
    teta_pH_MC(i_MC) = random('Extreme value',mu_teta_pH,sigma_teta_pH);
    k_pH_MC(i_MC) = 10.^(-3.7302*teta_pH_MC(i_MC)+8.8475 + random('Normal',0,err_k_pH));
    
    k_dark_MC(i_MC) = random('Extreme value',mu_k_dark,sigma_k_dark);
    if k_dark_MC(i_MC) > 0
        bern = randi([0 1]);
        teta_dark_MC(i_MC) = max(1,bern + (1-bern)*(1.2356*(k_dark_MC(i_MC)^(-0.051)) + random('Normal',0,err_teta_dark)));
    else
        teta_dark_MC(i_MC) = 1.31 + random('Normal',0,err_teta_dark);
    end
    alpha_sun_MC(i_MC) = randi([0 1])*random('Exponential',mu_alpha_sun);
   
    k_pH = k_pH_MC(i_MC);
    teta_pH = teta_pH_MC(i_MC);
    alpha_sun = alpha_sun_MC(i_MC);
    k_nat = k_dark_MC(i_MC);
    teta_nat = teta_dark_MC(i_MC);

    % In the following, modelled E. coli cells based on the fitted
    % parameters and input parameters/variables is performed as already
    % described in the Monte Carlo analysis to determine uncertainty on the
    % fitted parameters.
    coli_fit = [];
    N_decay_dark = 0;
    N_decay_pH = 0;
    N_decay_sun = 0;
      
    for i_model = 1:nExp
        coli_model_MC{i_model,i_MC} = NaN(nDataFile{i_model},1);

        k = 1;
        while isnan(coli_data_MC{i_model,i_MC}(k)) && k < nDataFile{i_model}
            k = k + 1;
        end
        coli_model_MC{i_model,i_MC}(k) = coli_data_MC{i_model,i_MC}(k);
        

        for m = k + 1:nDataFile{i_model}
            k_nat_inst = k_nat*teta_nat^(temp_data_MC{i_model,i_MC}(m-1)-20);
            k_pH_inst = k_pH*teta_pH^(temp_data_MC{i_model,i_MC}(m-1)-20)*10^(pH_data_MC{i_model,i_MC}(m-1)-14);
            k_sun_inst = alpha_sun*sun_data_MC{i_model,i_MC}(m-1)/(sigma_MC(i_MC,sigmaType(i_model))*d_MC(i_MC))*(1-exp(-sigma_MC(i_MC,sigmaType(i_model))*d_MC(i_MC)));
            
            % The total amount of E. coli cells decayed through each
            % mechanism are computed below.
            N_decay_dark = N_decay_dark + max(1,coli_model_MC{i_model,i_MC}(m-1)*k_nat_inst*(time_data_MC{i_model,i_MC}(m)-time_data_MC{i_model,i_MC}(m-1)));
            N_decay_pH = N_decay_pH + max(1,coli_model_MC{i_model,i_MC}(m-1)*k_pH_inst*(time_data_MC{i_model,i_MC}(m)-time_data_MC{i_model,i_MC}(m-1)));
            N_decay_sun = N_decay_sun + max(1,coli_model_MC{i_model,i_MC}(m-1)*k_sun_inst*(time_data_MC{i_model,i_MC}(m)-time_data_MC{i_model,i_MC}(m-1)));
        
            coli_model_MC{i_model,i_MC}(m) = max(1,coli_model_MC{i_model,i_MC}(m-1)*(1 - (k_nat_inst + k_pH_inst + k_sun_inst)*(time_data_MC{i_model,i_MC}(m)-time_data_MC{i_model,i_MC}(m-1))));
        end

        index_coli = find(~isnan(coli_data_MC{i_model,i_MC}));
        coli_fit = [coli_fit ; coli_data_MC{i_model,i_MC}(index_coli), coli_model_MC{i_model,i_MC}(index_coli), i_model*ones(length(index_coli),1)];
    end
    
    N_decay_dark_MC(i_MC) = N_decay_dark;
    N_decay_pH_MC(i_MC) = N_decay_pH;
    N_decay_sun_MC(i_MC) = N_decay_sun;
    
    Rel_contr_sun(i_MC) = N_decay_sun/(N_decay_sun+N_decay_pH+N_decay_dark)*100;
    Rel_contr_pH(i_MC) = N_decay_pH/(N_decay_sun+N_decay_pH+N_decay_dark)*100;
    Rel_contr_dark(i_MC) = N_decay_dark/(N_decay_sun+N_decay_pH+N_decay_dark)*100;

end

%% Filtration of results by filtrating parameters only in the 95% uncertainty range
% As described in the main manuscript, only results emanating from within
% the 95% uncertainty range of the fitted parameters are kept for
% interpretation.
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
