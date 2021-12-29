%% Main_modelFit_monteCarloUncertainty_sensitivityAnalysis
% This script aims at performing, as part of the publication "Study from
% microcosms and mesocosms reveals Escherichia coli removal in High Rate
% Alage Ponds is primarily caused by dark decay": 
    % 1) the fitting of E. coli decay as developed from lab assays to the
    % data from bench assays
    % 2) the uncertainty analysis of the fitting parameters
    % 3) the sensitivity analysis of the fitting parameters

% The fitting for base case scenario is performed using a Monte Carlo loop
% of only one loop, setting all errors to 0.

% The full uncertainty analysis is based on a Monte Carlo method by
% generating random data for the uncertain parameter or measurement used
% during model development and recompute fitting parameters following the
% base case scenatio routine (1).

% The sensitivity analsyis is performed using only one Monte Carlo loop
% setting all errors but one to zero, which is forced to the desired value,
% and retrieving the fitted values for each individual error tested.

% The uncertainty computed takes into account the uncertainty on data
% obtained for sunlight radiation, pH, temperature, E. coli cell counts,
% light extinction coefficient, and reactor depth. 

%% Monte Carlo Parameters
N_MC = 1000; % This number can be set to 1 if only the base case scenario is investigated
option_error = 1; % This number is set to 0 for base case scenario and 1 for Monte Carlo analysis

% Initialization parameters: based on lab experiment data
k_pH_i = 149000;
teta_pH_i = 1.14;
alpha_sun_i = 0.0626;

N_res_i = 50; % resolution used in the initial table of values
N_res_opt = 10; % Resolution used in the start of the second optimization loop

% Below are errors equaled to 0 to test for
if option_error == 0
    err_pH = 0;
    err_temp = 0;
    err_sun = 0;
    err_coli = 0;
    err_sigma = 0;
    err_depth = 0;
end

if option_error == 1
    err_pH = 0.01; % pH unit
    err_temp = 0.1; % °C
    err_sun = 0.1; % relative, *100%
    err_coli = 0.2; % relative, *100%
    err_sigma = 10; % m-1
    err_depth = 0.03; % m
end

%% Initialization
% Bench data results are loaded
load('../Datasets/bench_data_import_full_workspace.mat');

% Final tables of data generated for Monte Carlo Analysis
    % A list of tables for all outputs from Monte Carlo analysis is generated
    % below, in which is stored in each line the output calculated at each iteration

R2_MC = NaN(1,N_MC); % Coefficients of determination
SSR_MC_neut = NaN(1,N_MC); % Sum of Squared Residuals on the subset of neutral pH 
SSR_MC_alk = NaN(1,N_MC); % Sum of Squared Residuals on the subset of alkaline pH
SSR_MC = NaN(1,N_MC); % Sum of Squared Residuals
coli_MC_logmodel = NaN(nData,N_MC); % Table of the E. coli log transformed cell count modelled (i.e.) fitted at each iteration
k_nat_MC = NaN(1,N_MC); % Decay coefficient for dark decay at 20°C 
teta_nat_MC = NaN(1,N_MC); % Arrhenius coefficient for dark decay
k_pH_MC = NaN(1,N_MC); % Decay coefficient for alkaline pH induced decay at 20°C and pH 14
teta_pH_MC = NaN(1,N_MC); % Arrhenius coefficient for alkaline pH induced decay
alpha_sun_MC = NaN(1,N_MC); % Decay coefficient for direct sunlight decay
    
% Creation of the input tables which will be used for each iteration:
    % A cell array is created for each input necessary. The first index
    % corresponds to the experiment (i.e. individual reactor run) and the
    % seconf index correspond to the Monte Carlo simulation index.

coli_MC_logcounts = zeros(nData,N_MC); % E. coli cell count (log transformed MPN/100 mL) but randomized for each iteration.

time_data_MC = cell(nExp,N_MC); % time data (d)
sun_data_MC = cell(nExp,N_MC); % direct sunlight intensity data (W/m2) 
pH_data_MC = cell(nExp,N_MC); % pH data
temp_data_MC = cell(nExp,N_MC); % temperature data (°C)
coli_data_MC = cell(nExp,N_MC); % measured E. coli cell count (MPN/100 mL)
coli_model_MC = cell(nExp,N_MC); % prdicted E. coli cell count (MPN/100 mL)

% Parameters of Monte Carlo
sigma_MC = NaN(N_MC,4); % light extinction coefficient (m-1)
d_MC = NaN(N_MC,1); % reactor depth (m)

% In the following, all the data randomized for the purpose of Monte Carlo
% simulation is initialized according to the laws given in supplementary
% information S5. The cell array for modelled E. coli cell counts is
% initialized with NaN values.

for i_MC = 1:N_MC  

    for i = 1:nExp
        time_data_MC{i,i_MC} = timeDataFile{i}; % Time data is not randomized.
        sun_data_MC{i,i_MC} = sunDataFile{i}.*(1+err_sun/2*randn(nDataFile{i},1));
        pH_data_MC{i,i_MC} = pHDataFile{i}+err_pH/2*randn(nDataFile{i},1);
        temp_data_MC{i,i_MC} = tempDataFile{i}+err_temp/2*randn(nDataFile{i},1);
        coli_data_MC{i,i_MC} = 10.^(log10(coliDataFile{i}).*(1+err_coli/2*randn(nDataFile{i},1)));
        coli_model_MC{i,i_MC} = NaN(nDataFile{i},1); % creation of table of data output
    end    

    sigma_MC(i_MC,:) = max(0.01,sigma0 + err_sigma/2*randn);
    d_MC(i_MC) = d0 + err_depth/2*randn;
end

%% Monte Carlo Loop: fitting of the randomized data
for i_MC = 1:N_MC
% Values for k_pH and teta_pH are pre-assigned to generate a first set of
% modelled values based on lab experiments results

    k_pH = k_pH_i;
    teta_pH = teta_pH_i;
    alpha_sun = alpha_sun_i;

    %% Model fitting
    % The general fitting strategy is detailed in the manuscript.

    % The SSR is minimized by finding optimum values for k_pH, teta_pH,
    % k_nat, teta_nat, and alpha_sun. The optimum is considered reached
    % when SSR do not vary by more than 1%
    
    % To achieve so, a list of realistic values of the fitted parameters
    % are created with a resolution N_res_i. The set of values minimizing
    % the coefficient of determination are saved.
    
    % Initialization fo the resolution of the algorithm.
    N_res = N_res_i;


    %% Model pre-fitting
    
    % List of values for fitted parameters are created based on
    % realistically expected max, min, and the desited resolution.
    
    min_k_nat = 1;
    max_k_nat = 100;
    k_nat_table = min_k_nat:(max_k_nat-min_k_nat)/N_res:max_k_nat;


    min_teta_nat = 1;
    max_teta_nat = 2; %5
    teta_nat_table = min_teta_nat:(max_teta_nat-min_teta_nat)/N_res:max_teta_nat;


    min_k_pH = 1;
    max_k_pH = 250000;
    k_pH_table = min_k_pH:(max_k_pH-min_k_pH)/N_res:max_k_pH;


    min_teta_pH = 1;
    max_teta_pH = 2; %5
    teta_pH_table = min_teta_pH:(max_teta_pH-min_teta_pH)/N_res:max_teta_pH;

    min_alpha_sun = 0.00001;
    max_alpha_sun = 1;
    alpha_sun_table = [0,min_alpha_sun:(max_alpha_sun-min_alpha_sun)/N_res:max_alpha_sun];
    
    
    % A first screening to find minimum of SSR are performed
    % Impossible values are given to the SSR to make sure the while loop is
    % not ignored in the first time.
    SSR_tot_0 = -1; 
    SSR_tot = 0;

    while abs(SSR_tot_0 - SSR_tot)/SSR_tot > 0.01 % Acceptable condition to consider the best fit reached.

        SSR_tot_0 = SSR_tot;

        % Tables for future SSR values are generated corresponding to the
        % parameters that will be tested for each iteration of the next for loops.

        SSR_neut_table = NaN(N_res,N_res); % SSR on the neutral data set
        SSR_alk_table = NaN(N_res,N_res); % SSR on the alkalie data set
        SSR_tot_table = NaN(N_res,1); % SSR on the total data set

        % The minimal SSR the subset of neutral pH is first investigated by
        % varying k_nat and teta_nat on the predefined tables.
        
        for i = 1:N_res
            for j = 1:N_res
                
                % The fitted parameters are initialized: the ones not initialized here are equal to the initially set values (i.e. lab derived values) 
                k_nat = k_nat_table(i);
                teta_nat = teta_nat_table(j);

                coli_fit = []; % matrix which will gather all the measured (but randomized) values of E. coli cell counts, the fitted values of E. coli cell counts, and an identifier of the corresponding experiment 
                for i_model = 1:nExp % In this loop, all the existing values of measured (but randomized) as well as fitted values of E. coli cell counts are gathered in the order of experiments and time.
                    coli_model_MC{i_model,i_MC} = NaN(nDataFile{i_model},1);

                    k = 1; % This very small loop is necessary as E. coli cell count was not always available from the experiment begining
                    while isnan(coli_data_MC{i_model,i_MC}(k)) && k < nDataFile{i_model}
                        k = k + 1;
                    end
                    coli_model_MC{i_model,i_MC}(k) = coli_data_MC{i_model,i_MC}(k); % First available value is found and the first model is set equal to it. Importantly, this 
                    % value is not added to the coli_fit matrix meaning it is not used in the computation of the SSR. 
                    
                    % The modelled values are here computed based on fitted
                    % parameters values and first order decay law in batch
                    % reactors using Euler method (see manuscript)
                    for m = k + 1:nDataFile{i_model}
                        k_nat_inst = k_nat*teta_nat^(temp_data_MC{i_model,i_MC}(m-1)-20);
                        k_pH_inst = k_pH*teta_pH^(temp_data_MC{i_model,i_MC}(m-1)-20)*10^(pH_data_MC{i_model,i_MC}(m-1)-14);
                        k_sun_inst = alpha_sun*sun_data_MC{i_model,i_MC}(m-1)/(sigma_MC(i_MC,sigmaType(i_model))*d_MC(i_MC))*(1-exp(-sigma_MC(i_MC,sigmaType(i_model))*d_MC(i_MC)));

                        coli_model_MC{i_model,i_MC}(m) = max(1,coli_model_MC{i_model,i_MC}(m-1)*(1 - (k_nat_inst + k_pH_inst + k_sun_inst)*(time_data_MC{i_model,i_MC}(m)-time_data_MC{i_model,i_MC}(m-1))));
                    end

                    index_coli = find(~isnan(coli_data_MC{i_model,i_MC}));
                    coli_fit = [coli_fit ; coli_data_MC{i_model,i_MC}(index_coli), coli_model_MC{i_model,i_MC}(index_coli), i_model*ones(length(index_coli),1)];
                end

                % Sub data sets for neutral and alkaline pH are identified
                % and created.
                coli_fit_neut = [];
                coli_fit_alk = [];
                for i_fit = 1:nExp
                    if highPH(i_fit) == 0
                        row = find(coli_fit(:,3) == i_fit);
                        row = row(2:end);
                        coli_fit_neut = [coli_fit_neut ; coli_fit(row,1:2)];
                    else        
                        row = find(coli_fit(:,3) == i_fit);
                        row = row(2:end);
                        coli_fit_alk = [coli_fit_alk ; coli_fit(row,1:2)];
                    end
                end

                % Log cell count are log-transformed
                log_coli_fit_neut = log10(coli_fit_neut);
                log_coli_fit_alk = log10(coli_fit_alk);
                
                % SSR for neutral and alkaline pH are computed.
                sum_sq_res_neut = sum((log_coli_fit_neut(:,1)-log_coli_fit_neut(:,2)).^2);
                mean_obs_neut = mean(log_coli_fit_neut(:,1));
                sum_sq_tot_neut = sum((log_coli_fit_neut(:,1)-mean_obs_neut).^2);

                sum_sq_res_alk = sum((log_coli_fit_alk(:,1)-log_coli_fit_alk(:,2)).^2);
                mean_obs_alk = mean(log_coli_fit_alk(:,1));
                sum_sq_tot_alk = sum((log_coli_fit_alk(:,1)-mean_obs_alk).^2);

                SSR_neut_table(i,j) = sum_sq_res_neut;

            end
        end
        
        % Having computed all SSR on the neutral data set, the minimum is
        % found and corresponding fitted parameters are retrieved
        SSR_neut_table = real(SSR_neut_table); % In some diverging cases, complex numbers were obtained creating error in the algorithm. Selecting on real parts of the number solved this issue while still confidently selecting an existing minimum.

        [row,col] = find(SSR_neut_table == min(SSR_neut_table,[],'all'));

        if numel(row) == 1 && numel(col) == 1 % This condition is added to make sure a minimum exists.
            % k_nat and teta_nat triggering the minimum SSR are retrieved
            % and  replacing the "current" fitted parameter final value
            k_nat = k_nat_table(row);
            teta_nat = teta_nat_table(col);

            % The same algorithm is now repeated for the subset of high pH
            % when varying k_pH and teta_pH on the predefined tables.
            for i = 1:N_res
                for j = 1:N_res

                    k_pH = k_pH_table(i);
                    teta_pH = teta_pH_table(j);

                    coli_fit = [];
                    for i_model = 1:nExp
                        coli_model_MC{i_model,i_MC} = NaN(nDataFile{i_model},1);

                        k = 1;
                        while isnan(coli_data_MC{i_model,i_MC}(k)) && k < nDataFile{i_model}
                            k = k + 1;
                        end
                        coli_model_MC{i_model,i_MC}(k) = coli_data_MC{i_model,i_MC}(k);

                        for m = k + 1:nDataFile{i_model}
                            k_nat_inst = k_nat*teta_nat.^(temp_data_MC{i_model,i_MC}(m-1)-20);
                            k_pH_inst = k_pH*teta_pH.^(temp_data_MC{i_model,i_MC}(m-1)-20)*10^(pH_data_MC{i_model,i_MC}(m-1)-14);
                            k_sun_inst = alpha_sun*sun_data_MC{i_model,i_MC}(m-1)/(sigma_MC(i_MC,sigmaType(i_model))*d_MC(i_MC))*(1-exp(-sigma_MC(i_MC,sigmaType(i_model))*d_MC(i_MC)));

                            coli_model_MC{i_model,i_MC}(m) = max(1,coli_model_MC{i_model,i_MC}(m-1)*(1 - (k_nat_inst + k_pH_inst + k_sun_inst)*(time_data_MC{i_model,i_MC}(m)-time_data_MC{i_model,i_MC}(m-1))));
                        end

                        index_coli = find(~isnan(coli_data_MC{i_model,i_MC})); 
                        coli_fit = [coli_fit ; coli_data_MC{i_model,i_MC}(index_coli), coli_model_MC{i_model,i_MC}(index_coli), i_model*ones(length(index_coli),1)];
                    end

                    coli_fit_neut = [];
                    coli_fit_alk = [];
                    for i_fit = 1:nExp
                        if highPH(i_fit) == 0
                            row = find(coli_fit(:,3) == i_fit);
                            row = row(2:end);
                            coli_fit_neut = [coli_fit_neut ; coli_fit(row,1:2)];
                        else        
                            row = find(coli_fit(:,3) == i_fit);
                            row = row(2:end);
                            coli_fit_alk = [coli_fit_alk ; coli_fit(row,1:2)];
                        end
                    end

                    log_coli_fit_neut = log10(coli_fit_neut);
                    log_coli_fit_alk = log10(coli_fit_alk);

                    sum_sq_res_neut = sum((log_coli_fit_neut(:,1)-log_coli_fit_neut(:,2)).^2);
                    mean_obs_neut = mean(log_coli_fit_neut(:,1));
                    sum_sq_tot_neut = sum((log_coli_fit_neut(:,1)-mean_obs_neut).^2);

                    sum_sq_res_alk = sum((log_coli_fit_alk(:,1)-log_coli_fit_alk(:,2)).^2);
                    mean_obs_alk = mean(log_coli_fit_alk(:,1));
                    sum_sq_tot_alk = sum((log_coli_fit_alk(:,1)-mean_obs_alk).^2);

                    SSR_alk_table(i,j) = sum_sq_res_alk;
                    SSR_neut_table(i,j) = sum_sq_res_neut;

                end
            end

            SSR_alk_table = real(SSR_alk_table);
            SSR_neut_table = real(SSR_neut_table);

            [row,col] = find(SSR_alk_table == min(SSR_alk_table,[],'all'));

            if numel(row) == 1 && numel(col) == 1 
                k_pH = k_pH_table(row);
                teta_pH = teta_pH_table(col);
                SSR_alk = min(SSR_alk_table,[],'all');
                SSR_neut = SSR_neut_table(row,col);

                % The same algorithm is now repeated for the whole dataset
                % by varying only alpha sun

                for i = 1:N_res

                    alpha_sun = alpha_sun_table(i);

                    coli_fit = [];
                    for i_model = 1:nExp
                        coli_model_MC{i_model,i_MC} = NaN(nDataFile{i_model},1);

                        k = 1;
                        while isnan(coli_data_MC{i_model,i_MC}(k)) && k < nDataFile{i_model}
                            k = k + 1;
                        end
                        coli_model_MC{i_model,i_MC}(k) = coli_data_MC{i_model,i_MC}(k);

                        for m = k + 1:nDataFile{i_model}
                            k_nat_inst = k_nat*teta_nat.^(temp_data_MC{i_model,i_MC}(m-1)-20);
                            k_pH_inst = k_pH*teta_pH.^(temp_data_MC{i_model,i_MC}(m-1)-20)*10^(pH_data_MC{i_model,i_MC}(m-1)-14);
                            k_sun_inst = alpha_sun*sun_data_MC{i_model,i_MC}(m-1)/(sigma_MC(i_MC,sigmaType(i_model))*d_MC(i_MC))*(1-exp(-sigma_MC(i_MC,sigmaType(i_model))*d_MC(i_MC)));

                            coli_model_MC{i_model,i_MC}(m) = max(1,coli_model_MC{i_model,i_MC}(m-1)*(1 - (k_nat_inst + k_pH_inst + k_sun_inst)*(time_data_MC{i_model,i_MC}(m)-time_data_MC{i_model,i_MC}(m-1))));
                        end

                        index_coli = find(~isnan(coli_data_MC{i_model,i_MC})); 
                        coli_fit = [coli_fit ; coli_data_MC{i_model,i_MC}(index_coli), coli_model_MC{i_model,i_MC}(index_coli), i_model*ones(length(index_coli),1)];
                    end

                    coli_fit_neut = [];
                    coli_fit_alk = [];
                    for i_fit = 1:nExp
                        if highPH(i_fit) == 0
                            row = find(coli_fit(:,3) == i_fit);
                            row = row(2:end);
                            coli_fit_neut = [coli_fit_neut ; coli_fit(row,1:2)];
                        else        
                            row = find(coli_fit(:,3) == i_fit);
                            row = row(2:end);
                            coli_fit_alk = [coli_fit_alk ; coli_fit(row,1:2)];
                        end
                    end

                    log_coli_fit_neut = log10(coli_fit_neut);
                    log_coli_fit_alk = log10(coli_fit_alk);

                    sum_sq_res_tot = sum((log_coli_fit_neut(:,1)-log_coli_fit_neut(:,2)).^2) + sum((log_coli_fit_alk(:,1)-log_coli_fit_alk(:,2)).^2);
                    SSR_tot_table(i) = sum_sq_res_tot;

                end

                SSR_tot_table = real(SSR_tot_table);
                [row,~] = find(SSR_tot_table == min(SSR_tot_table,[],'all'));

                if numel(row) == 1
                    alpha_sun = alpha_sun_table(row);
                    SSR_tot = SSR_tot_table(row);
                else
                    display(strcat('resolution too high - no single minima was found at loop',num2str(i_MC)))
                    break
                end
            else
                display(strcat('resolution too high - no single minima was found at loop',num2str(i_MC)))
                break
            end
        else
            display(strcat('resolution too high - no single minima was found at loop',num2str(i_MC)))
            break
        end
    end
    % The first prior optimization on a large mesh grid is finishing here

    % A second optimization round is performed in order to precise the
    % local minima of the the SSR obtained by reducing drastically the
    % resolution of calculation in the vicinity of the values
    % pre-established until the fitting parameters don't move by more than
    % 1%, gradually reducing the resolution.

    k_nat_0 = 0; % these values (unrealistic) are set to make sure the while loop is not initially ignored
    teta_nat_0 = 1;
    k_pH_0 = 0;
    teta_pH_0 = 1;
    alpha_sun_0 = 0;

    N_res = N_res_opt;
    N_res_temoin = N_res*ones(N_MC);

    while max([abs(k_nat-k_nat_0)/k_nat,abs(teta_nat-teta_nat_0)/teta_nat,abs(k_pH-k_pH_0)/k_pH,abs(teta_pH-teta_pH_0)/teta_pH,abs(alpha_sun-alpha_sun_0)/alpha_sun]) > 0.001

        % The "current" optimized values for fitted parameters are assigned
        % as central values of the table of fitted parameters to be tested.
        k_nat_0 = k_nat;
        teta_nat_0 = teta_nat;
        k_pH_0 = k_pH;
        teta_pH_0 = teta_pH;
        alpha_sun_0 = alpha_sun;

        N_res = N_res_temoin(i_MC);

        % Tables of values to be tested are created by encompassing the
        % preestablished value by two times the definition used in the previous
        % step (the definition only increases from here. This step does put us
        % at risk of having missed the absolute minima but instead to focus on
        % a local minima. Using a relatively small mesh in the first step (i.e.
        % N_res >= 50) should prevent such issues.

        min_k_nat_0 = max(0,k_nat - 2*(max_k_nat-min_k_nat)/N_res);
        max_k_nat_0 = k_nat + 2*(max_k_nat-min_k_nat)/N_res;
        min_k_nat = min_k_nat_0;
        max_k_nat = max_k_nat_0;

        min_teta_nat_0 = max(1,teta_nat - 2*(max_teta_nat-min_teta_nat)/N_res);
        max_teta_nat_0 = teta_nat + 2*(max_teta_nat-min_teta_nat)/N_res;
        min_teta_nat = min_teta_nat_0;
        max_teta_nat = max_teta_nat_0;

        min_k_pH_0 = max(0,k_pH - 2*(max_k_pH-min_k_pH)/N_res);
        max_k_pH_0 = k_pH + 2*(max_k_pH-min_k_pH)/N_res;
        min_k_pH = min_k_pH_0;
        max_k_pH = max_k_pH_0;


        min_teta_pH_0 = max(1,teta_pH - 2*(max_teta_pH-min_teta_pH)/N_res);
        max_teta_pH_0 = teta_pH + 2*(max_teta_pH-min_teta_pH)/N_res;
        min_teta_pH = min_teta_pH_0;
        max_teta_pH = max_teta_pH_0;


        min_alpha_sun_0 = max(0,alpha_sun - 2*(max_alpha_sun-min_alpha_sun)/N_res);
        max_alpha_sun_0 = alpha_sun + 2*(max_alpha_sun-min_alpha_sun)/N_res;
        min_alpha_sun = min_alpha_sun_0;
        max_alpha_sun = max_alpha_sun_0;


        k_nat_table = min_k_nat:(max_k_nat-min_k_nat)/N_res:max_k_nat;
        teta_nat_table = min_teta_nat:(max_teta_nat-min_teta_nat)/N_res:max_teta_nat;
        k_pH_table = min_k_pH:(max_k_pH-min_k_pH)/N_res:max_k_pH;
        teta_pH_table = min_teta_pH:(max_teta_pH-min_teta_pH)/N_res:max_teta_pH;
        alpha_sun_table = min_alpha_sun:(max_alpha_sun-min_alpha_sun)/N_res:max_alpha_sun;

        % The same fitting routine is now repeated

        SSR_neut_table = NaN(N_res,N_res);
        SSR_alk_table = NaN(N_res,N_res);
        SSR_tot_table = NaN(N_res,1);

        for i = 1:N_res
            for j = 1:N_res

                k_nat = k_nat_table(i);
                teta_nat = teta_nat_table(j);

                coli_fit = [];
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

                        coli_model_MC{i_model,i_MC}(m) = max(1,coli_model_MC{i_model,i_MC}(m-1)*(1 - (k_nat_inst + k_pH_inst + k_sun_inst)*(time_data_MC{i_model,i_MC}(m)-time_data_MC{i_model,i_MC}(m-1))));
                    end

                    index_coli = find(~isnan(coli_data_MC{i_model,i_MC})); 
                    coli_fit = [coli_fit ; coli_data_MC{i_model,i_MC}(index_coli), coli_model_MC{i_model,i_MC}(index_coli), i_model*ones(length(index_coli),1)];
                end

                coli_fit_neut = [];
                coli_fit_alk = [];
                for i_fit = 1:nExp
                    if highPH(i_fit) == 0
                        row = find(coli_fit(:,3) == i_fit);
                        row = row(2:end);
                        coli_fit_neut = [coli_fit_neut ; coli_fit(row,1:2)];
                    else        
                        row = find(coli_fit(:,3) == i_fit);
                        row = row(2:end);
                        coli_fit_alk = [coli_fit_alk ; coli_fit(row,1:2)];
                    end
                end

                log_coli_fit_neut = log10(coli_fit_neut);
                log_coli_fit_alk = log10(coli_fit_alk);

                sum_sq_res_neut = sum((log_coli_fit_neut(:,1)-log_coli_fit_neut(:,2)).^2);
                mean_obs_neut = mean(log_coli_fit_neut(:,1));
                sum_sq_tot_neut = sum((log_coli_fit_neut(:,1)-mean_obs_neut).^2);

                sum_sq_res_alk = sum((log_coli_fit_alk(:,1)-log_coli_fit_alk(:,2)).^2);
                mean_obs_alk = mean(log_coli_fit_alk(:,1));
                sum_sq_tot_alk = sum((log_coli_fit_alk(:,1)-mean_obs_alk).^2);

                SSR_neut_table(i,j) = sum_sq_res_neut;
            end
        end

        SSR_neut_table = real(SSR_neut_table);
        [row,col] = find(SSR_neut_table == min(SSR_neut_table,[],'all'));

        if numel(row) == 1 && numel(col) == 1
            k_nat = k_nat_table(row);
            teta_nat = teta_nat_table(col);

            for i = 1:N_res
                for j = 1:N_res

                    k_pH = k_pH_table(i);
                    teta_pH = teta_pH_table(j);

                    coli_fit = [];
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
                            coli_model_MC{i_model,i_MC}(m) = max(1,coli_model_MC{i_model,i_MC}(m-1)*(1 - (k_nat_inst + k_pH_inst + k_sun_inst)*(time_data_MC{i_model,i_MC}(m)-time_data_MC{i_model,i_MC}(m-1))));
                        end

                        index_coli = find(~isnan(coli_data_MC{i_model,i_MC})); 
                        coli_fit = [coli_fit ; coli_data_MC{i_model,i_MC}(index_coli), coli_model_MC{i_model,i_MC}(index_coli), i_model*ones(length(index_coli),1)];
                    end

                    coli_fit_neut = [];
                    coli_fit_alk = [];
                    for i_fit = 1:nExp
                        if highPH(i_fit) == 0
                            row = find(coli_fit(:,3) == i_fit);
                            row = row(2:end);
                            coli_fit_neut = [coli_fit_neut ; coli_fit(row,1:2)];
                        else        
                            row = find(coli_fit(:,3) == i_fit);
                            row = row(2:end);
                            coli_fit_alk = [coli_fit_alk ; coli_fit(row,1:2)];
                        end
                    end

                    log_coli_fit_neut = log10(coli_fit_neut);
                    log_coli_fit_alk = log10(coli_fit_alk);

                    sum_sq_res_neut = sum((log_coli_fit_neut(:,1)-log_coli_fit_neut(:,2)).^2);
                    mean_obs_neut = mean(log_coli_fit_neut(:,1));
                    sum_sq_tot_neut = sum((log_coli_fit_neut(:,1)-mean_obs_neut).^2);

                    sum_sq_res_alk = sum((log_coli_fit_alk(:,1)-log_coli_fit_alk(:,2)).^2);
                    mean_obs_alk = mean(log_coli_fit_alk(:,1));
                    sum_sq_tot_alk = sum((log_coli_fit_alk(:,1)-mean_obs_alk).^2);

                    SSR_alk_table(i,j) = sum_sq_res_alk;
                    SSR_neut_table(i,j) = sum_sq_res_neut;
                end
            end

            SSR_alk_table = real(SSR_alk_table);
            SSR_neut_table = real(SSR_neut_table);

            [row,col] = find(SSR_alk_table == min(SSR_alk_table,[],'all'));

            if numel(row) == 1 && numel(col) == 1
                k_pH = k_pH_table(row);
                teta_pH = teta_pH_table(col);
                SSR_alk = min(SSR_alk_table,[],'all');
                SSR_neut = SSR_neut_table(row,col);
                
                for i = 1:N_res

                    alpha_sun = alpha_sun_table(i);

                    coli_fit = [];
                    for i_model = 1:nExp
                        coli_model_MC{i_model,i_MC} = NaN(nDataFile{i_model},1);

                        k = 1;
                        while isnan(coli_data_MC{i_model,i_MC}(k)) && k < nDataFile{i_model}
                            k = k + 1;
                        end
                        coli_model_MC{i_model,i_MC}(k) = coli_data_MC{i_model,i_MC}(k);

                        for m = k + 1:nDataFile{i_model}
                            k_nat_inst = k_nat*teta_nat.^(temp_data_MC{i_model,i_MC}(m-1)-20);
                            k_pH_inst = k_pH*teta_pH.^(temp_data_MC{i_model,i_MC}(m-1)-20)*10^(pH_data_MC{i_model,i_MC}(m-1)-14);
                            k_sun_inst = alpha_sun*sun_data_MC{i_model,i_MC}(m-1)/(sigma_MC(i_MC,sigmaType(i_model))*d_MC(i_MC))*(1-exp(-sigma_MC(i_MC,sigmaType(i_model))*d_MC(i_MC)));

                            coli_model_MC{i_model,i_MC}(m) = max(1,coli_model_MC{i_model,i_MC}(m-1)*(1 - (k_nat_inst + k_pH_inst + k_sun_inst)*(time_data_MC{i_model,i_MC}(m)-time_data_MC{i_model,i_MC}(m-1))));
                        end

                        index_coli = find(~isnan(coli_data_MC{i_model,i_MC})); 
                        coli_fit = [coli_fit ; coli_data_MC{i_model,i_MC}(index_coli), coli_model_MC{i_model,i_MC}(index_coli), i_model*ones(length(index_coli),1)];
                    end

                    coli_fit_neut = [];
                    coli_fit_alk = [];
                    for i_fit = 1:nExp
                        if highPH(i_fit) == 0
                            row = find(coli_fit(:,3) == i_fit);
                            row = row(2:end);
                            coli_fit_neut = [coli_fit_neut ; coli_fit(row,1:2)];
                        else        
                            row = find(coli_fit(:,3) == i_fit);
                            row = row(2:end);
                            coli_fit_alk = [coli_fit_alk ; coli_fit(row,1:2)];
                        end
                    end

                    log_coli_fit_neut = log10(coli_fit_neut);
                    log_coli_fit_alk = log10(coli_fit_alk);

                    log_coli_fit_neut = log10(coli_fit_neut);
                    log_coli_fit_alk = log10(coli_fit_alk);

                    sum_sq_res_neut = sum((log_coli_fit_neut(:,1)-log_coli_fit_neut(:,2)).^2);
                    mean_obs_neut = mean(log_coli_fit_neut(:,1));
                    sum_sq_tot_neut = sum((log_coli_fit_neut(:,1)-mean_obs_neut).^2);

                    sum_sq_res_alk = sum((log_coli_fit_alk(:,1)-log_coli_fit_alk(:,2)).^2);
                    mean_obs_alk = mean(log_coli_fit_alk(:,1));
                    sum_sq_tot_alk = sum((log_coli_fit_alk(:,1)-mean_obs_alk).^2);

                    log_coli_fit_tot = [log_coli_fit_neut;log_coli_fit_alk];
                    
                    sum_sq_res_tot = sum((log_coli_fit_tot(:,1)-log_coli_fit_tot(:,2)).^2);
                    SSR_tot_table(i) = sum_sq_res_tot;

                end

                SSR_tot_table = real(SSR_tot_table);

                [row,~] = find(SSR_tot_table == min(SSR_tot_table,[],'all'));
                if numel(row) == 1
                    alpha_sun = alpha_sun_table(row);
                    SSR_tot = SSR_tot_table(row);
                else
                    N_res_temoin(i_MC) = N_res_temoin(i_MC)*2; % If no optimal value is found, the resolution is reduced and calculation is repeated
                end
            else               
                N_res_temoin(i_MC) = N_res_temoin(i_MC)*2; % itou
            end
        else
            N_res_temoin(i_MC) = N_res_temoin(i_MC)*2; % itou
        end
        if N_res > 500 % This safety was created to avoid infinite loop if an optimized value cannot be found.
            display(strcat('resolution too high - might have entered infinite loop at loop n°',num2str(i_MC)))
            break
        end
    end

    log_coli_fit_tot = [log_coli_fit_neut;log_coli_fit_alk];
    sum_sq_res_tot = sum((log_coli_fit_tot(:,1)-log_coli_fit_tot(:,2)).^2);
    mean_obs_tot = mean(log_coli_fit_tot(:,1));
    sum_sq_tot_tot = sum((log_coli_fit_tot(:,1)-mean_obs_tot).^2);
    R2_tot = 1 - sum_sq_res_tot/sum_sq_tot_tot;
    SSR_tot = sum_sq_res_tot;

    %% Results storage
    SSR_MC(:,i_MC) = SSR_tot;
    R2_MC(:,i_MC) = R2_tot;
    SSR_MC_alk(:,i_MC) = SSR_alk;
    SSR_MC_neut(:,i_MC) = SSR_neut;

    coli_MC_logcounts(:,i_MC) = log_coli_fit_tot(:,1);
    coli_MC_logmodel(:,i_MC) = log_coli_fit_tot(:,2);

    k_pH_MC(i_MC) = k_pH;
    teta_pH_MC(i_MC) = teta_pH;
    k_nat_MC(i_MC) = k_nat;
    teta_nat_MC(i_MC) = teta_nat;
    alpha_sun_MC(i_MC) = alpha_sun;

end




