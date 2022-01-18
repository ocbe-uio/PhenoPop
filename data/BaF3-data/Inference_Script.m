% INFER SCRIPT
% 
% In this script, we infer 0,1,2,3,and 4 population models on Dagim's monoclonal
% and mixture data, MONOCLONAL_DATA  and Mixture_Data. Inference results are saved in
% Results_Monoclonal_max4pop_infer and Results_Mixture_max4pop_infer. Inferred models will be in form:
%  [ alpha_1, b_1, E_1, n_1, p_1,alpha_2, b_2, E_2, n_2, p_2, ... alpha_k, b_k, E_k, n_k, sigH, sigL ]
%
% Run time ~ 6 hours
% 
% I want to change...
% 
% 1) number of subpopulation models. 
%       Modify max_no_population and total_pop_models
%
% 2) Optimization Parameters.  
%       - upper bounds for nonnull model: alpha_ub, b_ub, E_ub, n_ub, sigma_ub, mix_ub
%       - upper bounds for null model: ub_null 
%       - lower_mixture_limit
%       - high and low variance threshold indices
%
% 3) Make script run faster. 
%       Decrease num_optim
%
% This file depends on functions
% vectorized_inference_two_noise_levels
% compute_bic
% compute_aic

%% Initializations

max_no_populations=4;
total_pop_models=5;    % models 0 to max_no_population

Conc = 10^(6)*[0  31.25*10^(-9) 62.5*10^(-9) 125*10^(-9) 250*10^(-9) 375*10^(-9) 500*10^(-9) 1.25*10^(-6) 2.5*10^(-6) 3.75*10^(-6) 5*10^(-6) ];
num_optim = 1000;

%set mixture parameter lower limit
lower_mixture_limit = 0.10;

%set upper bound optimization parameters; use for pop models with
% number of subpops > 0 
alpha_ub = 0.1;
b_ub = 1;
E_ub = 50;
n_ub = 20;
mix_ub=0.5;
sig_ub = 2500;


%set upper bound optimization parameters; use for pop model with 
% number of subpop = 0
ub_null = [10000; 5000];
lb_null = zeros(length(ub_null), 1);

%High Low variance thresholds
%high variance thresholds
%MONOCLONAL thresholds
ConcT_R250 = 6;     ConcT_R500 = 4; ConcT_S500= 3; ConcT_S1000= 3;
TimeT_R250 = 14;     TimeT_R500 = 14;   TimeT_S500 = 14;   TimeT_S1000 = 14;

%MIXTURE thresholds
ConcT = 6;      
TimeT = 14;  

%% Infer upto 4 population models on Mixture Data
load Mixture_Data %- use this to generate
%RESULTS_Mixture_Data_no_pops_ & RESULTS_Mixture_Data_one_two_pop_infer_aic_bic_values

% we obtain third mixture parameter by taking 1 - sum(p_i)


[NR, NC, NT] = size(BF_11); %all data sets are the same size.
%initialize structure to save parameters
xfinals=zeros(5*max_no_populations + 1, total_pop_models, 4); % (:,:,4) = [BF11 BF12 BF21 BF41]
fvals = zeros(1, total_pop_models ,4);

%init structures to save bic scores BIC
bic_array=zeros(3,total_pop_models, 4); %(:, :, 4) = [ bf11, bf12, bf21, bf41 ] 
bic_array_num_param_squared=zeros(3,total_pop_models, 4); %(:, :, 4) = [ bf11, bf12, bf21, bf41 ] 
aic_array=zeros(1,total_pop_models,4);

%infer
for ii = 0:max_no_populations  %ii marks the nth subpop
    if ii==0    %infer the null population model
        col=ii+1;
        lb=lb_null;
        ub=ub_null;
        num_params=length(ub);
        
        [fval_BF11, X_BF11 ] = vectorized_inference_k_subpopulations_two_noise_levels(BF_11, [], [], [], ub, lb, num_optim, ii, [], []);
        [fval_BF12, X_BF12 ] = vectorized_inference_k_subpopulations_two_noise_levels(BF_12, [], [], [], ub, lb, num_optim, ii, [], []);
        [fval_BF21, X_BF21 ] = vectorized_inference_k_subpopulations_two_noise_levels(BF_21, [], [], [], ub, lb, num_optim, ii, [], []);
        [fval_BF41, X_BF41 ] = vectorized_inference_k_subpopulations_two_noise_levels(BF_41, [], [], [], ub, lb, num_optim, ii, [], []);
        
        fvals(:, col, 1)=fval_BF11;
        fvals(:, col, 2)=fval_BF12;
        fvals(:, col, 3)=fval_BF21;
        fvals(:, col, 4)=fval_BF41;
        
        xfinals(1:length(X_BF11), col, 1 )=X_BF11;
        xfinals(1:length(X_BF12), col, 2 )=X_BF12;
        xfinals(1:length(X_BF21), col, 3)=X_BF21;
        xfinals(1:length(X_BF41), col, 4)=X_BF41;
        
        
        bic_array(:, col, 1) = compute_bic(BF_11, fval_BF11, num_params);
        bic_array(:, col, 2) = compute_bic(BF_12, fval_BF12, num_params);
        bic_array(:, col, 3) = compute_bic(BF_21, fval_BF21, num_params);
        bic_array(:, col, 4) = compute_bic(BF_41, fval_BF41, num_params);
        
        bic_array_num_param_squared(:, col, 1) = compute_bic(BF_11, fval_BF11, num_params^2);
        bic_array_num_param_squared(:, col, 2) = compute_bic(BF_12, fval_BF12, num_params^2);
        bic_array_num_param_squared(:, col, 3) = compute_bic(BF_21, fval_BF21, num_params^2);
        bic_array_num_param_squared(:, col, 4) = compute_bic(BF_41, fval_BF41, num_params^2);
        
        aic_array(:,col,1) = compute_aic(fval_BF11, num_params);
        aic_array(:,col,2) = compute_aic(fval_BF12, num_params);
        aic_array(:,col,3) = compute_aic(fval_BF21, num_params);
        aic_array(:,col,4) = compute_aic(fval_BF41, num_params);
        
        
    else            %infer model with (ii~=0) subpopulations
        col=ii+1;
        
        % set optimization constraints
        lb= zeros(5*ii + 1,1);
        ub=sig_ub*ones(size(lb));
        
        num_params=length(ub);
        
        ub(1:5:5*ii-4)=alpha_ub; % alpha
        ub(2:5:5*ii-3)=b_ub;  % b
        ub(3:5:5*ii-2)=E_ub; % E
        ub(4:5:5*ii-1)=n_ub; % n
        ub(5:5:5*ii-5)=mix_ub; % Mixture parameter
        
        % Ensure sum inferred mixture to be < 1
        A_inequality=zeros(1,length(ub));
        A_inequality(5:5:5*ii-5) = 1;
        b_inequality = 1 - lower_mixture_limit;
        
        
        [fval_BF11, X_BF11 ] = vectorized_inference_k_subpopulations_two_noise_levels(BF_11, Conc, ConcT, TimeT, ub, lb, num_optim, ii, A_inequality, b_inequality);
        [fval_BF12, X_BF12 ] = vectorized_inference_k_subpopulations_two_noise_levels(BF_12, Conc, ConcT, TimeT, ub, lb, num_optim, ii, A_inequality, b_inequality);
        [fval_BF21, X_BF21 ] = vectorized_inference_k_subpopulations_two_noise_levels(BF_21, Conc, ConcT, TimeT, ub, lb, num_optim, ii, A_inequality, b_inequality);
        [fval_BF41, X_BF41 ] = vectorized_inference_k_subpopulations_two_noise_levels(BF_41, Conc, ConcT, TimeT, ub, lb, num_optim, ii, A_inequality, b_inequality);
        
        fvals(:, col, 1)=fval_BF11;
        fvals(:, col, 2)=fval_BF12;
        fvals(:, col, 3)=fval_BF21;
        fvals(:, col, 4)=fval_BF41;
        
        xfinals(1:length(X_BF11), col, 1 )=X_BF11;
        xfinals(1:length(X_BF12), col, 2 )=X_BF12;
        xfinals(1:length(X_BF21), col, 3)=X_BF21;
        xfinals(1:length(X_BF41), col, 4)=X_BF41;
        
        
        bic_array(:, col, 1) = compute_bic(BF_11, fval_BF11, num_params);
        bic_array(:, col, 2) = compute_bic(BF_12, fval_BF12, num_params);
        bic_array(:, col, 3) = compute_bic(BF_21, fval_BF21, num_params);
        bic_array(:, col, 4) = compute_bic(BF_41, fval_BF41, num_params);
        
        bic_array_num_param_squared(:, col, 1) = compute_bic(BF_11, fval_BF11, num_params^2);
        bic_array_num_param_squared(:, col, 2) = compute_bic(BF_12, fval_BF12, num_params^2);
        bic_array_num_param_squared(:, col, 3) = compute_bic(BF_21, fval_BF21, num_params^2);
        bic_array_num_param_squared(:, col, 4) = compute_bic(BF_41, fval_BF41, num_params^2);
        
        aic_array(:,col,1) = compute_aic(fval_BF11, num_params);
        aic_array(:,col,2) = compute_aic(fval_BF12, num_params);
        aic_array(:,col,3) = compute_aic(fval_BF21, num_params);
        aic_array(:,col,4) = compute_aic(fval_BF41, num_params);
    end
end

[bic_min_value, bic_min_population]= min(bic_array(1,:,:)); %population 2 chosen for all datasets. 
[bic_square_num_params_min_value, bic_square_num_params_min_population]= min(bic_array_num_param_squared(1,:,:)); %population 2 chosen for all datasets. 
[aic_min_value, aic_min_population]= min(aic_array(1,:,:)); %population 2 chosen for all datasets. 

save('Results_Mixture_max4pop_infer', 'bic_array','bic_min_value', 'bic_min_population',...
                                                            'aic_array','aic_min_value', 'aic_min_population',...
                                                            'bic_array_num_param_squared', 'bic_square_num_params_min_value', 'bic_square_num_params_min_population',...
                                                            'xfinals' ,'fvals')
    

%% Infer upto 4 population models on Monoclonal Data
load MONOCLONAL_DATA

[NR, NC, NT] = size(RESISTANT_250_BF);



PLOT_FINAL_FIT=1;
xfinals=zeros(5*max_no_populations + 1, total_pop_models,4); % (:,:,4) = [R250 R500 S500 S1000]
fvals = zeros(1, total_pop_models,4);

%init BIC structures
bic_array=zeros(3,total_pop_models, 4); %(:, :, 4) = [R250 R500 S500 S1000]
bic_array_num_param_squared=zeros(3,total_pop_models, 4); %(:, :, 4) = [R250 R500 S500 S1000]
aic_array=zeros(1,total_pop_models,4);

%infer
for ii = 0: max_no_populations 
    if ii==0    %inferring the zero pop model, infer mean and variance
        col=ii+1;
        
        %set upper and lower bound for optimization. hardcode 
        lb=lb_null;
        ub=ub_null;
        num_params=length(ub);
        
        %infer a model with ii=0 subpopulations for each data set
        [fval_R250, X_R250 ] = vectorized_inference_k_subpopulations_two_noise_levels(RESISTANT_250_BF, [], [], [], ub, lb, num_optim, ii, [], []);
        [fval_R500, X_R500 ] = vectorized_inference_k_subpopulations_two_noise_levels(RESISTANT_500_BF, [], [], [], ub, lb, num_optim, ii, [], []);
        [fval_S500, X_S500 ] = vectorized_inference_k_subpopulations_two_noise_levels(SENSITIVE_1000_BF, [], [], [], ub, lb, num_optim, ii, [], []);
        [fval_S1000, X_S1000 ] = vectorized_inference_k_subpopulations_two_noise_levels(SENSITIVE_500_BF, [], [], [], ub, lb, num_optim, ii, [], []);
        
        %save scores and results into arrays
        fvals(:, col, 1)=fval_R250;
        fvals(:, col, 2)=fval_R500;
        fvals(:, col, 3)=fval_S500;
        fvals(:, col, 4)=fval_S1000;
        
        xfinals(1:length(X_R250), col, 1 )=X_R250;
        xfinals(1:length(X_R500), col, 2 )=X_R500;
        xfinals(1:length(X_S500), col, 3 )=X_S500;
        xfinals(1:length(X_S1000), col, 4 )=X_S1000;
        
        bic_array(:, col, 1) = compute_bic(RESISTANT_250_BF, fval_R250, num_params);
        bic_array(:, col, 2) = compute_bic(RESISTANT_500_BF, fval_R500, num_params);
        bic_array(:, col, 3) = compute_bic(SENSITIVE_500_BF, fval_S500, num_params);
        bic_array(:, col, 4) = compute_bic(SENSITIVE_1000_BF, fval_S1000, num_params);
        
        bic_array_num_param_squared(:, col, 1) = compute_bic(RESISTANT_250_BF, fval_R250, num_params^2);
        bic_array_num_param_squared(:, col, 2) = compute_bic(RESISTANT_500_BF, fval_R500, num_params^2);
        bic_array_num_param_squared(:, col, 3) = compute_bic(SENSITIVE_500_BF, fval_S500, num_params^2);
        bic_array_num_param_squared(:, col, 4) = compute_bic(SENSITIVE_1000_BF, fval_S1000, num_params^2);
        
        aic_array(:,col,1) = compute_aic(fval_R250, num_params);
        aic_array(:,col,2) = compute_aic(fval_R500, num_params);
        aic_array(:,col,3) = compute_aic(fval_S500, num_params);
        aic_array(:,col,4) = compute_aic(fval_S1000, num_params);
        
    else        %infer models with ii ~= 0 subpopulations 
        
        % set optimization constraints
        lb= zeros(5*ii + 1,1);
        ub=sig_ub*ones(size(lb));
        
        num_params=length(ub);
        
        ub(1:5:5*ii-4)=alpha_ub; % alpha
        ub(2:5:5*ii-3)=b_ub;  % b
        ub(3:5:5*ii-2)=E_ub; % E
        ub(4:5:5*ii-1)=n_ub; % n
        ub(5:5:5*ii-5)=mix_ub; % Mixture parameter
        
        
        % Ensure sum inferred mixture to be < 1
        A_inequality=zeros(1,length(ub));
        A_inequality(5:5:5*ii-5) = 1;
        b_inequality = 1 - lower_mixture_limit;
        
        %infer a model with ii subpopulations for each data set
        [fval_R250, X_R250 ] = vectorized_inference_k_subpopulations_two_noise_levels(RESISTANT_250_BF, Conc, ConcT_R250, TimeT_R250, ub, lb, num_optim, ii, A_inequality, b_inequality);
        [fval_R500, X_R500 ] = vectorized_inference_k_subpopulations_two_noise_levels(RESISTANT_500_BF, Conc, ConcT_R500, TimeT_R500, ub, lb, num_optim, ii, A_inequality, b_inequality);
        [fval_S500, X_S500 ] = vectorized_inference_k_subpopulations_two_noise_levels(SENSITIVE_1000_BF, Conc, ConcT_S500, TimeT_S500, ub, lb, num_optim, ii, A_inequality, b_inequality);
        [fval_S1000, X_S1000 ] = vectorized_inference_k_subpopulations_two_noise_levels(SENSITIVE_500_BF, Conc, ConcT_S1000, TimeT_S1000, ub, lb, num_optim, ii, A_inequality, b_inequality);
        
        fvals(:, col, 1)=fval_R250;
        fvals(:, col, 2)=fval_R500;
        fvals(:, col, 3)=fval_S500;
        fvals(:, col, 4)=fval_S1000;
        
        xfinals(1:length(X_R250), col, 1 )=X_R250;
        xfinals(1:length(X_R500), col, 2 )=X_R500;
        xfinals(1:length(X_S500), col, 3)=X_S500;
        xfinals(1:length(X_S1000), col, 4)=X_S1000;
        
        bic_array(:, col, 1) = compute_bic(RESISTANT_250_BF, fval_R250, num_params);
        bic_array(:, col, 2) = compute_bic(RESISTANT_500_BF, fval_R500, num_params);
        bic_array(:, col, 3) = compute_bic(SENSITIVE_500_BF, fval_S500, num_params);
        bic_array(:, col, 4) = compute_bic(SENSITIVE_1000_BF, fval_S1000, num_params);
        
        bic_array_num_param_squared(:, col, 1) = compute_bic(RESISTANT_250_BF, fval_R250, num_params^2);
        bic_array_num_param_squared(:, col, 2) = compute_bic(RESISTANT_500_BF, fval_R500, num_params^2);
        bic_array_num_param_squared(:, col, 3) = compute_bic(SENSITIVE_500_BF, fval_S500, num_params^2);
        bic_array_num_param_squared(:, col, 4) = compute_bic(SENSITIVE_1000_BF, fval_S1000, num_params^2);
        
        aic_array(:,col,1) = compute_aic(fval_R250, num_params);
        aic_array(:,col,2) = compute_aic(fval_R500, num_params);
        aic_array(:,col,3) = compute_aic(fval_S500, num_params);
        aic_array(:,col,4) = compute_aic(fval_S1000, num_params);
    end
end

[bic_min_value, bic_min_population]= min(bic_array(1,:,:)); %population 2 chosen for all datasets. 
[bic_square_num_params_min_value, bic_square_num_params_min_population]= min(bic_array_num_param_squared(1,:,:)); %population 2 chosen for all datasets. 
[aic_min_value, aic_min_population]= min(aic_array(1,:,:)); %population 2 chosen for all datasets. 


save('Results_Monoclonal_max4pop_infer', 'fvals', 'xfinals','bic_array', 'bic_min_value', 'bic_min_population',...
                                                  'aic_array', 'aic_min_value', 'aic_min_population', ...
                                                  'bic_array_num_param_squared', 'bic_square_num_params_min_population', 'bic_square_num_params_min_value')




