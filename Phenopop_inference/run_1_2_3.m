% Inference for different cases
% Noise 10, 200 optims. Elbow method used for model selection.
% Make folders named 1,2,3 for saving different numbers of populations
% Make folders named case1,case2,case3 within those for the different cases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are copied, not coded: 
base_rate_sensitive = 0.03;
base_rate_middle = 0.03;
base_rate_resistant = 0.03;
b1 = 0.3;
b2 = 0.4;
b3 = 0.5;
E_sensitive = 1.0e-4;
E_middle = 1.0e-2
E_resistant = 1.0e-1;
n = 3;

sensitive_cell_line = [base_rate_sensitive, b1, E_sensitive, n]
middle_cell_line = [base_rate_middle, b2, E_middle, n]
resistant_cell_line = [base_rate_resistant, b3, E_resistant, n]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_num_optim = 1000 % set globally and increased with no_populations to look for. 

% 1 population
N_true_populations = 1
Noise = 50
MixParams = []; % Length k-1
% Case 1
caseno = 1
Params = [sensitive_cell_line];
articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, false)

% Case 2
caseno = 2
Params = [middle_cell_line];
articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, false)

% Case 3
caseno = 3
Params = [resistant_cell_line];
articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, false)


% 2 populations
N_true_populations = 2
% Case 1
caseno = 1
MixParams = [0.5]; % Length k-1
Params = [sensitive_cell_line middle_cell_line];
articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, false)

% Case 2
caseno = 2
MixParams = [0.7]; % Length k-1
Params = [sensitive_cell_line middle_cell_line];
articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, false)

% Case 3
caseno = 3
MixParams = [0.3]; % Length k-1
Params = [sensitive_cell_line middle_cell_line];
articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, false)


% 3 populations 
N_true_populations = 3
Params = [sensitive_cell_line middle_cell_line resistant_cell_line];
% Case 1
caseno = 1
MixParams = [1/3 1/3]; % Length k-1
articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, false)

% Case 2
caseno = 2
MixParams = [0.4 0.3]; % Length k-1
articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, false)

% Case 3
caseno = 3
MixParams = [0.6 0.2]; % Length k-1
articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, false)
