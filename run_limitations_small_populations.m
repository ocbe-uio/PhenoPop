function [] = run_limitations_small_populations(slurm_array_task_id)

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
    Noise = 50
    
    switch slurm_array_task_id
    % 2 populations
    case 1
        % Case 1
        N_true_populations = 2
        MixParams = [0.70]; % Length k-1
        caseno = 1
        Params = [sensitive_cell_line middle_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    case 2
        % Case 2
        N_true_populations = 2
        MixParams = [0.80]; % Length k-1
        caseno = 2
        Params = [sensitive_cell_line middle_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)

    case 3
        % Case 3
        N_true_populations = 2
        MixParams = [0.90]; % Length k-1
        caseno = 3
        Params = [sensitive_cell_line middle_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    
    case 4
        % Case 1
        N_true_populations = 2
        MixParams = [0.95]; % Length k-1
        caseno = 1
        Params = [sensitive_cell_line middle_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    case 5
        % Case 2
        N_true_populations = 2
        MixParams = [0.97]; % Length k-1
        caseno = 2
        Params = [sensitive_cell_line middle_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)

    case 6
        % Case 3
        N_true_populations = 2
        MixParams = [0.99]; % Length k-1
        caseno = 3
        Params = [sensitive_cell_line middle_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)

    % 3 populations 
    % Two big, 1 small
    case 7
        % Case 1
        N_true_populations = 3
        MixParams = [0.50 0.45]; % Length k-1
        caseno = 1
        Params = [sensitive_cell_line middle_cell_line resistant_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    case 8
        % Case 2
        N_true_populations = 3
        MixParams = [0.50 0.47]; % Length k-1
        caseno = 2
        Params = [sensitive_cell_line middle_cell_line resistant_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    
    case 9
        % Case 3
        N_true_populations = 3
        MixParams = [0.50 0.49]; % Length k-1
        caseno = 3
        Params = [sensitive_cell_line middle_cell_line resistant_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    
    % 1 big, two small
    case 10
        % Case 1
        N_true_populations = 3
        MixParams = [0.8 0.1]; % Length k-1
        caseno = 1
        Params = [sensitive_cell_line middle_cell_line resistant_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    
    case 11
        % Case 2
        N_true_populations = 3
        MixParams = [0.9 0.05]; % Length k-1
        caseno = 2
        Params = [sensitive_cell_line middle_cell_line resistant_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    
    case 12
        % Case 3
        N_true_populations = 3
        MixParams = [0.95 0.02]; % Length k-1
        caseno = 3
        Params = [sensitive_cell_line middle_cell_line resistant_cell_line];
        articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    
    end

%% 1 population
%case 1
%    % Case 1
%    N_true_populations = 2
%    MixParams = [0.7]; % Length k-1
%    caseno = 1
%    Params = [sensitive_cell_line middle_cell_line];
%    articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
%case 2
%    % Case 2
%    N_true_populations = 2
%    MixParams = [0.8]; % Length k-1
%    caseno = 2
%    Params = [sensitive_cell_line middle_cell_line];
%    articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
%case 3
%    % Case 3
%    N_true_populations = 2
%    MixParams = [0.9]; % Length k-1
%    caseno = 3
%    Params = [sensitive_cell_line middle_cell_line];
%    articleplot_2(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
