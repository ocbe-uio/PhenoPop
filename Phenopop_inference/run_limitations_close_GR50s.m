function [] = run_limitations_close_GR50s(slurm_array_task_id)

    % Inference for different cases
    % Noise 50, 1000*Npop optimations. Elbow method used for model selection.
    % Make folders named 1,2,3 for saving different numbers of populations
    % Make folders named case1,case2,case3 within those for the different cases
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These are copied, not coded: 
    a1 = 0.03;
    a2 = 0.03;
    a3 = 0.03;
    a4 = 0.03;
    a5 = 0.03;
    a6 = 0.03;
    a7 = 0.03;
    b1 = 0.3;
    b2 = 0.4;
    b3 = 0.4;
    b4 = 0.4;
    b5 = 0.4;
    b6 = 0.4;
    b7 = 0.4;
    E_values = logspace(log10(1.0e-4), log10(1.0e-2), 7)
    E1 = E_values(1);
    E2 = E_values(2);
    E3 = E_values(3);
    E4 = E_values(4);
    E5 = E_values(5);
    E6 = E_values(6);
    E7 = E_values(7);
    n = 3;
    
    cell_line_1 = [a1, b1, E1, n]
    cell_line_2 = [a2, b2, E2, n]
    cell_line_3 = [a3, b3, E3, n]
    cell_line_4 = [a4, b4, E4, n]
    cell_line_5 = [a5, b5, E5, n]
    cell_line_6 = [a6, b6, E6, n]
    cell_line_7 = [a7, b7, E7, n]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    base_num_optim = 1000 % set globally and increased with no_populations to look for. 
    caseno = slurm_array_task_id
    
    % Difference in GR50 changes
    switch slurm_array_task_id
        % 50 / 50 mixture
        case 1
            N_true_populations = 2
            Noise = 50
            MixParams = [0.5]; % Length k-1
            Params = [cell_line_1 cell_line_2];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
        case 2
            N_true_populations = 2
            Noise = 50
            MixParams = [0.5]; % Length k-1
            Params = [cell_line_1 cell_line_3];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
        case 3
            N_true_populations = 2
            Noise = 50
            MixParams = [0.5]; % Length k-1
            Params = [cell_line_1 cell_line_4];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
        % 90 / 10 mixture
        case 4
            N_true_populations = 2
            Noise = 50
            MixParams = [0.9]; % Length k-1
            Params = [cell_line_1 cell_line_2];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
        case 5
            N_true_populations = 2
            Noise = 50
            MixParams = [0.9]; % Length k-1
            Params = [cell_line_1 cell_line_3];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
        case 6
            N_true_populations = 2
            Noise = 50
            MixParams = [0.9]; % Length k-1
            Params = [cell_line_1 cell_line_4];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
        % 95 / 5 mixture
        case 7
            N_true_populations = 2
            Noise = 50
            MixParams = [0.95]; % Length k-1
            Params = [cell_line_1 cell_line_2];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
        case 8
            N_true_populations = 2
            Noise = 50
            MixParams = [0.95]; % Length k-1
            Params = [cell_line_1 cell_line_3];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
        case 9
            N_true_populations = 2
            Noise = 50
            MixParams = [0.95]; % Length k-1
            Params = [cell_line_1 cell_line_4];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
        % 99 / 1 mixture
        case 10
            N_true_populations = 2
            Noise = 50
            MixParams = [0.99]; % Length k-1
            Params = [cell_line_1 cell_line_2];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
        case 11
            N_true_populations = 2
            Noise = 50
            MixParams = [0.99]; % Length k-1
            Params = [cell_line_1 cell_line_3];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
        case 12
            N_true_populations = 2
            Noise = 50
            MixParams = [0.99]; % Length k-1
            Params = [cell_line_1 cell_line_4];
            articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    end
    
    % Old version
    % Freaky cell lines. We plot 2 instead of 3 populations if that's the case
    
    % a) Just too much noise, 1 population
    
    % c) Narrow cell lines, GR50 close to ea.other. 2 populations, Noise 200
    
    % d) Narrow cell lines, GR50 close to ea.other. 3 populations, Noise 200
    
    % e) Small mixture proportions. Add a small third mixture and have the model select 2 populations and talk about model selection?
    
    %switch slurm_array_task_id
    %case 1 % Too much noise
    %    N_true_populations = 1
    %    Noise = 1000
    %    MixParams = []; % Length k-1
    %    Params = [cell_line_3];
    %    articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    %case 2 % Too close cell lines at some noise level
    %    N_true_populations = 2
    %    Noise = 200
    %    MixParams = [0.5]; % Length k-1
    %    Params = [cell_line_3 cell_line_4];
    %    articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    %case 3 % Too close cell lines at some noise level
    %    N_true_populations = 3
    %    Noise = 200
    %    MixParams = [1/3 1/3]; % Length k-1
    %    Params = [cell_line_3 cell_line_4 cell_line_5];
    %    articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    %case 4 % Small mixture parameter
    %    N_true_populations = 2
    %    Noise = 200
    %    MixParams = [0.95]; % Length k-1
    %    Params = [cell_line_3 cell_line_7];
    %    articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    %case 5 % Small mixture parameter
    %    N_true_populations = 3
    %    Noise = 200
    %    MixParams = [0.05, 0.90]; % Length k-1
    %    Params = [cell_line_1 cell_line_3 cell_line_7];
    %    articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, base_num_optim, caseno, true)
    %end
    %end
    