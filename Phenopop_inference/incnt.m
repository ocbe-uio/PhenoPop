function out = increase_nt(slurm_array_task_id, savebool)
    % Increase either Nc, Nc or R. 
    % 'Changing ',num2str(change_value), ', others fixed.'

    % This script runs MLE in parallell 
    % slurm_array_task_id determines: MixParam, 
    % ? n, height and width
    slurm_array_task_id = slurm_array_task_id
    
    % Specific to inc_nc_nt_r
    N_repetitions = 30
    base_case = 3
    no_populations = 2

    num_optim = 1000
    seednr_1 = 43
    generation_number = 4
    RUNNING_ON_CLUSTER = true;
    Noise = 100
    lower_limit_mixparam = 0 % Each inferred mixture parameter must be greater than this.
    lower_limit_E_ratio = 1 % The smallest ratio (E_i / E_j) among inferred E parameters must be greater than this.
    USE_E_THRESHOLD = false
    INFER_SIGMA = true % If true, sigma is estimated by MLE. If false, the true value is given.
    USE_TWO_NOISE_LEVELS = true
    if USE_TWO_NOISE_LEVELS %thresholds for sigma
        ConcT = 1.0e-3; %=E_sensitive
        TimeT = 48;
        if ~INFER_SIGMA
            error("If using two noise level, sigma must also be inferred or you must implement an estimate for sigma")
        end
    end
    %PLOT_DATA = true
    %PLOT_FINAL_FIT = true

    %% Define cell lines
    %base_rate_sensitive = 0.03;
    %base_rate_middle = 0.03;
    %base_rate_resistant = 0.03;
    %b1 = 0.3;
    %b2 = 0.4;
    %b3 = 0.5;
    %E_sensitive = 1.0e-4;
    %E_middle = 1.0e-2
    %E_resistant = 1.0e-1;
    %n = 3;
    %
    %sensitive_cell_line = [base_rate_sensitive, b1, E_sensitive, n]
    %middle_cell_line = [base_rate_middle, b2, E_middle, n]
    %resistant_cell_line = [base_rate_resistant, b3, E_resistant, n]
    
    %tlabels = strings(1,N_t);
    %for t_index=1:N_t
    %    tlabels(N_t-t_index+1) = strcat("T=", int2str(Time(t_index)));
    %end
    
    %change_value = 'R'
    %change_value = 'C'
    change_value = 'T'

    if true %generation_number == % 1 or 2
        % Based on Dagim's data. Concentrations around 1e-3 instead of 1e-6 due to numerical issues in find_the_right_concentrations.m
        base_rate_resistant = 0.03;
        base_rate_sensitive = 0.03;
        b_values = [0.3 0.5 0.7 0.9];
        outer_n_values = [1 3 5 7 9];
        outer_MixParams = 0.1:0.1:0.9;
        E_sensitive = 1.0e-3;
        E_factors = [4 16 64 256];

        gridlength = 4;
        N_c_values = 2.^(1:gridlength)+1; %[3 5 9 17];
        N_t_values = 2.^(1:gridlength)+1; %[3 5 9 17];
        R_values = 2.^(1:gridlength)+1; %[3 5 9 17];
    end

    max_N_c = max(N_c_values);
    max_N_t = max(N_t_values);
    max_R = max(R_values);

    gridlength_b = length(b_values);
    gridlength_E = length(E_factors);
    N_MixParams = length(outer_MixParams);
    N_n = length(outer_n_values);
    total_combinations = N_MixParams * N_n * gridlength_E * gridlength_b;
    
    % Find MixParam(s)
    arr_size_MixParam = total_combinations / N_MixParams;
    index_MixParam = fix(slurm_array_task_id / arr_size_MixParam) + 1;
    MixParams = [outer_MixParams(index_MixParam)] % Just to screw up compatibility with older versions
    %%effective_MixParam = min(MixParam, 1-MixParam);
    rest_MixParam = mod(slurm_array_task_id, arr_size_MixParam);

    % Find n value
    arr_size_n = arr_size_MixParam / N_n;
    index_n = fix(rest_MixParam / arr_size_n) + 1;
    n = outer_n_values(index_n) 
    rest_n = mod(rest_MixParam, arr_size_n);
    
    % Find b and E_sensitive, E_resistant
    b_index = fix(rest_n/gridlength_E) + 1; % whole integer division +1 is the row
    dE_index = mod(rest_n, gridlength_E) + 1; % The rest +1 is the column   
    b = b_values(b_index)
    E_factor = E_factors(dE_index)
    E_resistant = E_sensitive * E_factor;
    sensitive_cell_line = [base_rate_sensitive, b, E_sensitive, n]
    resistant_cell_line = [base_rate_resistant, b, E_resistant, n]
    Params = [sensitive_cell_line resistant_cell_line]

    sensitive_alpha = sensitive_cell_line(1);
    resistant_alpha = resistant_cell_line(1);
    
    Nested_time = linspace(0,96,max_N_t); % Min T = 0, max T = 96
    Nested_conc = zeros(1,max_N_c);
    Nested_conc(2:max_N_c) = logspace(log10(0.05*E_sensitive),log10(20*E_resistant),(max_N_c-1)); % Min C = 0, then logspace from 10-1 to 10^1
    true_gr50_sensitive = find_gr50(sensitive_cell_line, Nested_conc);
    true_gr50_resistant = find_gr50(resistant_cell_line, Nested_conc);

    %NOISELESS_NESTED_DATA = generateDataEven(Nested_conc,Nested_time,max_R,MixParam,sensitive_cell_line,resistant_cell_line);
    NOISELESS_NESTED_DATA = generate_data_k_populations(Nested_conc,Nested_time,max_R,MixParams,Params);
    NESTED_DATA = zeros(max_R,max_N_c,max_N_t, N_repetitions);
    for i_rep=1:N_repetitions
        rng(seednr_1 + i_rep - 1); % set seed
        NESTED_DATA(:,:,:,i_rep) = max(0, NOISELESS_NESTED_DATA + Noise*randn(max_R,max_N_c,max_N_t));
    end

    
        %%lb=zeros(10,1);
        %%lb(4) = 1e-3; lb(8) = 1e-3; lb(10) = 1e-3; %To avoid Inf or NaN in log() terms we require positive E1,E2 and Noise
        %%ub=1000*ones(10,1);
        %%ub(1)=0.5; % Ratio parameter
        %%ub(2)=0.1;ub(3)=1;ub(4)=5;ub(5)=10; % Hill parameters for 1st pop
        %%ub(6)=0.1;ub(7)=1;ub(8)=5;ub(9)=10; % Hill parameters for 2nd pop
        %%ub(10)=5000; % Noise
    % xx = [p, alpha1, b1, E1, n1, alpha2, b2, E2, n2, sig]
    % Set upper and lower bounds for [alpha, b, E, n, sigma(both Hi, Lo)]
    % The E limits are then converted to log values for log spaced sampling
    %           a     b     E     n   sigma
    lb_all = [   0, 0.27, 1e-6, 0.01, 1e-6]; %[-0.1, 0, 1e-6,   0, 1e-6];
    ub_all = [ 0.1,    1,   max(Nested_conc),   10,  5e4];

    lb_all(3) = log10(lb_all(3));
    ub_all(3) = log10(ub_all(3));

    [lb, ub, A_inequality, b_inequality, nonlcon] = optimization_constraints_real_data(no_populations, INFER_SIGMA, USE_TWO_NOISE_LEVELS, lower_limit_mixparam, lower_limit_E_ratio, lb_all, ub_all);

    options = optimoptions(@fmincon,'FiniteDifferenceType','central','SpecifyObjectiveGradient',false,'Display','off','MaxFunctionEvaluations',5990,'MaxIterations',5990,'OptimalityTolerance',1.0e-10,'StepTolerance',1.0000e-10);
    mydate = sprintf('%s', datestr(now,'yyyy-mm-dd'));

    % Setting the base case
    N_c = base_case;
    N_t = base_case;
    R = base_case;
    
    p_estimates = zeros(N_repetitions, gridlength);
    gr50_sensitive_estimates = zeros(N_repetitions, gridlength);
    gr50_resistant_estimates = zeros(N_repetitions, gridlength);
    sensitive_alpha_estimates = zeros(N_repetitions, gridlength);
    resistant_alpha_estimates = zeros(N_repetitions, gridlength);
    estimates_for_histogram = zeros(gridlength, N_repetitions, num_optim, length(ub));
    tStart = tic;
    x_finals = zeros(N_repetitions, gridlength, length(lb'));
    f_vals = zeros(N_repetitions, gridlength);
    for i_rep=1:N_repetitions
        tStart_inner = tic;
        i_rep = i_rep
        for change_value_index=1:gridlength
            switch change_value
                case 'R'
                    R = R_values(change_value_index);
                case 'C'
                    N_c = N_c_values(change_value_index);
                case 'T'
                    N_t = N_t_values(change_value_index);
            end
            R_N_c_N_t = [R, N_c, N_t]
            
            stepsize_c = (max_N_c - 1) / (N_c - 1);
            stepsize_t = (max_N_t - 1) / (N_t - 1);
            Conc = Nested_conc(1:stepsize_c:max_N_c);
            Time = Nested_time(1:stepsize_t:max_N_t);
        
            DATA = NESTED_DATA(1:R, 1:stepsize_c:max_N_c, 1:stepsize_t:max_N_t, i_rep);
            
            %f=@(x)vectorized_objective_function(x,DATA,Conc,Time,R);
            if INFER_SIGMA
                if USE_TWO_NOISE_LEVELS
                    f=@(x)vectorized_objective_function_k_subpop_and_two_noise_levels(x,no_populations,DATA,Conc,Time,R,ConcT,TimeT);
                else
                    f=@(x)vectorized_objective_function_k_subpopulations(x,no_populations,DATA,Conc,Time,R);
                end
            else
                if USE_TWO_NOISE_LEVELS
                    error("Again: If using two noise level, sigma must also be inferred or you must implement an estimate for sigma")
                else 
                    % Estimate the standard deviation 
                    std_estimate = std(DATA); % estimates the standard deviation for all (time, concentration) combinations
                    Noise_estimate = mean(std_estimate, 'all'); % Finds the average noise estimate instead of using particulars
                    f=@(x)vectorized_objective_function_k_subpopulations([x' Noise_estimate]',no_populations,DATA,Conc,Time,R); % Drop Noise parameter
                    %f=@(x)vectorized_objective_function_k_subpopulations([x' Noise]',no_populations,DATA,Conc,Time,R); % Drop Noise parameter
                end
            end

            fval=inf;
            clear xx;
            clear ff;
            clear x_final;
            rng(42); %set seed to start iterations the same place for every version of R
            for nn=1:num_optim
                %%x0=rand(10,1).*(ub-lb) + lb;
                %%[xx,ff]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
                x0=rand(length(ub),1).*(ub-lb) + lb;
                x0(3:5:5*no_populations-2)=10.^(x0(3:5:5*no_populations-2));
                lb(3:5:5*no_populations-2)=1e-6; %To avoid Inf or NaN in log() terms we require positive E values 
                ub(3:5:5*no_populations-2)=max(Conc); % E
                if USE_E_THRESHOLD
                    [xx,ff]=fmincon(f,x0,A_inequality,b_inequality,[],[],lb,ub,nonlcon,options);
                else
                    [xx,ff]=fmincon(f,x0,A_inequality,b_inequality,[],[],lb,ub,[],options);
                end
        
                %%% Figure out which cell line estimate is sensitive: Minimize the error in the E parameter
                %%abs_sum_1_sensitive = abs(xx(4) - E_sensitive) + abs(xx(8) - E_resistant);
                %%abs_sum_2_sensitive = abs(xx(8) - E_sensitive) + abs(xx(4) - E_resistant);
                %%if abs_sum_1_sensitive < abs_sum_2_sensitive 
                %%    estimates_for_histogram(change_value_index, i_rep, nn, :) = xx;
                %%else
                %%    estimates_for_histogram(change_value_index, i_rep, nn, :) = [xx(1)' xx(6:9)' xx(2:5)' xx(10)'];
                %%end

                if ff<fval
                    x_final=xx;
                    fval=ff;
                end    
            end

            % Sort the parameters at once by E values.
            % xx = [p, alpha1, b1, E1, n1, alpha2, b2, E2, n2, sig]
            if x_final(4) > x_final(8)
                x_final(1) = 1-x_final(1);
                temp = x_final(2:5);
                x_final(2:5) = x_final(6:9);
                x_final(6:9) = temp;
            end
            x_finals(i_rep, change_value_index,:) = x_final';
            f_vals(i_rep, change_value_index,:) = fval;

            p_estimates(i_rep,change_value_index) = x_final(1);
            gr50_sensitive_estimates(i_rep,change_value_index) = find_gr50(x_final(2:5), Conc);
            gr50_resistant_estimates(i_rep,change_value_index) = find_gr50(x_final(6:9), Conc);
            sensitive_alpha_estimates(i_rep,change_value_index) = x_final(2);
            resistant_alpha_estimates(i_rep,change_value_index) = x_final(6);

            %%% Figure out which cell line estimate is sensitive: Minimize the error in the E parameter
            %%abs_sum_1_sensitive = abs(x_final(4) - E_sensitive) + abs(x_final(8) - E_resistant);
            %%abs_sum_2_sensitive = abs(x_final(8) - E_sensitive) + abs(x_final(4) - E_resistant);
            %%if abs_sum_1_sensitive < abs_sum_2_sensitive 
            %%    sensitive_alpha_estimates(i_rep,change_value_index) = x_final(2)
            %%    resistant_alpha_estimates(i_rep,change_value_index) = x_final(6)
            %%else
            %%    sensitive_alpha_estimates(i_rep,change_value_index) = x_final(6)
            %%    resistant_alpha_estimates(i_rep,change_value_index) = x_final(2)
            %%    % Correct this in x_finals too! 
            %%    temp = x_finals(i_rep, change_value_index, 2:5);
            %%    x_finals(i_rep, change_value_index, 2:5) = x_finals(i_rep, change_value_index, 6:9);
            %%    x_finals(i_rep, change_value_index, 6:9) = temp;
            %%end

        end % first for loop
        tEnd_inner = toc(tStart_inner)
    end % loop over N_repetitions
    tEnd = toc(tStart)

    % Compare with truth: 
    L1_p = abs(p_estimates - MixParams(1));
    logmaxratio_gr50_sensitive = log(max(gr50_sensitive_estimates./true_gr50_sensitive, true_gr50_sensitive./gr50_sensitive_estimates));
    logmaxratio_gr50_resistant = log(max(gr50_resistant_estimates./true_gr50_resistant, true_gr50_resistant./gr50_resistant_estimates));
    L1_sensitive_alpha = abs(sensitive_alpha_estimates - sensitive_alpha);
    L1_resistant_alpha = abs(resistant_alpha_estimates - resistant_alpha);
    
    if savebool
        if RUNNING_ON_CLUSTER
            save(strcat('/cluster/projects/nn9705k/evenmm/difference_data/inc-', num2str(change_value), '-gen-', num2str(generation_number), 'base_case', int2str(base_case), '-num_optim-', int2str(num_optim), '-N-rep-', int2str(N_repetitions), '-Noise-', int2str(Noise), '-seed-', int2str(seednr_1), '-MixParams-', num2str(MixParams(1)), '-n-', num2str(n), '-b-', num2str(b_index), '-dE-', num2str(dE_index), '.mat'), 'p_estimates', 'L1_p', 'sensitive_alpha_estimates', 'resistant_alpha_estimates', 'L1_sensitive_alpha', 'L1_resistant_alpha', 'MixParams', 'estimates_for_histogram', 'x_finals', 'f_vals', 'outer_MixParams', 'outer_n_values', 'gridlength_E', 'gridlength_b', 'b_values', 'base_rate_resistant', 'N_t', 'N_c', 'R', 'generation_number', 'resistant_cell_line', 'sensitive_cell_line', 'Conc', 'Time', 'logmaxratio_gr50_sensitive', 'logmaxratio_gr50_resistant', 'gr50_sensitive_estimates', 'gr50_resistant_estimates')
        else
            save(strcat('./data/difference/inc-', num2str(change_value), '-gen-', num2str(generation_number), 'base_case', int2str(base_case), '-num_optim-', int2str(num_optim), '-N-rep-', int2str(N_repetitions), '-Noise-', int2str(Noise), '-seed-', int2str(seednr_1), '-MixParams-', num2str(MixParams(1)), '-n-', num2str(n), '-b-', num2str(b_index), '-dE-', num2str(dE_index), '.mat'), 'p_estimates', 'L1_p', 'sensitive_alpha_estimates', 'resistant_alpha_estimates', 'L1_sensitive_alpha', 'L1_resistant_alpha', 'MixParams', 'estimates_for_histogram', 'x_finals', 'f_vals', 'outer_MixParams', 'outer_n_values', 'gridlength_E', 'gridlength_b', 'b_values', 'base_rate_resistant', 'N_t', 'N_c', 'R', 'generation_number', 'resistant_cell_line', 'sensitive_cell_line', 'Conc', 'Time', 'logmaxratio_gr50_sensitive', 'logmaxratio_gr50_resistant', 'gr50_sensitive_estimates', 'gr50_resistant_estimates')
        end
    end
    if ~RUNNING_ON_CLUSTER
        figure
        hold on
        title(strcat('Changing ',num2str(change_value), ', others fixed.'))
        xlabel('Noise')
        yyaxis left
        ylabel('Absolute error (MixParam)')
        ylim([0 0.3])
        for i_rep=1:N_repetitions
            plot(N_c_values, L1_p(i_rep, :), '-ko')
        end
        yyaxis right
        ylabel('Absolute error (rates)')
        ylim([0 0.025])
        for i_rep=1:N_repetitions
            plot(N_c_values, L1_sensitive_alpha(i_rep, :), '-bo')
            plot(N_c_values, L1_resistant_alpha(i_rep, :), '-ro')
        end
        %legend('Mixture param','alpha 1', 'alpha 2')
        set(gca, 'xtick', N_c_values)
        %set(gca, 'XDir','reverse')
        grid

        if savebool
            saveas(gcf,[pwd strcat('/plots/difference/gen-', num2str(generation_number), '-Nc-', int2str(N_c), '-Nt-', int2str(N_t), '-R-', int2str(R), '-num_optim-', int2str(num_optim), '-N-rep-', int2str(N_repetitions), '-seed-', int2str(seednr_1), '-MixParams-', num2str(MixParams(1)), '-n-', num2str(n), '-bindx-', num2str(b_index), '-dEindx-', num2str(dE_index), '-', mydate, '.png')])
        end
    end
    'finished'
end
