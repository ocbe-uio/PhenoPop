function [] = articleplot_2_limitations(N_true_populations, MixParams, Params, Noise, num_optim, caseno, RUNNING_ON_CLUSTER)
    set(0,'DefaultFigureVisible','off')
    % Generate data, infer and do model selection
    % 01.09.2021
    
    N_c = 17;
    Time = [0 12 24 36 48 60 72 84 96];
    N_t = length(Time);
    R = 4;
    
    %highest_rate = 0.1 % Not used when using optimization_constraints_real_data   % Highest exponential rate we assume there may be
    max_no_populations = 4 % The highest number of populations to try fitting to the data.
    seednr_1 = 45 % Seed for the optimization initialization
    lower_limit_mixparam = 0 % Each inferred mixture parameter must be greater than this.
    lower_limit_E_ratio = 1 % The smallest ratio (E_i / E_j) among inferred E parameters must be greater than this.
    USE_E_THRESHOLD = false
    INFER_SIGMA = true % If true, sigma is estimated by MLE. If false, the true value is given.
    USE_TWO_NOISE_LEVELS = true
    if USE_TWO_NOISE_LEVELS
        ConcT = 0.1; %thresholds for sigma
        TimeT = 48;
        if ~INFER_SIGMA
            error("If using two noise level, sigma must also be inferred or you must implement an estimate for sigma")
        end
    end
    PLOT_DATA = false
    PLOT_FINAL_FIT = true
    
    if RUNNING_ON_CLUSTER
        save_location = '/cluster/projects/nn9705k/evenmm/article/limitations/'
    else 
        save_location = './article/limitations/'
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These are copied, not coded: 
    % Old cell lines
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
    b3 = 0.5;
    b4 = 0.3;
    b5 = 0.4;
    b6 = 0.5;
    b7 = 0.3;
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
    
    tlabels = strings(1,N_t);
    for t_index=1:N_t
        tlabels(N_t-t_index+1) = strcat("T=", int2str(Time(t_index)));
    end
    stepsize_c = 2^4 / (N_c - 1);
    max_R = R;
    max_N_t = N_t;
    Nested_time = Time; 
    max_N_c = N_c;
    Nested_conc = zeros(1,max_N_c);
    Nested_conc(2:max_N_c) = logspace(log10(0.05*E_sensitive),log10(0.5),(max_N_c-1)); %log10(20*E_resistant),(max_N_c-1)); % Min C = 0, then logspace from 10-1 to 10^1
    Conc = Nested_conc(1:stepsize_c:max_N_c)
    
    % 2 Generate data
    NOISELESS_NESTED_DATA = generate_data_k_populations(Nested_conc,Nested_time,max_R,MixParams,Params);
    rng(seednr_1); % set seed for generating noise
    NESTED_DATA = max(0, NOISELESS_NESTED_DATA + Noise*randn(max_R,max_N_c,max_N_t));
    stepsize_c = 1; 
    stepsize_t = 1; 
    DATA = NESTED_DATA(1:R, 1:stepsize_c:max_N_c, 1:stepsize_t:max_N_t);
    
    colors_true = [
        0.847656250000000   0.371093750000000   0.007812500000000
        0.648437500000000   0.460937500000000   0.113281250000000
        0.2578 0.5273 0.9570
        0.457031250000000   0.437500000000000   0.699218750000000
        0.902343750000000   0.160156250000000   0.539062500000000
        0.398437500000000   0.648437500000000   0.117187500000000
        0.898437500000000   0.667968750000000   0.007812500000000
        0.398437500000000   0.398437500000000   0.398437500000000
    ];
    colors_estimated = [        
        0.9570 0.6640 0.2578
        0.105468750000000   0.617187500000000   0.464843750000000
        0.9570 0.2578 0.5039
        0.2578 0.9570 0.6172 
        0.7578 0.2578 0.9570 
        0                   0   0.726562000000000
        0.957000000000000   0.257800000000000   0.503900000000000
    ];
    switch N_true_populations
    case 1
        savestr = strcat('-Noise-', num2str(Noise), '-num_optim-', num2str(num_optim));
    case 2
        savestr = strcat('-Noise-', num2str(Noise), '-mix-', num2str(MixParams(1)), '-num_optim-', num2str(num_optim));
    case 3
        savestr = strcat('-Noise-', num2str(Noise), '-mix-', num2str(MixParams(1)), '-', num2str(MixParams(2)), '-num_optim-', num2str(num_optim));
    end
    
    plot_N_c = 1000;
    x = zeros(1,plot_N_c);
    x(2:plot_N_c) = logspace(log10(0.1*min(Conc(2:N_c))),log10(10*max(Conc(2:N_c))),(plot_N_c-1));
    min_x = 0.1*min(Conc(2:N_c));
    max_x = 10*max(Conc(2:N_c));
    
    %populations = generateDataEven(x,Time,R,MixParam,sensitive_cell_line,resistant_cell_line);
    populations = generate_data_k_populations(x,Time,R,MixParams,Params);
    max_cell_count = max(populations, [], "all");
    % N_c = 17
    % Conc = zeros(1,N_c)
    % Conc(2:N_c) = logspace(-1,1,(N_c-1)) % Min C = 0, then logspace from 10-1 to 10^1
    
    %fig = figure %('Position',[800 800 400 300]);
    %movegui(fig,[1275 630]); % x y positions of bottom left corner
    %h = axes;
    %set(h,'xscale','log')
    %%hold on
    %for t_index=N_t:-1:1
    %    semilogx(x, 100*populations(:,:,t_index)/max_cell_count, "linewidth",2)
    %    hold on
    %end
    %for c_index=2:N_c
    %    xline(Conc(c_index), "k")
    %end
    %ylabel('Tissue response (% max)')
    %xlabel('Drug dose')
    %legend(tlabels)
    %%title(['Cell counts' newline '             alpha          b             E          n' newline 'Sensitive   ', num2str(sensitive_cell_line), newline 'Resistant   ', num2str(resistant_cell_line)])
    %saveas(gcf,[pwd strcat(save_location, num2str(N_true_populations), '/case', num2str(caseno), '/0presentation-cell_counts_at_some_times.png')])
    
    x_finals_temp = zeros(1, max_no_populations, 5*max_no_populations);
    f_vals_temp = Inf(1, max_no_populations);
    negative_loglikelihood_values = Inf(1,max_no_populations);
    gr50values = zeros(max_no_populations,max_no_populations);
    true_gr50values = zeros(1,max_no_populations);
    
    if PLOT_DATA
        repeated_colors = [
            repmat([0.9570 0.2578 0.5039],R+1,1)
            repmat([0.2578 0.5273 0.9570],R+1,1)
            repmat([0.9570 0.6640 0.2578],R+1,1)
            repmat([0.2578 0.9570 0.6172 ],R+1,1)
            repmat([0.7578 0.2578 0.9570 ],R+1,1)
            repmat([0                   0   0.726562000000000],R+1,1)
            repmat([0.957000000000000   0.257800000000000   0.503900000000000],R+1,1)
            repmat([0.257800000000000   0.527300000000000   0.957000000000000],R+1,1)
            repmat([0.957000000000000   0.664000000000000   0.257800000000000],R+1,1)
        ];
    
        % Plot data points and true viability curves
        % a) train data
        plot_N_c = 1000;
        increment = 1; % Step in time. Increase to plot less of the data points
        plot_Time = Time(1:increment:N_t); %[0 24 48 72 96]; %[0 12 24 26 48 60 72 84 96];
        plot_N_t = length(plot_Time);
        tlabels = strings(1,N_t);
        for t_index=1:plot_N_t
            tlabels(plot_N_t-t_index+1) = strcat("T=", int2str(plot_Time(t_index)));
        end
        x = logspace(log10(0.1*min(Conc(2:N_c))),log10(10*max(Conc(2:N_c))),plot_N_c);
        plot_populations = generate_data_k_populations(x,plot_Time,1,MixParams,Params);
    
        fig = figure; %('Position',[800 800 400 300]);
        colororder(repeated_colors);
        movegui(fig,[1275 630]); % x y positions of bottom left corner
        h = axes;
        %set(h,'xscale','log')
        for t_index=plot_N_t:-1:1
            semilogx(x, plot_populations(:,:,t_index),"linewidth",2)
            hold on
            for r_index = 1:R
                semilogx(Conc(2:N_c), DATA(r_index,2:N_c,increment*t_index - (increment-1)), '.', 'MarkerSize', 10,'HandleVisibility','off')
            end
        end
        for c_index=2:N_c
            xline(Conc(c_index), "k",'HandleVisibility','off')
        end
        ylabel('Cell count')
        xlabel('Drug dose')
        legend(tlabels)
        title(['Observations, 1 dot for each replicate)']) % newline '             alpha          b             E          n' newline 'Sensitive   ', num2str(sensitive_cell_line), newline 'Resistant   ', num2str(resistant_cell_line)])
        saveas(gcf, [save_location, 'case-', num2str(caseno), '-Noise-', num2str(Noise), '.png'])
    
    %    repeated_colors = [
    %        repmat([0.9570 0.2578 0.5039],R,1)
    %        repmat([0.2578 0.5273 0.9570],R,1)
    %        repmat([0.9570 0.6640 0.2578],R,1)
    %        repmat([0.2578 0.9570 0.6172 ],R,1)
    %        repmat([0.7578 0.2578 0.9570 ],R,1)
    %        repmat([0                   0   0.726562000000000],R,1)
    %        repmat([0.957000000000000   0.257800000000000   0.503900000000000],R,1)
    %        repmat([0.257800000000000   0.527300000000000   0.957000000000000],R,1)
    %        repmat([0.957000000000000   0.664000000000000   0.257800000000000],R,1)
    %    ]; 
    %    
    %    % Plot observationss
    %    fig = figure; %('Position',[800 800 400 300]);
    %    movegui(fig,[1275 630]); % x y positions of bottom left corner
    %    h = axes;
    %    set(h,'xscale','log')
    %    colororder(repeated_colors);
    %    
    %    % Set time and conc
    %    increment = 1; % Step in time. Increase to plot less of the data points
    %    plot_Time = Time(1:increment:N_t); %[0 24 48 72 96]; %[0 12 24 26 48 60 72 84 96];
    %    plot_N_t = length(plot_Time);
    %    tlabels = strings(1,N_t);
    %    for t_index=1:plot_N_t
    %        tlabels(plot_N_t-t_index+1) = strcat("T=", int2str(plot_Time(t_index)));
    %    end
    %    plot_Conc = [min_x Conc(2:N_c)];
    %    for t_index=plot_N_t:-1:1
    %        hold on
    %        for r_index = 1:R
    %            if r_index == 1
    %                % Show handle
    %                semilogx(plot_Conc, DATA(r_index,1:N_c,increment*t_index - (increment-1)), '.', 'MarkerSize', 10,'HandleVisibility','on')
    %            else
    %                %semilogx(Conc(2:N_c), DATA(r_index,2:N_c,increment*t_index - (increment-1)), '.', 'MarkerSize', 10,'HandleVisibility','off')
    %                % Include Conc 0 data
    %                semilogx(plot_Conc, DATA(r_index,1:N_c,increment*t_index - (increment-1)), '.', 'MarkerSize', 10,'HandleVisibility','off')
    %                %semilogx(min_x, DATA(r_index,1,increment*t_index - (increment-1)), '.', 'MarkerSize', 10,'HandleVisibility','off')
    %            end
    %        end
    %    end
    %    for c_index=2:N_c
    %        xline(Conc(c_index), "k",'HandleVisibility','off')
    %    end
    %    ylabel('Cell count')
    %    xlabel('Drug dose')
    %    xlim([min_x, max_x])
    %    legend(tlabels)
    %    title(['Cell counts' '1 dot per replicate, dose 0 included'])
    %    saveas(gcf, [save_location, num2str(N_true_populations), '/case', num2str(caseno), '/logl-data.png'])
    end % if PLOT_DATA
    
    % Fit with 1,2,3,+++ subpopulations
    options = optimoptions(@fmincon,'FiniteDifferenceType','central','SpecifyObjectiveGradient',false,'Display','off','MaxFunctionEvaluations',5990,'MaxIterations',5990,'OptimalityTolerance',1.0e-10,'StepTolerance',1.0000e-10);
    mydate = sprintf('%s', datestr(now,'yyyy-mm-dd'));
    
    for jj = 1:N_true_populations
        true_gr50values(jj) = find_gr50(Params(4*jj-3 : 4*jj), Conc);
    end
    
    
    for ii = 1:max_no_populations
        tStart_inner = tic;
        % Infer ii populations
        no_populations = ii
    
        % Set upper and lower bounds for [alpha, b, E, n, sigma(both Hi, Lo)]
        % The E limits are then converted to log values for log spaced sampling
        %           a     b     E     n   sigma
        lb_all = [   0, 0.27, 1e-6, 0.01, 1e-6]; %[-0.1, 0, 1e-6,   0, 1e-6];
        ub_all = [ 0.1,    1,   max(Conc),   10,  5e4];
    
        lb_all(3) = log10(lb_all(3));
        ub_all(3) = log10(ub_all(3));
    
        [lb, ub, A_inequality, b_inequality, nonlcon] = optimization_constraints_real_data(no_populations, INFER_SIGMA, USE_TWO_NOISE_LEVELS, lower_limit_mixparam, lower_limit_E_ratio, lb_all, ub_all);
        
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
        
        rng(42); %set seed to start iterations the same place for every version of R
        for nn=1:(num_optim*no_populations)
            x0=rand(length(ub),1).*(ub-lb) + lb;
            x0(3:5:5*no_populations-2)=10.^(x0(3:5:5*no_populations-2));
            lb(3:5:5*no_populations-2)=1e-6; %To avoid Inf or NaN in log() terms we require positive E values 
            ub(3:5:5*no_populations-2)=max(Conc); % E
            if USE_E_THRESHOLD
                [xx,ff]=fmincon(f,x0,A_inequality,b_inequality,[],[],lb,ub,nonlcon,options);
            else
                [xx,ff]=fmincon(f,x0,A_inequality,b_inequality,[],[],lb,ub,[],options);
            end
            if ff<fval
                x_final=xx;
                fval=ff;
            end    
        end
        x_finals_temp(1, ii, 1:length(ub)) = x_final';
        f_vals_temp(1, ii) = fval;
        negative_loglikelihood_values(ii) = fval
        % BIC
        N_observations = N_c*N_t*R;
        Nparams = length(ub);
        % BIC:         %         k          *        ln(N)        + 2 * negative loglikelihood
        %bic_values_temp(ii) = Nparams * log(N_observations) + 2 * fval;
    
        inferred_mixtures = squeeze(x_finals_temp(1, no_populations, 5:5:5*(no_populations-1)))';
        inferred_mixtures = [inferred_mixtures, 1-sum(inferred_mixtures)];
        sorted_inferred_parameters = squeeze(x_finals_temp(1, no_populations,1:5*no_populations))';
        true_mixtures = [MixParams 1-sum(MixParams)];
        inferred_GR50s = zeros(1, no_populations);
    
        for jj = 1:no_populations
            gr50values(ii,jj) = find_gr50(sorted_inferred_parameters(5*jj-4:5*jj-1), Conc);
        end
    
        if PLOT_FINAL_FIT && ((no_populations <= N_true_populations) || (no_populations <= 2))
            plot_N_c = 200;
            x = zeros(1,plot_N_c);
            x(2:plot_N_c) = logspace(log10(min_x),log10(max_x),(plot_N_c-1));    
            % Plot all the inferred cell lines 
            newcolors0 = [colors_true(1:N_true_populations,:); colors_estimated(1:no_populations,:)];
            fig = figure;
            colororder(newcolors0);
            %movegui(fig,[1275 100]); % x y positions of bottom left corner
            h = axes;
            set(h,'xscale','log')
            hold on
            % Plot true and inferred
            for jj=1:N_true_populations
                true_parameters_jj = Params(4*jj-3 : 4*jj);
                y_true = ratefunc(true_parameters_jj', x);
                semilogx(x,y_true, '-', 'LineWidth', 3)
            end
            for jj=1:no_populations
                inferred_parameters = sorted_inferred_parameters(5*jj-4:5*jj-1);
                y_inferred = ratefunc(inferred_parameters', x);
                semilogx(x,y_inferred, '--', 'LineWidth', 3)
                inferred_GR50s(jj) = find_gr50(inferred_parameters, Conc);
            end
            ylabel('Growth rate')
            xlabel('Drug concentration')  
            title([strcat("Inferred growth rates")])
            true_legends = {};
            for ii = 1:N_true_populations
                str = strjoin({'True ' num2str(ii) '      Mixture ' num2str(true_mixtures(ii), '%.2f') ' GR50 ' num2str(true_gr50values(ii), '%.2e') });
                true_legends = [true_legends str];
            end
            inferred_legends = {};
            for ii = 1:no_populations
                str = strjoin({'Inferred ' num2str(ii) ' Mixture ' num2str(inferred_mixtures(ii), '%.2f') ' GR50 ' num2str(inferred_GR50s(ii), '%.2e') });
                inferred_legends = [inferred_legends str];
            end
            % These concentration placements are not informative; Look at data plotting
            %for c_index=2:N_c
            %    xline(Conc(c_index), "k",'HandleVisibility','off')
            %end
            legend([true_legends inferred_legends],'Location','southwest')
            saveas(gcf, [save_location, 'case-', num2str(caseno), '-logl-inferred_populations-no-pop-', num2str(no_populations), num2str(savestr), '.png'])
        end % if PLOT_FINAL_FIT
        tEnd_inner = toc(tStart_inner)
    end % of inference for ii populations loop
    
    true_mixtures
    inferred_mixtures
    true_gr50values
    inferred_GR50s
    
    negative_loglikelihood_values
    " Parameters 1 pop"
    squeeze(x_finals_temp(1,1,1:5))
    " Parameters 2 pops"
    squeeze(x_finals_temp(1,2,1:10))
    
    % Save data
    save([strcat(save_location, 'case-', num2str(caseno), '-negative_loglikelihood_values', savestr, '.mat')], 'negative_loglikelihood_values', 'x_finals_temp', 'true_mixtures', 'inferred_mixtures', 'true_gr50values', 'inferred_GR50s')
    
    % Plot the negative loglikelihood values
    fig = figure;
    hold on
    xlabel('Number of inferred populations')
    ylabel('Negative loglikelihood')
    title([strcat('Negative loglikelihood')])
    plot(1:max_no_populations, negative_loglikelihood_values, '.-k')
    saveas(gcf, [save_location, '-case-', num2str(caseno), '-negative_loglikelihood_values', num2str(savestr), '.png'])
    
    end % function
    