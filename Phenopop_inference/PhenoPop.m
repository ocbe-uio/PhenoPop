% PhenoPop
% This function runs PhenoPop for up to "max_no_populations_optional" subpopulations and plots the estimated dose-response curves for each number of subpopulations. 
% The negative loglikelihoods, AIC and BIC values are plotted and saved to the current repository. 

% Required data format: 
% The user specifies the number of replicates NR and provides a list of observed concentrations (Concentration_array) and a list of observed timepoints (Time_array).
% The "datafile" argument must point to a csv file of dimensions (NR*NT, NC) containing the experiment data:
%   The first column must contain all data for the first concentration in Concentration_array; 
%   first all replicates for the first timepoint followed by all replicates for the second timepoint, and so on.

% ----- Explanation o function arguments: -----
% experiment_name: a string for plot titles 
% NR: Number of replicates
% Concentration_array:  Array containing the concentrations in the drug screen data. Concentrations can be specified in the user's preferred unit. 
%                       Their GR50 result will be in that unit. The threshold ConcT should be changed to a cencentration dose somewhere between the min and the max of the concentrations. 
% Time_array: Array containing the time points (in hours) in the drug screen data
% datafile: a string specifying location of the data file with dimensions (NR*NT, NC)
% ----- Optional arguments: -----
% max_no_populations_optional: The maximum number of populations to fit the model with
% num_optim_optional: The number of individual local optimizations performed. Increase if lack of convergence is suspected
% lower_bounds_optional: Lower bounds on parameters [alpha, b, E, n, sigma(both Hi, Lo)]. The default are [-0.1, 0, 1e-6,   0, 1e-6]
% upper_bounds_optional: Higher bounds on parameters [alpha, b, E, n, sigma(both Hi, Lo)]. The default are [ 0.1, 1,  2e3, 100, 5e4]
% USE_TWO_NOISE_LEVELS_optional: If false, the noise is assumed to be the same for all data. 
%                                If true, two noise levels are used: the higher "sigH" is used when time >= TimeT && Conc <= ConcT, otherwise "sigL" is used. 
% ConcT_optional: Threshold value
% TimeT_optional: Threshold value

function [x_finals_temp, f_vals_temp, negative_loglikelihood_values, AIC_values, BIC_values] = PhenoPop(experiment_name, NR, Concentration_array, Time_array, datafile, max_no_populations_optional, num_optim_optional, lower_bounds_optional, upper_bounds_optional, USE_TWO_NOISE_LEVELS_optional, ConcT_optional, TimeT_optional) 
    disp("Starting PhenoPop")
    disp("Loading data...")
    tStart = tic;
    NC = length(Concentration_array);
    NT = length(Time_array);
    Conc = Concentration_array;
    Time = Time_array;
    if nargin > 5
        max_no_populations = max_no_populations_optional;
    else 
        max_no_populations = 4;
    end
    if nargin > 6
        num_optim = num_optim_optional;
    else 
        num_optim = 1000;
    end
    if nargin > 7
        lb_all = lower_bounds_optional;
        ub_all = upper_bounds_optional;
    else 
        lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        ub_all = [ 0.1, 1,  2e3, 100, 5e4];
    end
    if nargin > 9
        PLOT_DATA = PLOT_DATA_optional;
    else 
        PLOT_DATA = true;
    end
    if nargin > 10
        USE_TWO_NOISE_LEVELS = USE_TWO_NOISE_LEVELS_optional;
        ConcT_optional = ConcT;
        TimeT_optional = TimeT;
    else
        USE_TWO_NOISE_LEVELS = true;
        ConcT = 10;
        TimeT = 48;
    end
    lb_all(3) = log10(lb_all(3));
    ub_all(3) = log10(ub_all(3));
    % Read data
    flat_data = table2array(readtable(datafile));

    % Reshape flat data to 3D
    DATA = zeros(NR,NC,NT);
    for time_index = 0:NT-1 
        for conc_index = 1:NC
            DATA(1:NR,conc_index,time_index+1) = squeeze(flat_data(NR*time_index+1:NR*(time_index+1),conc_index));
        end 
    end 
    R = size(DATA, 1)
    N_c = size(DATA, 2)
    N_t = size(DATA, 3)

    seednr_1 = 45; % Seed for the optimization initializations
    lower_limit_mixparam = 0; % Each inferred mixture parameter must be greater than this.
    lower_limit_E_ratio = 1; % The smallest ratio (E_i / E_j) among inferred E parameters must be greater than this.
    USE_E_THRESHOLD = false;
    INFER_SIGMA = true; % If true, sigma is estimated by MLE. If false, the true value is given.
    colors_estimated = [
        0.9570 0.6640 0.2578
        0.105468750000000   0.617187500000000   0.464843750000000
        0.9570 0.2578 0.5039
        0.2578 0.9570 0.6172 
        0.7578 0.2578 0.9570 
        0                   0   0.726562000000000
        0.957000000000000   0.257800000000000   0.503900000000000
    ];
        
    x_finals_temp = zeros(1, max_no_populations, 5*max_no_populations);
    f_vals_temp = Inf(1, max_no_populations);
    negative_loglikelihood_values = Inf(1,max_no_populations);
    all_gr50_values = zeros(max_no_populations, max_no_populations);

    conclabels = strings(1,N_c);
    for c_index=1:N_c
        conclabels(c_index) = num2str(Conc(c_index));
    end
    newcolors = [0.267004 0.004874 0.329415
    0.282623 0.140926 0.457517
    0.253935 0.265254 0.529983
    0.206756 0.371758 0.553117
    0.163625 0.471133 0.558148
    0.127568 0.566949 0.550556
    0.134692 0.658636 0.517649
    0.266941 0.748751 0.440573
    0.477504 0.821444 0.318195
    0.741388 0.873449 0.149561
    0.993248 0.906157 0.143936];
    newcolors = [
        repmat([0.267004 0.004874 0.329415],R,1)
        repmat([0.253935 0.265254 0.529983],R,1)
        repmat([0.206756 0.371758 0.553117],R,1)
        repmat([0.163625 0.471133 0.558148],R,1)
        repmat([0.127568 0.566949 0.550556],R,1)
        repmat([0.134692 0.658636 0.517649],R,1)
        repmat([0.266941 0.748751 0.440573],R,1)
        repmat([0.477504 0.821444 0.318195],R,1)
        repmat([0.993248 0.906157 0.143936],R,1)
        ]; 

    if PLOT_DATA
        repeated_colors = [
            repmat([0.9570 0.2578 0.5039],R,1)
            repmat([0.2578 0.5273 0.9570],R,1)
            repmat([0.9570 0.6640 0.2578],R,1)
            repmat([0.2578 0.9570 0.6172 ],R,1)
            repmat([0.7578 0.2578 0.9570 ],R,1)
            repmat([0                   0   0.726562000000000],R,1)
            repmat([0.957000000000000   0.257800000000000   0.503900000000000],R,1)
            repmat([0.257800000000000   0.527300000000000   0.957000000000000],R,1)
            repmat([0.957000000000000   0.664000000000000   0.257800000000000],R,1)
        ]; 
        
        % Plot observationss
        fig = figure; %('Position',[800 800 400 300]);
        movegui(fig,[1275 630]); % x y positions of bottom left corner
        h = axes;
        set(h,'xscale','log');
        colororder(repeated_colors);
        min_x = 0.1*min(Conc(2:N_c));
        max_x = 10*max(Conc(2:N_c));
        
        % Set time and conc
        increment = 1; % Step in time. Increase to plot less of the data points
        plot_Time = Time(1:increment:N_t); %[0 24 48 72 96]; %[0 12 24 26 48 60 72 84 96];
        plot_N_t = length(plot_Time);
        tlabels = strings(1,N_t);
        for t_index=1:plot_N_t
            tlabels(plot_N_t-t_index+1) = strcat("T=", int2str(plot_Time(t_index)));
        end
        plot_Conc = [min_x Conc(2:N_c)];
        for t_index=plot_N_t:-1:1
            hold on
            for r_index = 1:R
                if r_index == 1
                    % Show handle
                    semilogx(plot_Conc, DATA(r_index,1:N_c,increment*t_index - (increment-1)), '.', 'MarkerSize', 10,'HandleVisibility','on')
                else
                    %semilogx(Conc(2:N_c), DATA(r_index,2:N_c,increment*t_index - (increment-1)), '.', 'MarkerSize', 10,'HandleVisibility','off')
                    % Include Conc 0 data
                    semilogx(plot_Conc, DATA(r_index,1:N_c,increment*t_index - (increment-1)), '.', 'MarkerSize', 10,'HandleVisibility','off')
                    %semilogx(min_x, DATA(r_index,1,increment*t_index - (increment-1)), '.', 'MarkerSize', 10,'HandleVisibility','off')
                end
            end
        end
        for c_index=2:N_c
            xline(Conc(c_index), "k",'HandleVisibility','off')
        end
        ylabel('Cell count')
        xlabel('Drug dose')
        xlim([min_x, max_x])
        legend(tlabels)
        title([strcat(newline, experiment_name) 'Cell counts' '1 dot per replicate, dose 0 included'])
        saveas(gcf, [pwd, '/logl-data-', num2str(experiment_name), '.png'])
    
        % Plot in time 
        figure
        colororder(newcolors);
        DATA(r_index,c_index,1:N_t);
        hold on
        for c_index=1:N_c
            for r_index=1:R
                plot(Time, squeeze(DATA(r_index,c_index,1:N_t)), '.-', 'MarkerSize', 10,'HandleVisibility','off')
            end
        end
        ylabel('Cell count')
        xlabel('Time (h)')
        ylim([0 inf])
        %legend(conclabels(1:N_c), 'location', 'northeast') % did not work, empty
        title([strcat(newline, experiment_name) 'All replicates for each concentration'])
        saveas(gcf, [pwd, '/logl-data-every-replicate-', num2str(experiment_name), '.png'])
    end % if PLOT_DATA

    % Fit with 1,2,3,+++ subpopulations
    options = optimoptions(@fmincon,'FiniteDifferenceType','central','SpecifyObjectiveGradient',false,'Display','off','MaxFunctionEvaluations',5990,'MaxIterations',5990,'OptimalityTolerance',1.0e-10,'StepTolerance',1.0000e-10);
    mydate = sprintf('%s', datestr(now,'yyyy-mm-dd'));

    for ii = 1:max_no_populations
        tStart_inner = tic;
        % Infer ii populations
        no_populations = ii;
        disp(['Fitting model with ' num2str(no_populations) ' populations...'])

        % Set upper and lower bounds for [alpha, b, E, n, sigma(both Hi, Lo)]
        % The E limits are then converted to log values for log spaced sampling
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
        for nn=1:num_optim
            x0=rand(length(ub),1).*(ub-lb) + lb;
            x0(3:5:5*no_populations-2)=10.^(x0(3:5:5*no_populations-2));
            lb(3:5:5*no_populations-2)=1e-6; %To avoid Inf or NaN in log() terms we require positive E values 
            ub(3:5:5*no_populations-2)=1e4; % E
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
        negative_loglikelihood_values(ii) = fval;

        inferred_mixtures = squeeze(x_finals_temp(1, no_populations, 5:5:5*(no_populations-1)))';
        inferred_mixtures = [inferred_mixtures, 1-sum(inferred_mixtures)];
        sorted_inferred_parameters = squeeze(x_finals_temp(1, no_populations,1:5*no_populations))';
        inferred_GR50s = zeros(1, no_populations);
        all_gr50_values(ii,1:no_populations) = inferred_GR50s;
        tEnd_inner = toc(tStart_inner);
    end % of inference for ii populations loop

    for no_populations = 1:max_no_populations
        disp(['------------ ' num2str(no_populations) ' populations model' '------------'])
        disp(['Negative loglikelihood ' num2str(negative_loglikelihood_values(no_populations)) ':'])
        for m = 1:no_populations
            disp(['-- Parameters and GR50 for population ' num2str(m) ':'])
            parameters_m = squeeze(x_finals_temp(1,no_populations,5*m-4:5*m-1));
            if no_populations == 1
                disp(['alpha: ' num2str(parameters_m(1))])
                disp(['b: ' num2str(parameters_m(2))])
                disp(['E: ' num2str(parameters_m(3))])
                disp(['n: ' num2str(parameters_m(4))])
            else
                if m < no_populations
                    mixparam_m = squeeze(x_finals_temp(1,no_populations,5*m));
                else
                    mixparam_m = 1-sum(x_finals_temp(1,no_populations,5:5:5*(no_populations-1)));
                end
                disp(['alpha: ' num2str(parameters_m(1))])
                disp(['b: ' num2str(parameters_m(2))])
                disp(['E: ' num2str(parameters_m(3))])
                disp(['n: ' num2str(parameters_m(4))])
                disp(['pi: ' num2str(mixparam_m)])
            end
            disp(['GR50 value ' num2str(all_gr50_values(no_populations,m)) ':'])
        end
        disp(['-- Estimated noise levels:'])
        if USE_TWO_NOISE_LEVELS
            disp(['sigH: ' num2str(squeeze(x_finals_temp(1,no_populations,5*no_populations)))])
            disp(['sigL: ' num2str(squeeze(x_finals_temp(1,no_populations,5*no_populations+1)))])
        else
            disp(['sigma: ' num2str(squeeze(x_finals_temp(1,no_populations,5*no_populations)))])
        end
    end

    % Save data
    savestr = strcat(num2str(experiment_name), '-num_optim-', num2str(num_optim));
    save([strcat(pwd, '/results-', savestr, '.mat')], 'negative_loglikelihood_values', 'x_finals_temp', 'inferred_GR50s')

    % BIC and AIC
    N_observations = N_c*N_t*R;
    Nparams = length(ub);
    BIC_values = Nparams * log(N_observations) + 2 * negative_loglikelihood_values;
    AIC_values = 2*Nparams + 2 * negative_loglikelihood_values;

    disp(['All negative loglikelihood values for models with 1,2,etc populations:'])
    disp([num2str(negative_loglikelihood_values)])
    disp(['All AIC values for models with 1,2,etc populations:'])
    disp([num2str(AIC_values)])
    disp(['All BIC values for models with 1,2,etc populations:'])
    disp([num2str(BIC_values)])

    % Plot the negative loglikelihood values
    fig = figure;
    hold on
    xlabel('Number of inferred populations')
    ylabel('Negative loglikelihood')
    title(['Negative loglikelihood ', experiment_name])
    plot(1:max_no_populations, negative_loglikelihood_values, '.-k')
    [min_nll, min_pos] = min(negative_loglikelihood_values);
    plot(min_pos, negative_loglikelihood_values(min_pos), '.', 'MarkerSize', 30)
    saveas(gcf, [pwd, '/negative_loglikelihood_values-', savestr, '.png'])

    % Plot the AIC values
    fig = figure;
    hold on
    xlabel('Number of inferred populations')
    ylabel('AIC')
    title(['AIC ', experiment_name])
    plot(1:max_no_populations, AIC_values, '.-k')
    [min_nll, min_pos] = min(AIC_values);
    plot(min_pos, AIC_values(min_pos), '.', 'MarkerSize', 30)
    saveas(gcf, [pwd, '/AIC_values-', savestr, '.png'])

    % Plot the BIC values
    fig = figure;
    hold on
    xlabel('Number of inferred populations')
    ylabel('BIC')
    title(['BIC ', experiment_name])
    plot(1:max_no_populations, BIC_values, '.-k')
    [min_nll, min_pos] = min(BIC_values);
    plot(min_pos, BIC_values(min_pos), '.', 'MarkerSize', 30)
    saveas(gcf, [pwd, '/BIC_values-', savestr, '.png'])
    tEnd = toc(tStart);
    disp("Total time spent (seconds): ")
    disp(tEnd)
    disp(['Results (negative loglikelihood values, GR50 values and parameter estimates) are saved in:'])
    disp(['"results-' savestr '.mat"'])
end