set(0,'DefaultFigureVisible','off')
tStartOuter = tic;

% Compute negative loglikelihood for real data, fitting 1,2,3,+++ populations
% 01.07.2021

S = load(strcat("./shannonTREATS.mat"));

case_array = ["TREAT E11", "TREAT E14", "TREAT E19", "TREAT E1975", "TREAT E41", "TREAT E827", "TREAT P11", "TREAT P14", "TREAT P19", "TREAT P1975", "TREAT P41", "TREAT P827"];
savestr_array = ["TREAT_E_11", "TREAT_E_14", "TREAT_E_19", "TREAT_E_1975", "TREAT_E_41", "TREAT_E_827", "TREAT_P_11", "TREAT_P_14", "TREAT_P_19", "TREAT_P_1975", "TREAT_P_41", "TREAT P827"];

for casenumber = 1:length(case_array) %%%%%%%%%%%%%%%% parfor loop over drugs
casename = case_array(casenumber);
savestr1 = savestr_array(casenumber);
tStart = tic;
%parfor casenumber = 1:4 %%%%%%%%%%%%%%%% parfor loop over drugs

num_optim = 3000 % Number of starting points in maximum likelihood optimization.
%highest_rate = 0.1 % Not used when using optimization_constraints_real_data   % Highest exponential rate we assume there may be
max_no_populations = 5 % The highest number of populations to try fitting to the data.
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
PLOT_DATA = true
PLOT_FINAL_FIT = true

%x_finals = zeros(N_settings, 1, max_no_populations, 5*max_no_populations); % Estimated parameters on training sets (k-1 folds)
%f_vals = Inf(N_settings, 1, max_no_populations); % negative loglikelihood values on training sets (k-1 folds)
%best_no_pops_from_average_loglikelihood = zeros(1, N_settings); % The number of populations with highest average loglikelihood score across k folds, for N_settings
%best_no_pops_from_bic = zeros(1, N_settings);
%average_loglikelihoods = zeros(N_settings, max_no_populations);
%bic_values = zeros(N_settings, max_no_populations);

colors_estimated = [        
    0.9570 0.6640 0.2578
    0.105468750000000   0.617187500000000   0.464843750000000
    0.9570 0.2578 0.5039
    0.2578 0.9570 0.6172 
    0.7578 0.2578 0.9570 
    0                   0   0.726562000000000
    0.957000000000000   0.257800000000000   0.503900000000000
];

switch casenumber
case 1
    DATA = S.TREAT_E_11;
case 2
    DATA = S.TREAT_E_14;
case 3
    DATA = S.TREAT_E_19;
case 4
    DATA = S.TREAT_E_1975;
case 5
    DATA = S.TREAT_E_41;
case 6
    DATA = S.TREAT_E_827;
case 7
    DATA = S.TREAT_P_11;
case 8
    DATA = S.TREAT_P_14;
case 9
    DATA = S.TREAT_P_19;
case 10
    DATA = S.TREAT_P_1975;
case 11
    DATA = S.TREAT_P_41;
case 12
    DATA = S.TREAT_P_827;
end
savestr = strcat(savestr1, '-num_optim-', num2str(num_optim));

Time = [0 24 48 72];
Conc = [0 0.1 1 10];

R = size(DATA, 1)
N_c = size(DATA, 2)
N_t = size(DATA, 3)

x_finals_temp = zeros(1, max_no_populations, 5*max_no_populations);
f_vals_temp = Inf(1, max_no_populations);
negative_loglikelihood_values = Inf(1,max_no_populations);

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
    set(h,'xscale','log')
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
    title([strcat(newline, casename) 'Cell counts' '1 dot per replicate, dose 0 included'])
    saveas(gcf, [pwd, '/plots/shannon_data/logl-data-', num2str(savestr1), '.png'])
end % if PLOT_DATA

% Fit with 1,2,3,+++ subpopulations
options = optimoptions(@fmincon,'FiniteDifferenceType','central','SpecifyObjectiveGradient',false,'Display','off','MaxFunctionEvaluations',5990,'MaxIterations',5990,'OptimalityTolerance',1.0e-10,'StepTolerance',1.0000e-10);
mydate = sprintf('%s', datestr(now,'yyyy-mm-dd'));

for ii = 1:max_no_populations
    tStart_inner = tic;
    % Infer ii populations
    no_populations = ii
    % Set upper and lower bounds for [alpha, b, E, n, sigma(both Hi, Lo)]
    % The E limits are then converted to log values for log spaced sampling

    remaining_cell_fraction = 0.001;
    ub_alpha = 0.06;
    delta_t = Time(2) - Time(1);
    lb_b = find_lower_limit_b(remaining_cell_fraction, ub_alpha, delta_t);
    %            alpha,    b,                  E,   n, sigma(both Hi, Lo)
    lb_all = [       0, lb_b,               1e-6, 0.1, 1e-6]; %[-0.1, 0, 1e-6,   0, 1e-6];
    ub_all = [ub_alpha,    1, Conc(length(Conc)),  10,  5e4];

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
    negative_loglikelihood_values(ii) = fval
    % BIC
    N_observations = N_c*N_t*R;
    Nparams = length(ub);
    % BIC:         %         k          *        ln(N)        + 2 * negative loglikelihood
    %bic_values_temp(ii) = Nparams * log(N_observations) + 2 * fval;

    inferred_mixtures = squeeze(x_finals_temp(1, no_populations, 5:5:5*(no_populations-1)))';
    inferred_mixtures = [inferred_mixtures, 1-sum(inferred_mixtures)];
    sorted_inferred_parameters = squeeze(x_finals_temp(1, no_populations,1:5*no_populations))';

    if PLOT_FINAL_FIT
        plot_N_c = 1000;
        x = zeros(1,plot_N_c);
        x(2:plot_N_c) = logspace(log10(min_x),log10(max_x),(plot_N_c-1));    
        % Plot all the inferred cell lines 
        newcolors0 = colors_estimated(1:no_populations,:);
        fig = figure;
        colororder(newcolors0);
        %movegui(fig,[1275 100]); % x y positions of bottom left corner
        h = axes;
        set(h,'xscale','log')
        hold on
        for jj=1:no_populations
            inferred_parameters = sorted_inferred_parameters(5*jj-4:5*jj-1);
            y_inferred = ratefunc(inferred_parameters', x);
            semilogx(x,y_inferred, '--', 'LineWidth', 3)
        end
        ylabel('Growth rate')
        xlabel('Drug concentration')  
        title([strcat(newline, casename) strcat("Inferred growth rates")])
        inferred_legends = {};
        for ii = 1:no_populations
            str = strjoin({'Inferred ' num2str(ii) ' Mixture ' num2str(inferred_mixtures(ii), '%.2f')});
            inferred_legends = [inferred_legends str];
        end
        % These concentration placements are not informative; Look at data plotting
        %for c_index=2:N_c
        %    xline(Conc(c_index), "k",'HandleVisibility','off')
        %end
        legend([inferred_legends],'Location','southwest')
        saveas(gcf, [pwd, '/plots/shannon_data/logl-inferred_populations-', num2str(savestr1), '-no-pop-', num2str(no_populations), '-num_optim-', num2str(num_optim), '.png'])
    end % if PLOT_FINAL_FIT
    tEnd_inner = toc(tStart_inner)
end % of inference for ii populations loop

negative_loglikelihood_values
" Parameters 1 pop"
squeeze(x_finals_temp(1,1,1:6))
" Parameters 2 pops"
squeeze(x_finals_temp(1,2,1:11))

% Save data
save([strcat('./data/shannon_data/negative_loglikelihood_values-', savestr, '.mat')], 'negative_loglikelihood_values', 'x_finals_temp')

% Plot the negative loglikelihood values
fig = figure;
hold on
xlabel('Number of inferred populations')
ylabel('Negative loglikelihood')
title([strcat(newline, casename) strcat('Negative loglikelihood')])
plot(1:max_no_populations, negative_loglikelihood_values, '.-k')
saveas(gcf, [pwd, '/plots/shannon_data/negative_loglikelihood_values-', num2str(savestr), '.png'])

tEnd = toc(tStart)
end %%%% case number loop
tEndOuter = toc(tStartOuter)
