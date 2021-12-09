% Turn off figures popping up 
set(0,'DefaultFigureVisible','off')
tStartOuter = tic;

% Compute negative loglikelihood for real data, fitting 1,2,3,+++ populations
% 11.06.2021

%patient_number = "MM210819"
%patient_number = "MM130120" % aka MM720
%patient_number = "MM19520"
%patient_number = "MM3620"
%patient_number = "MM1420"

patno_array = ["MM210819", "MM130120", "MM19520", "MM3620", "MM1420"];
for patient_index = 1:length(patno_array)
patient_number = patno_array(patient_index)

for drugnumber = 1:4 %%%%%%%%%%%%%%%% parfor loop over drugs
drugnumber
tStart = tic;
%parfor drugnumber = 1:4 %%%%%%%%%%%%%%%% parfor loop over drugs

num_optim = 3000 % Number of starting points in maximum likelihood optimization.
% 300 worked with current find_gr50 version (30.08.2021)
%highest_rate = 0.1 % Not used when using optimization_constraints_real_data   % Highest exponential rate we assume there may be
max_no_populations = 5 % The highest number of populations to try fitting to the data.
seednr_1 = 45 % Seed for the optimization initialization
lower_limit_mixparam = 0 % Each inferred mixture parameter must be greater than this.
lower_limit_E_ratio = 1 % The smallest ratio (E_i / E_j) among inferred E parameters must be greater than this.
USE_E_THRESHOLD = false
INFER_SIGMA = true % If true, sigma is estimated by MLE. If false, the true value is given.
USE_TWO_NOISE_LEVELS = true
if USE_TWO_NOISE_LEVELS
    ConcT = 10;
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

S = load(strcat("./data/real-MM-cellines/", patient_number, ".mat"));
Conc1 = S.Conc1;
Conc2 = S.Conc2;
Conc3 = S.Conc3;
Conc4 = S.Conc4;
Drug1 = S.Drug1;
Drug2 = S.Drug2;
Drug3 = S.Drug3;
Drug4 = S.Drug4;
Time = S.Time;

colors_estimated = [        
    0.9570 0.6640 0.2578
    0.105468750000000   0.617187500000000   0.464843750000000
    0.9570 0.2578 0.5039
    0.2578 0.9570 0.6172 
    0.7578 0.2578 0.9570 
    0                   0   0.726562000000000
    0.957000000000000   0.257800000000000   0.503900000000000
];

switch drugnumber
case 1
    DATA = Drug1;
    Conc = Conc1;
    % Time = Time;
case 2  
    DATA = Drug2;
    Conc = Conc2;
    % Time = Time;
case 3
    DATA = Drug3;
    Conc = Conc3;
    % Time = Time;
case 4
    DATA = Drug4;
    Conc = Conc4;
    % Time = Time;
end

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
    title([strcat(newline, patient_number, ' drug ', num2str(drugnumber)) 'Cell counts' '1 dot per replicate, dose 0 included'])
    saveas(gcf, [pwd, '/plots/real-MM-cellines/', num2str(patient_number), '/logl-data-', num2str(patient_number), '-drug-', num2str(drugnumber), '.png'])
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
    switch patient_number
    case "MM210819"
        lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        ub_all = [ 0.1, 1,  2e3, 100, 5e4];
    case "MM130120" % aka MM720
        lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        ub_all = [ 0.1, 1,  5e4, 100, 1e6];
        %if no_populations == 1
        %    lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        %    ub_all = [ 0.1, 1,  5e4, 100, 1e6];
        %else
        %    lb_all = [-0.1, 0,  1e-6,   0, 10^-6];
        %    ub_all = [ 0.1, 1, 5e3, 100, 1.2*10^6];
        %end
    case "MM19520"
        lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        ub_all = [ 0.1, 1,  5e3, 100, 1.5e5];
        %if no_populations == 1
        %    lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        %    ub_all = [ 0.1, 1,  5e3, 100, 1.5e5];
        %else
        %    lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        %    ub_all = [ 0.1, 1,  2e3, 100, 1.5e5];
        %end
    case "MM3620"
        lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        ub_all = [ 0.1, 1,  5e3, 100, 2.5e5];
        %if no_populations == 1
        %    lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        %    ub_all = [ 0.1, 1,  5e3, 100, 2.5e5];
        %else
        %    lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        %    ub_all = [ 0.1, 1,  2e3, 100, 2.5e5];
        %end
    case "MM1420"
        lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        ub_all = [ 0.1, 1,  1e5, 100, 1.5e5];
        %if no_populations == 1
        %    lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        %    ub_all = [ 0.1, 1,  1e5, 100, 1.5e5];
        %else
        %    lb_all = [-0.1, 0, 1e-6,   0, 1e-6];
        %    ub_all = [ 0.1, 1,  1e4, 100, 1.5e5];
        %end
    end
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
    inferred_GR50s = zeros(1, no_populations);

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
            inferred_GR50s(jj) = find_gr50(inferred_parameters, Conc);
        end
        ylabel('Growth rate')
        xlabel('Drug concentration')  
        title([strcat(newline, patient_number, ' drug ', num2str(drugnumber)) strcat("Inferred growth rates")])
        inferred_legends = {};
        for jj = 1:no_populations
            str = strjoin({'Clone ' num2str(jj) ' Mixture ' num2str(inferred_mixtures(jj), '%.2f') ' GR50 ' num2str(inferred_GR50s(jj), '%.2e') });
            inferred_legends = [inferred_legends str];
        end
        % These concentration placements are not informative; Look at data plotting
        %for c_index=2:N_c
        %    xline(Conc(c_index), "k",'HandleVisibility','off')
        %end
        legend([inferred_legends],'Location','southwest')
        saveas(gcf, [pwd, '/plots/real-MM-cellines/', num2str(patient_number), '/logl-inferred_populations-', num2str(patient_number), '-drug-', num2str(drugnumber), '-no-pop-', num2str(no_populations), '-num_optim-', num2str(num_optim), '.png'])
    end % if PLOT_FINAL_FIT
    tEnd_inner = toc(tStart_inner);
end % of inference for ii populations loop

negative_loglikelihood_values
inferred_GR50s
" Parameters 1 pop"
squeeze(x_finals_temp(1,1,1:6))
" Parameters 2 pops"
squeeze(x_finals_temp(1,2,1:11))

savestr = strcat(num2str(patient_number), '-drug-', num2str(drugnumber), '-num_optim-', num2str(num_optim));

% Save data
save([strcat('./plots/real-MM-cellines/negLL/negative_loglikelihood_and_GR50-', savestr, '.mat')], 'negative_loglikelihood_values', 'x_finals_temp', 'inferred_GR50s')

% Plot the negative loglikelihood values
fig = figure;
hold on
xlabel('Number of inferred populations')
ylabel('Negative loglikelihood')
title(strcat(patient_number, ' drug ', num2str(drugnumber)))
plot(1:max_no_populations, negative_loglikelihood_values, '.-k')
saveas(gcf, [pwd, '/plots/real-MM-cellines/', num2str(patient_number), '/negative_loglikelihood_values-', savestr, '.png'])

tEnd = toc(tStart)
end %%%% drug number loop
end %%%% patient number loop
tEndOuter = toc(tStartOuter)
