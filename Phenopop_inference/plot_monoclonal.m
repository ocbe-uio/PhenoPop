% Unpack and plot Dagim's monoclonal data
% The data comes in 4 different settings:   - Sensitive 500 and 1000 initial cells
%                                           - Resistant 250 and 500 initial cells
% All 4 carried out with Bright field imaging (BF)
% Both of the sensitive are also done with Green fluorescent protein imaging (GFP)
% ... for 7 replicates, 25 time points (0 to 72 hours with 3 hour intervals), and 11 concentrations.
% 
% In this script:
% A) Data cleaning by outlier removal
% B) Cell line inference on the cleaned data
INFERENCE = false
% C) Plot fitted vs estimated rates
PLOT_FITTED_VS_ESTIMATED_RATES = false
% D) Plotting the clean data
PLOT_CLEAN_DATA = true

% Each replicate has a different initial cell count. Therefore the first datapoint of each replicate is used to predict. 
FIX_INITIAL_CELL_COUNT_FOR_FITTING_RATES = true

R = 7
Time = 0:3:72;
Conc= [0, 31.25*10^-9, 62.5*10^-9, 125*10^-9, 250*10^-9, 375*10^-9, 500*10^-9, 1.25*10^-6, 2.5*10^-6, 3.75*10^-6, 5*10^-6];
N_t = length(Time) % 25
N_c = length(Conc) % 11
max_T = Time(N_t);
max_C = Conc(N_c);

% To navigate xslx file and csv we introduce these
first_letters = ["B", "C", "D", "E", "F", "G", "H"];
second_letters = ["I", "J", "K", "L", "M", "N", "O"]; 
letters = [first_letters, second_letters];
csv_conc_indexes_BF = [10:23,2:9];
csv_conc_indexes_GFP = [10:12,2:9];
num_col_csv_BF = length(csv_conc_indexes_BF);
num_col_csv_GFP = length(csv_conc_indexes_GFP);

% Bright field: All 4:
SENSITIVE_500_BF = zeros(R, N_c, N_t);
SENSITIVE_1000_BF = zeros(R, N_c, N_t);
RESISTANT_250_BF = zeros(R, N_c, N_t);
RESISTANT_500_BF = zeros(R, N_c, N_t);

% Green fluorescent protein: Only sensitives
SENSITIVE_500_GFP = zeros(R, N_c, N_t);
SENSITIVE_1000_GFP = zeros(R, N_c, N_t);

% Import csv
BF_data = readtable('Fwd__Results_of_Monoculture_experiment/Summary_BF.csv'); % Bright field
GFP_data = readtable('Fwd__Results_of_Monoculture_experiment/Summary_GFP.csv'); % Green fluorescent protein

% BF data
N_rows_BF = N_t * N_c * R * 4; % 25*11*7*4 rows in the csv file
for row_index=0:(N_rows_BF-1)
    time_index = mod(row_index, 25) + 1;
    % Indexes in xslx file
    char_index = fix(row_index / (2*N_c*N_t)); % 550 = 2*25*11, # columns equal 2 * # concentrations, # rows equal 2 * # replicates
    proxy_conc_index = mod( fix(row_index / N_t), num_col_csv_BF); % 
    
    % These match with the filenames in the csv: 
    %BF_data{row_index+1,1} % filename
    %letter = letters(char_index + 1)
    %csv_conc_index = csv_conc_indexes_BF(proxy_conc_index + 1)

    conc_index_22 = csv_conc_indexes_BF(proxy_conc_index+1) - 2;
    conc_index_11 = mod(conc_index_22, 11) + 1;

    r_index = mod(char_index, 7) + 1;
    %[row_index, proxy_conc_index, conc_index_22, conc_index_11] 

    if (char_index+1) <= 7 & (conc_index_22+1) <= 11
        SENSITIVE_500_BF(r_index, conc_index_11, time_index) = BF_data{row_index+1,2}; % Sensitive cells, 500 initial

    elseif (char_index+1) > 7 & (conc_index_22+1) <= 11
        SENSITIVE_1000_BF(r_index, conc_index_11, time_index) = BF_data{row_index+1,2}; % Sensitive cells, 1000 initial

    elseif (char_index+1) <= 7 & (conc_index_22+1) > 11
        RESISTANT_250_BF(r_index, conc_index_11, time_index) = BF_data{row_index+1,2}; % Resistant cells, 250 initial

    elseif (char_index+1) > 7 & (conc_index_22+1) > 11
        RESISTANT_500_BF(r_index, conc_index_11, time_index) = BF_data{row_index+1,2}; % Resistant cells, 500 initial

    end
end

% GFP data: Only for sensitive? 
N_rows_GFP = N_t * N_c * R * 2; % Only sensitive
for row_index=0:(N_rows_GFP-1)
    time_index = mod(row_index, 25) + 1;
    % Indexes in xslx file
    char_index = fix(row_index / (N_c*N_t)); % # columns equal # concentrations, # rows equal 2 * # replicates
    proxy_conc_index = mod( fix(row_index / N_t), num_col_csv_GFP); % 

    % These match with the filenames in the csv: 
    %GFP_data{row_index+1,1}; % filename
    %letter = letters(char_index + 1);
    %csv_conc_index = csv_conc_indexes_GFP(proxy_conc_index + 1);

    conc_index_22 = csv_conc_indexes_GFP(proxy_conc_index+1) - 2;
    conc_index_11 = mod(conc_index_22, 11) + 1;

    r_index = mod(char_index, 7) + 1;

    if (char_index+1) <= 7
        SENSITIVE_500_GFP(r_index, conc_index_11, time_index) = GFP_data{row_index+1,2}; % Sensitive cells, 500 initial

    elseif (char_index+1) > 7
        SENSITIVE_1000_GFP(r_index, conc_index_11, time_index) = GFP_data{row_index+1,2}; % Sensitive cells, 1000 initial

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outlier removal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Censor obvious measurement errors: i) cells that seem to rise from the dead and ii) RESISTANT_250_BF, zero concentration:
% 2) Take the logarithm to transform exponential data to linear data. Find the slope for each replicate and concentration individually (rate in exponential growth) 
%    and calculate residuals w.r.t this line
% 3) Remove outliers based on residuals using Matlab function "isoutlier" with interquartile range = 2.
% 4) Remove times {0,3,6} for all data, and time {72} for resistant 500 due to violation of exponential growth assumption
% For the inference, replicates with any NaN values are removed.

% Censor obvious measurement errors: i) cells that seem to rise from the dead and ii) RESISTANT_250_BF, zero concentration:
SENSITIVE_500_GFP(:,9,14:25) = NaN;
SENSITIVE_500_GFP(:,10,11:25) = NaN;
SENSITIVE_500_GFP(:,11,10:25) = NaN;

SENSITIVE_1000_GFP(:,9:11,12:25) = NaN;

RESISTANT_250_BF(:,1,8:14) = NaN;

% Take the logarithm to transform exponential data to linear data
logSENSITIVE_500_BF = log(SENSITIVE_500_BF);
logSENSITIVE_1000_BF = log(SENSITIVE_1000_BF);
logRESISTANT_250_BF = log(RESISTANT_250_BF);
logRESISTANT_500_BF = log(RESISTANT_500_BF);
logSENSITIVE_500_GFP = log(SENSITIVE_500_GFP);
logSENSITIVE_1000_GFP = log(SENSITIVE_1000_GFP);

% Find the slope of each replicate for each concentration (rate in exponential growth) and calculate residuals w.r.t this line
% Remove outliers using Matlab function "isoutlier" on the residuals
x = ones(N_t,2);
x(:,2) = 3*(0:(N_t-1))';

OUTLIER_BOOL = zeros(R, N_c, N_t);
% For each replicate evolving in time we want to remove measurement erros. 
for rep_index=1:R
    for conc_index=1:N_c
        squeeze(SENSITIVE_500_GFP(rep_index,conc_index,:));
        % Fit a linear model to the data
        y = squeeze(logSENSITIVE_500_GFP(rep_index,conc_index,:));
        coef_beta = x\y;
        ypred = x*coef_beta;

        %figure
        %scatter(x(:,2),y)
        %hold on
        %plot(x(:,2), ypred, 'r')
        %ylim([0 inf])
        %grid on
        %hold off
        %saveas(gcf,[pwd strcat('/plots/dagim_monoclonal/logMEAN_concentration_', int2str(conc_index), '.png')])

        y_residuals = y - ypred;
        TF = isoutlier(y_residuals,'quartiles', 'ThresholdFactor', 2);
        %TF = isoutlier(y_residuals, 'grubbs') %, 'movmedian', 9, 'ThresholdFactor', 3)
        OUTLIER_BOOL(rep_index,conc_index,:) = TF;
    end
end
% Apply boolean outlier filter
SENSITIVE_500_GFP(OUTLIER_BOOL == 1) = NaN;

OUTLIER_BOOL = zeros(R, N_c, N_t);
% For each replicate evolving in time we want to remove measurement erros. 
for rep_index=1:R
    for conc_index=1:N_c
        squeeze(SENSITIVE_1000_GFP(rep_index,conc_index,:));
        % Fit a linear model to the data
        y = squeeze(logSENSITIVE_1000_GFP(rep_index,conc_index,:));
        coef_beta = x\y;
        ypred = x*coef_beta;

        %figure
        %scatter(x(:,2),y)
        %hold on
        %plot(x(:,2), ypred, 'r')
        %ylim([0 inf])
        %grid on
        %hold off
        %saveas(gcf,[pwd strcat('/plots/dagim_monoclonal/logMEAN_concentration_', int2str(conc_index), '.png')])

        y_residuals = y - ypred;
        TF = isoutlier(y_residuals,'quartiles', 'ThresholdFactor', 2);
        %TF = isoutlier(y_residuals, 'grubbs') %, 'movmedian', 9, 'ThresholdFactor', 3)
        OUTLIER_BOOL(rep_index,conc_index,:) = TF;
    end
end
% Apply boolean outlier filter
SENSITIVE_1000_GFP(OUTLIER_BOOL == 1) = NaN;

OUTLIER_BOOL = zeros(R, N_c, N_t);
% For each replicate evolving in time we want to remove measurement erros. 
for rep_index=1:R
    for conc_index=1:N_c
        squeeze(SENSITIVE_500_BF(rep_index,conc_index,:));
        % Fit a linear model to the data
        y = squeeze(logSENSITIVE_500_BF(rep_index,conc_index,:));
        coef_beta = x\y;
        ypred = x*coef_beta;

        %figure
        %scatter(x(:,2),y)
        %hold on
        %plot(x(:,2), ypred, 'r')
        %ylim([0 inf])
        %grid on
        %hold off
        %saveas(gcf,[pwd strcat('/plots/dagim_monoclonal/logMEAN_concentration_', int2str(conc_index), '.png')])

        y_residuals = y - ypred;
        TF = isoutlier(y_residuals,'quartiles', 'ThresholdFactor', 2);
        %TF = isoutlier(y_residuals, 'grubbs') %, 'movmedian', 9, 'ThresholdFactor', 3)
        OUTLIER_BOOL(rep_index,conc_index,:) = TF;
    end
end
% Apply boolean outlier filter
SENSITIVE_500_BF(OUTLIER_BOOL == 1) = NaN;

OUTLIER_BOOL = zeros(R, N_c, N_t);
% For each replicate evolving in time we want to remove measurement erros. 
for rep_index=1:R
    for conc_index=1:N_c
        squeeze(SENSITIVE_1000_BF(rep_index,conc_index,:));
        % Fit a linear model to the data
        y = squeeze(logSENSITIVE_1000_BF(rep_index,conc_index,:));
        coef_beta = x\y;
        ypred = x*coef_beta;

        %figure
        %scatter(x(:,2),y)
        %hold on
        %plot(x(:,2), ypred, 'r')
        %ylim([0 inf])
        %grid on
        %hold off
        %saveas(gcf,[pwd strcat('/plots/dagim_monoclonal/logMEAN_concentration_', int2str(conc_index), '.png')])

        y_residuals = y - ypred;
        TF = isoutlier(y_residuals,'quartiles', 'ThresholdFactor', 2);
        %TF = isoutlier(y_residuals, 'grubbs') %, 'movmedian', 9, 'ThresholdFactor', 3)
        OUTLIER_BOOL(rep_index,conc_index,:) = TF;
    end
end
% Apply boolean outlier filter
SENSITIVE_1000_BF(OUTLIER_BOOL == 1) = NaN;


OUTLIER_BOOL = zeros(R, N_c, N_t);
% For each replicate evolving in time we want to remove measurement erros. 
for rep_index=1:R
    for conc_index=1:N_c
        squeeze(RESISTANT_250_BF(rep_index,conc_index,:));
        % Fit a linear model to the data
        y = squeeze(logRESISTANT_250_BF(rep_index,conc_index,:));
        coef_beta = x\y;
        ypred = x*coef_beta;

        %figure
        %scatter(x(:,2),y)
        %hold on
        %plot(x(:,2), ypred, 'r')
        %ylim([0 inf])
        %grid on
        %hold off
        %saveas(gcf,[pwd strcat('/plots/dagim_monoclonal/logMEAN_concentration_', int2str(conc_index), '.png')])

        y_residuals = y - ypred;
        TF = isoutlier(y_residuals,'quartiles', 'ThresholdFactor', 2);
        %TF = isoutlier(y_residuals, 'grubbs') %, 'movmedian', 9, 'ThresholdFactor', 3)
        OUTLIER_BOOL(rep_index,conc_index,:) = TF;
    end
end
% Apply boolean outlier filter
RESISTANT_250_BF(OUTLIER_BOOL == 1) = NaN;

OUTLIER_BOOL = zeros(R, N_c, N_t);
% For each replicate evolving in time we want to remove measurement erros. 
for rep_index=1:R
    for conc_index=1:N_c
        squeeze(RESISTANT_500_BF(rep_index,conc_index,:));
        % Fit a linear model to the data
        y = squeeze(logRESISTANT_500_BF(rep_index,conc_index,:));
        coef_beta = x\y;
        ypred = x*coef_beta;

        %figure
        %scatter(x(:,2),y)
        %hold on
        %plot(x(:,2), ypred, 'r')
        %ylim([0 inf])
        %grid on
        %hold off
        %saveas(gcf,[pwd strcat('/plots/dagim_monoclonal/logMEAN_concentration_', int2str(conc_index), '.png')])

        y_residuals = y - ypred;
        TF = isoutlier(y_residuals,'quartiles', 'ThresholdFactor', 2);
        %TF = isoutlier(y_residuals, 'grubbs') %, 'movmedian', 9, 'ThresholdFactor', 3)
        OUTLIER_BOOL(rep_index,conc_index,:) = TF;
    end
end
% Apply boolean outlier filter
RESISTANT_500_BF(OUTLIER_BOOL == 1) = NaN;

% Remove times {0,3,6}
SENSITIVE_500_BF = SENSITIVE_500_BF(:,:,4:N_t);
SENSITIVE_1000_BF = SENSITIVE_1000_BF(:,:,4:N_t);
RESISTANT_250_BF = RESISTANT_250_BF(:,:,4:N_t);
RESISTANT_500_BF = RESISTANT_500_BF(:,:,4:N_t);
SENSITIVE_500_GFP = SENSITIVE_500_GFP(:,:,4:N_t);
SENSITIVE_1000_GFP = SENSITIVE_1000_GFP(:,:,4:N_t);
Time = Time(4:N_t);
N_t = N_t - 3;
x = ones(N_t,2);
x(:,2) = 3*(0:(N_t-1))';

% For resistant 500 remove time {72}
N_t_RES_500 = N_t - 1;
Time_res500bf = Time(1:N_t_RES_500);
RESISTANT_500_BF = RESISTANT_500_BF(:,:,1:N_t_RES_500);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recompute log values
logSENSITIVE_500_BF = log(SENSITIVE_500_BF);
logSENSITIVE_1000_BF = log(SENSITIVE_1000_BF);
logRESISTANT_250_BF = log(RESISTANT_250_BF);
logRESISTANT_500_BF = log(RESISTANT_500_BF);
logSENSITIVE_500_GFP = log(SENSITIVE_500_GFP);
logSENSITIVE_1000_GFP = log(SENSITIVE_1000_GFP);

%% Checking if we need to remove the last timepoints
%yar = 7;
%for rep_index=yar:yar %R
%    for conc_index=1:N_c
%        squeeze(SENSITIVE_500_BF(rep_index,conc_index,:))
%        % Fit a linear model to the data
%        y = squeeze(logSENSITIVE_500_BF(rep_index,conc_index,:));
%        coef_beta = x\y;
%        ypred = x*coef_beta;
%
%        figure
%        scatter(x(:,2),y)
%        hold on
%        plot(x(:,2), ypred, 'r')
%        ylim([0 inf])
%        grid on
%        hold off
%        saveas(gcf,[pwd strcat('/plots/dagim_monoclonal/logMEAN_concentration_', int2str(conc_index), '.png')])
%
%        y_residuals = y - ypred;
%        TF = isoutlier(y_residuals,'quartiles', 'ThresholdFactor', 2)
%        %TF = isoutlier(y_residuals, 'grubbs') %, 'movmedian', 9, 'ThresholdFactor', 3)
%        %OUTLIER_BOOL(rep_index,conc_index,:) = TF;
%    end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average data and find std
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MEAN_SENSITIVE_500_BF = squeeze(mean(SENSITIVE_500_BF, 1, 'omitnan'));
MEAN_SENSITIVE_1000_BF = squeeze(mean(SENSITIVE_1000_BF, 1, 'omitnan'));
MEAN_RESISTANT_250_BF = squeeze(mean(RESISTANT_250_BF, 1, 'omitnan'));
MEAN_RESISTANT_500_BF = squeeze(mean(RESISTANT_500_BF, 1, 'omitnan'));
MEAN_SENSITIVE_500_GFP = squeeze(mean(SENSITIVE_500_GFP, 1, 'omitnan'));
MEAN_SENSITIVE_1000_GFP = squeeze(mean(SENSITIVE_1000_GFP, 1, 'omitnan'));

STD_SENSITIVE_500_BF = squeeze(std(SENSITIVE_500_BF, 1, 'omitnan'));
STD_SENSITIVE_1000_BF = squeeze(std(SENSITIVE_1000_BF, 1, 'omitnan'));
STD_RESISTANT_250_BF = squeeze(std(RESISTANT_250_BF, 1, 'omitnan'));
STD_RESISTANT_500_BF = squeeze(std(RESISTANT_500_BF, 1, 'omitnan'));
STD_SENSITIVE_500_GFP = squeeze(std(SENSITIVE_500_GFP, 1, 'omitnan'));
STD_SENSITIVE_1000_GFP = squeeze(std(SENSITIVE_1000_GFP, 1, 'omitnan'));

% Compute average and min, max std across times and concentrations
average_STD_SENSITIVE_500_BF = mean(STD_SENSITIVE_500_BF, 'all', 'omitnan');
average_STD_SENSITIVE_1000_BF = mean(STD_SENSITIVE_1000_BF, 'all', 'omitnan');
average_STD_RESISTANT_250_BF = mean(STD_RESISTANT_250_BF, 'all', 'omitnan');
average_STD_RESISTANT_500_BF = mean(STD_RESISTANT_500_BF, 'all', 'omitnan');
average_STD_SENSITIVE_500_GFP = mean(STD_SENSITIVE_500_GFP, 'all', 'omitnan');
average_STD_SENSITIVE_1000_GFP = mean(STD_SENSITIVE_1000_GFP, 'all', 'omitnan');

min_STD_SENSITIVE_500_BF = min(STD_SENSITIVE_500_BF, [], 'all');
min_STD_SENSITIVE_1000_BF = min(STD_SENSITIVE_1000_BF, [], 'all');
min_STD_RESISTANT_250_BF = min(STD_RESISTANT_250_BF, [], 'all');
min_STD_RESISTANT_500_BF = min(STD_RESISTANT_500_BF, [], 'all');
min_STD_SENSITIVE_500_GFP = min(STD_SENSITIVE_500_GFP, [], 'all');
min_STD_SENSITIVE_1000_GFP = min(STD_SENSITIVE_1000_GFP, [], 'all');

max_STD_SENSITIVE_500_BF = max(STD_SENSITIVE_500_BF, [], 'all');
max_STD_SENSITIVE_1000_BF = max(STD_SENSITIVE_1000_BF, [], 'all');
max_STD_RESISTANT_250_BF = max(STD_RESISTANT_250_BF, [], 'all');
max_STD_RESISTANT_500_BF = max(STD_RESISTANT_500_BF, [], 'all');
max_STD_SENSITIVE_500_GFP = max(STD_SENSITIVE_500_GFP, [], 'all');
max_STD_SENSITIVE_1000_GFP = max(STD_SENSITIVE_1000_GFP, [], 'all');

logMEAN_SENSITIVE_500_BF = log(MEAN_SENSITIVE_500_BF);
logMEAN_SENSITIVE_1000_BF = log(MEAN_SENSITIVE_1000_BF);
logMEAN_RESISTANT_250_BF = log(MEAN_RESISTANT_250_BF);
logMEAN_RESISTANT_500_BF = log(MEAN_RESISTANT_500_BF);
logMEAN_SENSITIVE_500_GFP = log(MEAN_SENSITIVE_500_GFP);
logMEAN_SENSITIVE_1000_GFP = log(MEAN_SENSITIVE_1000_GFP);

logSTD_SENSITIVE_500_BF = log(STD_SENSITIVE_500_BF);
logSTD_SENSITIVE_1000_BF = log(STD_SENSITIVE_1000_BF);
logSTD_RESISTANT_250_BF = log(STD_RESISTANT_250_BF);
logSTD_RESISTANT_500_BF = log(STD_RESISTANT_500_BF);
logSTD_SENSITIVE_500_GFP = log(STD_SENSITIVE_500_GFP);
logSTD_SENSITIVE_1000_GFP = log(STD_SENSITIVE_1000_GFP);

% Find fitted growth rate for each replicate, then average
% In cases where the initial y value is NaN, we skip the replicate
fitted_rates_sens500bf_per_replicate = zeros(R,N_c);
for rep_index=1:R
    for conc_index=1:N_c
        % Fit a linear model to the data
        y = squeeze(logSENSITIVE_500_BF(rep_index,conc_index,:));
        if FIX_INITIAL_CELL_COUNT_FOR_FITTING_RATES
            fitted_rates_sens500bf_per_replicate(rep_index, conc_index) = fit_with_fixed_point(x,y,N_t);
        else
            coef_beta = x\y;
            fitted_rates_sens500bf_per_replicate(rep_index, conc_index) = coef_beta(2);
        end
    end
end

fitted_rates_sens1000bf_per_replicate = zeros(R,N_c);
for rep_index=1:R
    for conc_index=1:N_c
        % Fit a linear model to the data
        y = squeeze(logSENSITIVE_1000_BF(rep_index,conc_index,:));
        if FIX_INITIAL_CELL_COUNT_FOR_FITTING_RATES
            fitted_rates_sens1000bf_per_replicate(rep_index, conc_index) = fit_with_fixed_point(x,y,N_t);
        else
            coef_beta = x\y;
            fitted_rates_sens1000bf_per_replicate(rep_index, conc_index) = coef_beta(2);
        end
    end
end

fitted_rates_res250bf_per_replicate = zeros(R,N_c);
for rep_index=1:R
    for conc_index=1:N_c
        % Fit a linear model to the data
        y = squeeze(logRESISTANT_250_BF(rep_index,conc_index,:));
        x_res_250 = x;
        if conc_index == 1
            % Remove data that is corrupted across all replicates
            y = [y(1:4)', y(12:N_t)']';
            x_res_250 = [x(1:4,:)', x(12:N_t,:)']';
        end
        if FIX_INITIAL_CELL_COUNT_FOR_FITTING_RATES
            fitted_rates_res250bf_per_replicate(rep_index, conc_index) = fit_with_fixed_point(x_res_250,y,N_t);
        else
            coef_beta = x_res_250\y;
            fitted_rates_res250bf_per_replicate(rep_index, conc_index) = coef_beta(2);
        end
    end
end

fitted_rates_res500bf_per_replicate = zeros(R,N_c);
for rep_index=1:R
    for conc_index=1:N_c
        % Fit a linear model to the data
        y = squeeze(logRESISTANT_500_BF(rep_index,conc_index,:));
        if FIX_INITIAL_CELL_COUNT_FOR_FITTING_RATES
            fitted_rates_res500bf_per_replicate(rep_index, conc_index) = fit_with_fixed_point(x,y,N_t);
        else
            coef_beta = x(1:N_t-1,:)\y;
            fitted_rates_res500bf_per_replicate(rep_index, conc_index) = coef_beta(2);
        end
    end
end


fitted_rates_sens500gfp_per_replicate = zeros(R,N_c);
for rep_index=1:R
    for conc_index=1:N_c
        % Fit a linear model to the data
        y = squeeze(logSENSITIVE_500_GFP(rep_index,conc_index,:));
        if FIX_INITIAL_CELL_COUNT_FOR_FITTING_RATES
            fitted_rates_sens500gfp_per_replicate(rep_index, conc_index) = fit_with_fixed_point(x,y,N_t);
        else
            coef_beta = x\y;
            fitted_rates_sens500gfp_per_replicate(rep_index, conc_index) = coef_beta(2);
        end
    end
end

fitted_rates_sens1000gfp_per_replicate = zeros(R,N_c);
for rep_index=1:R
    for conc_index=1:N_c
        % Fit a linear model to the data
        y = squeeze(logSENSITIVE_1000_GFP(rep_index,conc_index,:));
        if FIX_INITIAL_CELL_COUNT_FOR_FITTING_RATES
            fitted_rates_sens1000gfp_per_replicate(rep_index, conc_index) = fit_with_fixed_point(x,y,N_t);
        else
            coef_beta = x\y;
            fitted_rates_sens1000gfp_per_replicate(rep_index, conc_index) = coef_beta(2);
        end
    end
end

% Find experimental growth rate by fitting a line to the mean
x = ones(N_t,2);
x(:,2) = 3*(0:(N_t-1))';

fitted_rates_sens500bf = zeros(1,N_c);
for i=1:N_c
    y = squeeze(logMEAN_SENSITIVE_500_BF(i,:))';
    coef_beta = x\y;
    fitted_rates_sens500bf(i) = coef_beta(2);
end
fitted_rates_sens1000bf = zeros(1,N_c);
for i=1:N_c
    y = squeeze(logMEAN_SENSITIVE_1000_BF(i,:))';
    coef_beta = x\y;
    fitted_rates_sens1000bf(i) = coef_beta(2);
end

fitted_rates_res250bf = zeros(1,N_c);
for i=1:N_c
    y = squeeze(logMEAN_RESISTANT_250_BF(i,:))';
    x_res_250 = x;
    if i == 1
        % Remove data that is corrupted across all replicates
        y = [y(1:4)', y(12:N_t)']';
        x_res_250 = [x(1:4,:)', x(12:N_t,:)']';
    end
    coef_beta = x_res_250\y;
    fitted_rates_res250bf(i) = coef_beta(2);
end

fitted_rates_res500bf = zeros(1,N_c);
for i=1:N_c
    y = squeeze(logMEAN_RESISTANT_500_BF(i,:))';
    coef_beta = x(1:N_t_RES_500,:)\y;
    fitted_rates_res500bf(i) = coef_beta(2);
end
fitted_rates_sens500gfp = zeros(1,N_c);
for i=1:N_c
    y = squeeze(logMEAN_SENSITIVE_500_GFP(i,:))';
    coef_beta = x\y;
    fitted_rates_sens500gfp(i) = coef_beta(2);
end
fitted_rates_sens1000gfp = zeros(1,N_c);
for i=1:N_c
    y = squeeze(logMEAN_SENSITIVE_1000_GFP(i,:))';
    coef_beta = x\y;
    fitted_rates_sens1000gfp(i) = coef_beta(2);
end
% save data
save('./data/dagim_monoclonal/CLEAN_DATA.mat', 'SENSITIVE_500_BF', 'SENSITIVE_1000_BF', 'RESISTANT_250_BF', 'RESISTANT_500_BF', 'SENSITIVE_500_GFP', 'SENSITIVE_1000_GFP', 'MEAN_SENSITIVE_500_BF', 'MEAN_SENSITIVE_1000_BF', 'MEAN_RESISTANT_250_BF', 'MEAN_RESISTANT_500_BF', 'MEAN_SENSITIVE_500_GFP', 'MEAN_SENSITIVE_1000_GFP', 'STD_SENSITIVE_500_BF', 'STD_SENSITIVE_1000_BF', 'STD_RESISTANT_250_BF', 'STD_RESISTANT_500_BF', 'STD_SENSITIVE_500_GFP', 'STD_SENSITIVE_1000_GFP')
save('./data/dagim_monoclonal/monoclonal_fitted_rates.mat', 'fitted_rates_sens500bf', 'fitted_rates_sens1000bf', 'fitted_rates_res250bf', 'fitted_rates_res500bf', 'fitted_rates_sens500gfp', 'fitted_rates_sens1000gfp', 'fitted_rates_sens500bf_per_replicate', 'fitted_rates_sens1000bf_per_replicate', 'fitted_rates_res250bf_per_replicate', 'fitted_rates_res500bf_per_replicate', 'fitted_rates_sens500gfp_per_replicate', 'fitted_rates_sens1000gfp_per_replicate')

%% Plot mean logarithm values
%figure
%for c_index=1:N_c
%    plot(Time, logMEAN_SENSITIVE_500_GFP(c_index,:), '-x')
%    %errorbar(Time, logMEAN_SENSITIVE_500_GFP(c_index,:), logSTD_SENSITIVE_500_GFP(c_index, :))
%    if c_index == 1
%        hold on
%    end
%end
%ylabel('Cell count')
%xlabel('Time (h)')
%ylim([0 inf])
%%legend(conclabels, 'location', 'northwest')
%title("Sensitive 500 GFP (averaged)")
%saveas(gcf, "plots/dagim_monoclonal/mono_GFP_sensitive_500.png")
%
%figure
%for c_index=1:N_c
%    plot(Time, logMEAN_SENSITIVE_1000_GFP(c_index,:), '-x')
%    %errorbar(Time, logMEAN_SENSITIVE_500_GFP(c_index,:), logSTD_SENSITIVE_500_GFP(c_index, :))
%    if c_index == 1
%        hold on
%    end
%end
%ylabel('Cell count')
%xlabel('Time (h)')
%ylim([0 inf])
%%legend(conclabels, 'location', 'northwest')
%title("Sensitive 1000 GFP (averaged)")
%saveas(gcf, "plots/dagim_monoclonal/mono_GFP_sensitive_1000.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Infer and plot estimated cell lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if INFERENCE
% xx = [p, alpha1, coef_beta, E1, n1, alpha2, b2, E2, n2, sig]
lb=zeros(5,1);
ub=1000*ones(length(lb),1);
%%%% ub(1)=0.5; Mixture parameter excluded for mono culture
ub(1)=0.1;ub(2)=1;ub(3)=50;ub(4)=10; % Hill parameters for 1st pop
ub(5)=2000; % Noise
%%%%%%ub(6)=0.1;ub(7)=1;ub(8)=50;ub(9)=10; % Hill parameters for 2nd pop
%%%%%%ub(10)=5000; % Noise
num_optim = 200;
options = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',5990,'Display','off');

DATA = SENSITIVE_500_BF;
f=@(x)mono_vectorized_objective_function(x,DATA,Conc,Time,R);
fval=inf;
clear xx;
clear ff;
clear x_final;
rng(42);
for n=1:num_optim
    x0=rand(length(ub),1).*(ub-lb) + lb;
    [xx,ff]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
    if ff<fval
        x_final=xx;
        fval=ff;
    end    
end
x_final_sensitive_500_bf = x_final

DATA = SENSITIVE_1000_BF;
f=@(x)mono_vectorized_objective_function(x,DATA,Conc,Time,R);
fval=inf;
clear xx;
clear ff;
clear x_final;
rng(42);
for n=1:num_optim
    x0=rand(length(ub),1).*(ub-lb) + lb;
    [xx,ff]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
    if ff<fval
        x_final=xx;
        fval=ff;
    end    
end
x_final_sensitive_1000_bf = x_final

Conc_res_250 = [31.25*10^-9, 62.5*10^-9, 125*10^-9, 250*10^-9, 375*10^-9, 500*10^-9, 1.25*10^-6, 2.5*10^-6, 3.75*10^-6, 5*10^-6];
DATA = RESISTANT_250_BF(:,2:11,:); % We remove the zero drug data because of NaN values
f=@(x)mono_vectorized_objective_function(x,DATA,Conc_res_250,Time,R);
fval=inf;
clear xx;
clear ff;
clear x_final;
rng(42);
for n=1:num_optim
    x0=rand(length(ub),1).*(ub-lb) + lb;
    [xx,ff]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
    if ff<fval
        x_final=xx;
        fval=ff;
    end    
end
x_final_resistant_250_bf = x_final

DATA = RESISTANT_500_BF;
f=@(x)mono_vectorized_objective_function(x,DATA,Conc,Time_res500bf,R);
fval=inf;
clear xx;
clear ff;
clear x_final;
rng(42);
for n=1:num_optim
    x0=rand(length(ub),1).*(ub-lb) + lb;
    [xx,ff]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
    if ff<fval
        x_final=xx;
        fval=ff;
    end    
end
x_final_resistant_500_bf = x_final

DATA = SENSITIVE_500_GFP(2:7,:,:);
R_sens_500_gfp = 6;
f=@(x)mono_vectorized_objective_function(x,DATA,Conc,Time,R_sens_500_gfp);
fval=inf;
clear xx;
clear ff;
clear x_final;
rng(42);
for n=1:num_optim
    x0=rand(length(ub),1).*(ub-lb) + lb;
    [xx,ff]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
    if ff<fval
        x_final=xx;
        fval=ff;
    end    
end
x_final_sensitive_500_gfp = x_final

R_sens_1000_gfp = 6;
Conc_sens_1000_gfp = [0, 31.25*10^-9, 62.5*10^-9, 125*10^-9, 250*10^-9, 375*10^-9, 500*10^-9, 1.25*10^-6, 2.5*10^-6, 3.75*10^-6];;
DATA = SENSITIVE_1000_GFP(1:6,1:10,:);
f=@(x)mono_vectorized_objective_function(x,DATA,Conc_sens_1000_gfp,Time,R_sens_1000_gfp);
fval=inf;
clear xx;
clear ff;
clear x_final;
rng(42);
for n=1:num_optim
    x0=rand(length(ub),1).*(ub-lb) + lb;
    [xx,ff]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
    if ff<fval
        x_final=xx;
        fval=ff;
    end    
end
x_final_sensitive_1000_gfp = x_final

save('./data/dagim_monoclonal/monoclonal_x_finals.mat', 'x_final_sensitive_500_bf', 'x_final_sensitive_1000_bf', 'x_final_resistant_250_bf', 'x_final_resistant_500_bf', 'x_final_sensitive_500_gfp', 'x_final_sensitive_1000_gfp')
end % if INFERENCE loop

if PLOT_FITTED_VS_ESTIMATED_RATES
load('./data/dagim_monoclonal/monoclonal_x_finals.mat')
% Plot rates as function of concentration
x = zeros(1,1000);
x(2:1000) = logspace(-8,-4.86,999);
% Experiment concentrations tested: 
Conc = [0, 31.25*10^-9, 62.5*10^-9, 125*10^-9, 250*10^-9, 375*10^-9, 500*10^-9, 1.25*10^-6, 2.5*10^-6, 3.75*10^-6, 5*10^-6]; 

inferred_rate_sens_500_bf = ratefunc(x_final_sensitive_500_bf,x);
inferred_rate_sens_1000_bf = ratefunc(x_final_sensitive_1000_bf,x);
inferred_rate_res_250_bf = ratefunc(x_final_resistant_250_bf,x);
inferred_rate_res_500_bf = ratefunc(x_final_resistant_500_bf,x);
inferred_rate_sens_500_gfp = ratefunc(x_final_sensitive_500_gfp,x);
inferred_rate_sens_1000_gfp = ratefunc(x_final_sensitive_1000_gfp,x);

newcolors0 = [
    0 0 0.726562
    0 0 0.726562
    0.477504 0.821444 0.318195
    0.477504 0.821444 0.318195
    0.792968 0 0
    0.792968 0 0
];
figure
colororder(newcolors0)
semilogx(x,inferred_rate_sens_500_bf, '--', 'linewidth', 2) %, '-o')
hold on 
semilogx(x,inferred_rate_sens_1000_bf, '-', 'linewidth', 2) %, '-+')
semilogx(x,inferred_rate_sens_500_gfp, '--', 'linewidth', 2) %, '-^')
semilogx(x,inferred_rate_sens_1000_gfp, '-', 'linewidth', 2) %, '-p')
semilogx(x,inferred_rate_res_250_bf, '--', 'linewidth', 2) %, '->')
semilogx(x,inferred_rate_res_500_bf, '-', 'linewidth', 2) %, '-v')
ylabel('Growth rate')
xlabel('Concentration')
ylim([-0.25, 0.1])
xlim([10^(-8) 10^(-5.2)])
for i=2:length(Conc)
    line([Conc(i), Conc(i)], [-0.25, 0.1], 'Color', 'k');
end
% Fitted to means
semilogx(Conc(2:N_c), fitted_rates_sens500bf(2:N_c), 'xk-')
semilogx(Conc(2:N_c), fitted_rates_sens1000bf(2:N_c), 'xk-')
semilogx(Conc(2:N_c), fitted_rates_res250bf(2:N_c), 'xk-')
semilogx(Conc(2:N_c), fitted_rates_res500bf(2:N_c), 'xk-')
semilogx(Conc(2:N_c), fitted_rates_sens500gfp(2:N_c), 'xk-')
semilogx(Conc(2:N_c), fitted_rates_sens1000gfp(2:N_c), 'xk-')

title("Inferred growth rates")
legend('Sensitive 500 BF', 'Sensitive 1000 BF', 'Sensitive 500 GFP', 'Sensitive 1000 GFP', 'Resistant 250 BF', 'Resistant 500 BF', 'location', 'southwest')
saveas(gcf, "plots/dagim_monoclonal/inferred_and_fitted_rates_from_plot_monoclonal.png")

% Plot for each dataset separately
% Sens 500 BF
singlecolor = [
    0 0 0.726562
];
figure
colororder(singlecolor)
semilogx(x,inferred_rate_sens_500_bf, '--', 'linewidth', 2) %, '-o')
hold on 
ylabel('Growth rate')
xlabel('Concentration')
ylim([-0.25, 0.1])
xlim([10^(-8) 10^(-5.2)])
for i=2:length(Conc)
    line([Conc(i), Conc(i)], [-0.25, 0.1], 'Color', 'k');
end
% Fitted to means
for r_index=1:R
    semilogx(Conc(2:N_c), fitted_rates_sens500bf_per_replicate(r_index,2:N_c), 'xk-')
end
mean_rate_across_replicates = squeeze(mean(fitted_rates_sens500bf_per_replicate(:,2:N_c), 'omitnan'))
semilogx(Conc(2:N_c), mean_rate_across_replicates, 'xr-') %, 'Linewidth', 2)
title("Inferred growth rates")
legend('Sensitive 500 BF') %, 'Sensitive 1000 BF', 'Sensitive 500 GFP', 'Sensitive 1000 GFP', 'Resistant 250 BF', 'Resistant 500 BF', 'location', 'southwest')
saveas(gcf, "plots/dagim_monoclonal/inferred_and_fitted_rates_from_plot_monoclonal_per_replicate_sens500bf.png")

% Sens 1000 BF
singlecolor = [
    0 0 0.726562
];
figure
colororder(singlecolor)
semilogx(x,inferred_rate_sens_1000_bf, '-', 'linewidth', 2) %, '-+')
hold on 
ylabel('Growth rate')
xlabel('Concentration')
ylim([-0.25, 0.1])
xlim([10^(-8) 10^(-5.2)])
for i=2:length(Conc)
    line([Conc(i), Conc(i)], [-0.25, 0.1], 'Color', 'k');
end
% Fitted to means
for r_index=1:R
    semilogx(Conc(2:N_c), fitted_rates_sens1000bf_per_replicate(r_index,2:N_c), 'xk-')
end
mean_rate_across_replicates = squeeze(mean(fitted_rates_sens1000bf_per_replicate(:,2:N_c), 'omitnan'))
semilogx(Conc(2:N_c), mean_rate_across_replicates, 'xr-') %, 'Linewidth', 2)
title("Inferred growth rates")
legend('Sensitive 1000 BF') %, 'Sensitive 1000 BF', 'Sensitive 500 GFP', 'Sensitive 1000 GFP', 'Resistant 250 BF', 'Resistant 500 BF', 'location', 'southwest')
saveas(gcf, "plots/dagim_monoclonal/inferred_and_fitted_rates_from_plot_monoclonal_per_replicate_sens1000bf.png")

% Res 250 BF
singlecolor = [
    0.477504 0.821444 0.318195
];
figure
colororder(singlecolor)
semilogx(x,inferred_rate_res_250_bf, '--', 'linewidth', 2) %, '->')
hold on 
ylabel('Growth rate')
xlabel('Concentration')
ylim([-0.25, 0.1])
xlim([10^(-8) 10^(-5.2)])
for i=2:length(Conc)
    line([Conc(i), Conc(i)], [-0.25, 0.1], 'Color', 'k');
end
% Fitted to means
for r_index=1:R
    semilogx(Conc(2:N_c), fitted_rates_res250bf_per_replicate(r_index,2:N_c), 'xk-')
end
mean_rate_across_replicates = squeeze(mean(fitted_rates_res250bf_per_replicate(:,2:N_c), 'omitnan'))
semilogx(Conc(2:N_c), mean_rate_across_replicates, 'xr-') %, 'Linewidth', 2)
title("Inferred growth rates")
legend('Resistant 250 BF') %, 'Sensitive 1000 BF', 'Sensitive 500 GFP', 'Sensitive 1000 GFP', 'Resistant 250 BF', 'Resistant 500 BF', 'location', 'southwest')
saveas(gcf, "plots/dagim_monoclonal/inferred_and_fitted_rates_from_plot_monoclonal_per_replicate_res250bf.png")

% Res 500 BF
singlecolor = [
    0.477504 0.821444 0.318195
];
figure
colororder(singlecolor)
semilogx(x,inferred_rate_res_500_bf, '-', 'linewidth', 2) %, '-v')
hold on 
ylabel('Growth rate')
xlabel('Concentration')
ylim([-0.25, 0.1])
xlim([10^(-8) 10^(-5.2)])
for i=2:length(Conc)
    line([Conc(i), Conc(i)], [-0.25, 0.1], 'Color', 'k');
end
% Fitted to means
for r_index=1:R
    semilogx(Conc(2:N_c), fitted_rates_res500bf_per_replicate(r_index,2:N_c), 'xk-')
end
mean_rate_across_replicates = squeeze(mean(fitted_rates_res500bf_per_replicate(:,2:N_c), 'omitnan'))
semilogx(Conc(2:N_c), mean_rate_across_replicates, 'xr-') %, 'Linewidth', 2)
title("Inferred growth rates")
legend('Resistant 500 BF') %, 'Sensitive 1000 BF', 'Sensitive 500 GFP', 'Sensitive 1000 GFP', 'Resistant 250 BF', 'Resistant 500 BF', 'location', 'southwest')
saveas(gcf, "plots/dagim_monoclonal/inferred_and_fitted_rates_from_plot_monoclonal_per_replicate_res500bf.png")

% Sens 500 GFP
singlecolor = [
    0.792968 0 0
];
figure
colororder(singlecolor)
semilogx(x,inferred_rate_sens_500_gfp, '--', 'linewidth', 2) %, '-^')
hold on 
ylabel('Growth rate')
xlabel('Concentration')
ylim([-0.25, 0.1])
xlim([10^(-8) 10^(-5.2)])
for i=2:length(Conc)
    line([Conc(i), Conc(i)], [-0.25, 0.1], 'Color', 'k');
end
% Fitted to means
for r_index=1:R
    semilogx(Conc(2:N_c), fitted_rates_sens500gfp_per_replicate(r_index,2:N_c), 'xk-')
end
mean_rate_across_replicates = squeeze(mean(fitted_rates_sens500gfp_per_replicate(:,2:N_c), 'omitnan'))
semilogx(Conc(2:N_c), mean_rate_across_replicates, 'xr-') %, 'Linewidth', 2)
title("Inferred growth rates")
legend('Sensitive 500 GFP') %, 'Sensitive 1000 BF', 'Sensitive 500 GFP', 'Sensitive 1000 GFP', 'Resistant 250 BF', 'Resistant 500 BF', 'location', 'southwest')
saveas(gcf, "plots/dagim_monoclonal/inferred_and_fitted_rates_from_plot_monoclonal_per_replicate_sens500gfp.png")

% Sens 1000 GFP
singlecolor = [
    0.792968 0 0
];
figure
colororder(singlecolor)
semilogx(x,inferred_rate_sens_1000_gfp, '-', 'linewidth', 2) %, '-p')
hold on 
ylabel('Growth rate')
xlabel('Concentration')
ylim([-0.25, 0.1])
xlim([10^(-8) 10^(-5.2)])
for i=2:length(Conc)
    line([Conc(i), Conc(i)], [-0.25, 0.1], 'Color', 'k');
end
% Fitted to means
for r_index=1:R
    semilogx(Conc(2:N_c), fitted_rates_sens1000gfp_per_replicate(r_index,2:N_c), 'xk-')
end
mean_rate_across_replicates = squeeze(mean(fitted_rates_sens1000gfp_per_replicate(:,2:N_c), 'omitnan'))
semilogx(Conc(2:N_c), mean_rate_across_replicates, 'xr-') %, 'Linewidth', 2)
title("Inferred growth rates")
legend('Sensitive 1000 GFP') %, 'Sensitive 1000 BF', 'Sensitive 500 GFP', 'Sensitive 1000 GFP', 'Resistant 250 BF', 'Resistant 500 BF', 'location', 'southwest')
saveas(gcf, "plots/dagim_monoclonal/inferred_and_fitted_rates_from_plot_monoclonal_per_replicate_sens1000gfp.png")

end % if PLOT_FITTED_VS_ESTIMATED_RATES loop

if PLOT_CLEAN_DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting in time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tlabels = strings(1,N_t);
for t_index=1:N_t
    tlabels(t_index) = strcat("T=", int2str(Time(t_index)));
end

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

newcolors2 = [
0.546157 0.038954 0.64701
0.640959 0.116492 0.602065
0.723444 0.196158 0.538981
0.798216 0.280197 0.469538
0.85975 0.360588 0.406917
0.9128 0.444029 0.346251
0.95547 0.533093 0.28549
];

figure
colororder(newcolors2)
hAx=axes;
hAx.XScale='log'
for t_index=N_t:-4:N_t
    %semilogx(Conc, MEAN_SENSITIVE_500_BF(:,t_index)', '-x')
    hold on
    errorbar(Conc, MEAN_SENSITIVE_500_BF(:,t_index)', STD_SENSITIVE_500_BF(:,t_index))
end
ylabel('Cell count')
xlabel('Concentration (M)')
ylim([0 inf])
legend(tlabels(N_t-1:-4:N_t), 'location', 'northeast')
title("Sensitive 500 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/conc_axis_mono_BF_sensitive_500_one_time.png")
ewfefwqewqdewd

figure
colororder(newcolors2)
hAx=axes;
hAx.XScale='log'
for t_index=N_t:-4:1
    %semilogx(Conc, MEAN_SENSITIVE_500_BF(:,t_index)', '-x')
    hold on
    errorbar(Conc, MEAN_SENSITIVE_500_BF(:,t_index)', STD_SENSITIVE_500_BF(:,t_index))
end
ylabel('Cell count')
xlabel('Concentration (M)')
ylim([0 inf])
legend(tlabels(N_t-1:-4:1), 'location', 'northeast')
title("Sensitive 500 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/conc_axis_mono_BF_sensitive_500.png")

figure
colororder(newcolors2)
hAx=axes;
hAx.XScale='log'
for t_index=N_t:-4:1
    %semilogx(Conc, MEAN_SENSITIVE_1000_BF(:,t_index)', '-x')
    hold on
    errorbar(Conc, MEAN_SENSITIVE_1000_BF(:,t_index)', STD_SENSITIVE_1000_BF(:,t_index))
end
ylabel('Cell count')
xlabel('Concentration (M)')
ylim([0 inf])
legend(tlabels(N_t-1:-4:1), 'location', 'northeast')
title("Sensitive 1000 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/conc_axis_mono_BF_sensitive_1000.png")

figure
colororder(newcolors2)
hAx=axes;
hAx.XScale='log'
for t_index=N_t:-4:1
    %semilogx(Conc, MEAN_RESISTANT_250_BF(:,t_index)', '-x')
    hold on
    errorbar(Conc, MEAN_RESISTANT_250_BF(:,t_index)', STD_RESISTANT_250_BF(:,t_index))
end
ylabel('Cell count')
xlabel('Concentration (M)')
ylim([0 inf])
legend(tlabels(N_t-1:-4:1), 'location', 'northeast')
title("Resistant 250 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/conc_axis_mono_BF_resistant_250.png")

figure
colororder(newcolors2)
hAx=axes;
hAx.XScale='log'
for t_index=N_t_RES_500:-4:1
    %semilogx(Conc, MEAN_RESISTANT_500_BF(:,t_index)', '-x')
    hold on
    errorbar(Conc, MEAN_RESISTANT_500_BF(:,t_index)', STD_RESISTANT_500_BF(:,t_index))
end
ylabel('Cell count')
xlabel('Concentration (M)')
ylim([0 inf])
legend(tlabels(N_t-1:-4:1), 'location', 'northeast')
title("Resistant 500 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/conc_axis_mono_BF_resistant_500.png")

figure
colororder(newcolors2)
hAx=axes;
hAx.XScale='log'
for t_index=N_t:-4:1
    %semilogx(Conc, MEAN_SENSITIVE_500_GFP(:,t_index)', '-x')
    hold on
    errorbar(Conc, MEAN_SENSITIVE_500_GFP(:,t_index)', STD_SENSITIVE_500_BF(:,t_index))
end
ylabel('Cell count')
xlabel('Concentration (M)')
ylim([0 inf])
legend(tlabels(N_t-1:-4:1), 'location', 'northeast')
title("Sensitive 500 GFP (averaged)")
saveas(gcf, "plots/dagim_monoclonal/conc_axis_mono_GFP_sensitive_500.png")

figure
colororder(newcolors2)
hAx=axes;
hAx.XScale='log'
for t_index=N_t:-4:1
    %semilogx(Conc, MEAN_SENSITIVE_1000_GFP(:,t_index)', '-x')
    hold on
    errorbar(Conc, MEAN_SENSITIVE_1000_GFP(:,t_index)', STD_SENSITIVE_1000_BF(:,t_index))
end
ylabel('Cell count')
xlabel('Concentration (M)')
ylim([0 inf])
legend(tlabels(N_t-1:-4:1), 'location', 'northeast')
title("Sensitive 1000 GFP (averaged)")
saveas(gcf, "plots/dagim_monoclonal/conc_axis_mono_GFP_sensitive_1000.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time on x axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time x axis
figure
colororder(newcolors)
for c_index=1:N_c
    plot(Time, MEAN_SENSITIVE_500_BF(c_index,:), '-x')
    %errorbar(Time, MEAN_SENSITIVE_500_BF(c_index,:), STD_SENSITIVE_500_BF(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Sensitive 500 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_BF_sensitive_500.png")

figure
colororder(newcolors)
for c_index=1:N_c
    plot(Time, MEAN_SENSITIVE_1000_BF(c_index,:), '-x')
    %errorbar(Time, MEAN_SENSITIVE_1000_BF(c_index,:), STD_SENSITIVE_1000_BF(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Sensitive 1000 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_BF_sensitive_1000.png")


figure
colororder(newcolors)
for c_index=1:N_c
    plot(Time, MEAN_RESISTANT_250_BF(c_index,:), '-x')
    %errorbar(Time, MEAN_RESISTANT_250_BF(c_index,:), STD_RESISTANT_250_BF(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Resistant 250 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_BF_resistant_250.png")


figure
colororder(newcolors)
for c_index=1:N_c
    plot(Time_res500bf, MEAN_RESISTANT_500_BF(c_index,:), '-x')
    %errorbar(Time_res500bf, MEAN_RESISTANT_500_BF(c_index,:), STD_RESISTANT_500_BF(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Resistant 500 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_BF_resistant_500.png")


figure
colororder(newcolors)
for c_index=1:N_c
    plot(Time, MEAN_SENSITIVE_500_GFP(c_index,:), '-x')
    %errorbar(Time, MEAN_SENSITIVE_500_GFP(c_index,:), STD_SENSITIVE_500_GFP(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Sensitive 500 GFP (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_GFP_sensitive_500.png")


figure
colororder(newcolors)
for c_index=1:N_c
    plot(Time, MEAN_SENSITIVE_1000_GFP(c_index,:), '-x')
    %errorbar(Time, MEAN_SENSITIVE_1000_GFP(c_index,:), STD_SENSITIVE_1000_GFP(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Sensitive 1000 GFP (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_GFP_sensitive_1000.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% with errorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tlabels = strings(1,N_t);
for t_index=1:N_t
    tlabels(N_t-t_index+1) = strcat("T=", int2str(Time(t_index)));
end

% Time x axis
figure
colororder(newcolors)
for c_index=1:N_c
    errorbar(Time, MEAN_SENSITIVE_500_BF(c_index,:), STD_SENSITIVE_500_BF(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Sensitive 500 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_errorbar_BF_sensitive_500.png")


figure
colororder(newcolors)
for c_index=1:N_c
    errorbar(Time, MEAN_SENSITIVE_1000_BF(c_index,:), STD_SENSITIVE_1000_BF(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Sensitive 1000 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_errorbar_BF_sensitive_1000.png")


figure
colororder(newcolors)
for c_index=1:N_c
    errorbar(Time, MEAN_RESISTANT_250_BF(c_index,:), STD_RESISTANT_250_BF(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Resistant 250 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_errorbar_BF_resistant_250.png")


figure
colororder(newcolors)
for c_index=1:N_c
    errorbar(Time_res500bf, MEAN_RESISTANT_500_BF(c_index,:), STD_RESISTANT_500_BF(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Resistant 500 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_errorbar_BF_resistant_500.png")


figure
colororder(newcolors)
for c_index=1:N_c
    errorbar(Time, MEAN_SENSITIVE_500_GFP(c_index,:), STD_SENSITIVE_500_GFP(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Sensitive 500 GFP (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_errorbar_GFP_sensitive_500.png")


figure
colororder(newcolors)
for c_index=1:N_c
    errorbar(Time, MEAN_SENSITIVE_1000_GFP(c_index,:), STD_SENSITIVE_1000_GFP(c_index, :))
    if c_index == 1
        hold on
    end
end
ylabel('Cell count')
xlabel('Time (h)')
ylim([0 inf])
legend(conclabels, 'location', 'northwest')
title("Sensitive 1000 GFP (averaged)")
saveas(gcf, "plots/dagim_monoclonal/mono_errorbar_GFP_sensitive_1000.png")

end %(if PLOT_CLEAN_DATA loop)

function [a] = fit_with_fixed_point(x,y,N_t)
    % Function y = y1 + a(x-x1) fitted using least squares
    % This fixes initial measurement
    % Nan values are skipped
    x1 = x(1,2);
    y1 = y(1);
    if isnan(y1)
        a = NaN;
    else
        xdata = x(2:length(y),2);
        ydata = y(2:length(y));
        fun = @(a,xdata) y1 + a*(xdata-x1);
        initial_a = 2;
        %options = optimoptions(@lsqcurvefit, 'Display','off');
        x_mask = ~isnan(ydata);
        y_mask = ~isnan(ydata);
        a = lsqcurvefit(fun, initial_a, xdata(x_mask), ydata(y_mask)); %, options);
    end
    end
