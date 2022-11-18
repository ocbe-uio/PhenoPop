% Unpack and plot Dagim's monoclonal data
% 4 different settings: - Sensitive 500 and 1000 initial cells
%                       - Resistant 250 and 500 initial cells
% All 4 carried out with Bright field imaging (BF)
% Both sensitive also done with Green fluorescent protein imaging (GFP)
% ... for 7 replicates, 25 time points (0 to 72 hours with 3 hour intervals), and 11 concentrations.

PLOTTING = false

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
% 1) Take the logarithm to transform exponential data to linear data
% 2) Find the slope of each replicate for each concentration (rate in exponential growth) and calculate residuals w.r.t this line
% 3) Remove outliers using Matlab function "isoutlier" on the residuals

% All data:
%SENSITIVE_500_BF(:,:,:)
%SENSITIVE_1000_BF(:,:,:)
%RESISTANT_250_BF(:,:,:)
%RESISTANT_500_BF(:,:,:)
%SENSITIVE_500_GFP(:,:,:)
%SENSITIVE_1000_GFP(:,:,:)

% 1) Take the logarithm to transform exponential data to linear data
logSENSITIVE_500_BF = log(SENSITIVE_500_BF);
logSENSITIVE_1000_BF = log(SENSITIVE_1000_BF);
logRESISTANT_250_BF = log(RESISTANT_250_BF);
logRESISTANT_500_BF = log(RESISTANT_500_BF);
logSENSITIVE_500_GFP = log(SENSITIVE_500_GFP);
logSENSITIVE_1000_GFP = log(SENSITIVE_1000_GFP);

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
%saveas(gcf, "plots/dagim_monoclonal/mono_GFP_sensitive_500.png")

% For each concentration and timepoint there are 7 replicates. Among these we remove the outlier. 
%for time_index=1:3 %N_t
%    SENSITIVE_500_GFP(:,:,time_index)
%    TF = isoutlier(SENSITIVE_500_GFP(:,:,time_index), 1)
%    SENSITIVE_500_GFP(:,:,time_index)(TF) = NaN;
%end

%BOOL_ARRAY = zeros(R, N_c, N_t);
%% For each replicate evolving in time we want to remove measurement erros. 
%for rep_index=1:R
%    for conc_index=1:8 %N_c
%        TF = squeeze(isoutlier(logSENSITIVE_500_GFP(rep_index,conc_index,:), 'movmedian', 9, 'ThresholdFactor', 3))';
%        % We leave the first and last measurement as they are since they are harder to detect as outliers:
%        %TF(1) = 0;
%        %TF(25) = 0;
%        if sum(TF) > 0
%            rep_index
%            conc_index
%            squeeze(SENSITIVE_500_GFP(rep_index,conc_index,:))'
%            TF
%        end
%        BOOL_ARRAY(rep_index,conc_index,:) = TF;
%    %    SENSITIVE_500_GFP(:,:,time_index)(TF) = NaN;
%    end
%end
%%SENSITIVE_500_GFP(BOOL_ARRAY) = NaN;

% 2) Find the slope of each replicate for each concentration (rate in exponential growth) and calculate residuals w.r.t this line
% 3) Remove outliers using Matlab function "isoutlier" on the residuals
x = ones(N_t,2);
x(:,2) = 3*(1:N_t)'

BOOL_ARRAY = zeros(R, N_c, N_t);
% For each replicate evolving in time we want to remove measurement erros. 
for rep_index=1:1
    for conc_index=1:N_c
        %squeeze(SENSITIVE_500_BF(rep_index,conc_index,:))
        squeeze(MEAN_RESISTANT_250_BF(conc_index,:))
        % Fit a linear model to the data
        %y = squeeze(logSENSITIVE_500_BF(rep_index,conc_index,:));
        y = squeeze(logMEAN_RESISTANT_250_BF(conc_index,:))'
        coef_beta = x\y;
        ypred = x*coef_beta;

        figure
        scatter(x(:,2),y)
        hold on
        plot(x(:,2), ypred, 'r')
        ylim([0 inf])
        grid on
        hold off
        saveas(gcf,[pwd strcat('/plots/dagim_monoclonal/logMEAN_concentration_', int2str(conc_index), '.png')])

%        return
%
%        y_residuals = y - ypred;
%        TF = isoutlier(y_residuals, 'grubbs'); %, 'movmedian', 9, 'ThresholdFactor', 3)
%
%        return
%
%        TF = squeeze(isoutlier(logSENSITIVE_500_BF(rep_index,conc_index,:), 'movmedian', 9, 'ThresholdFactor', 3))';
%        % We leave the first and last measurement as they are since they are harder to detect as outliers:
%        %TF(1) = 0;
%        %TF(25) = 0;
%        if sum(TF) > 0
%            rep_index
%            conc_index
%            squeeze(SENSITIVE_1000_GFP(rep_index,conc_index,:))'
%            TF
%        end
%        BOOL_ARRAY(rep_index,conc_index,:) = TF;
%    %    SENSITIVE_500_GFP(:,:,time_index)(TF) = NaN;
    end
end

return

figure
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


return

% Observed oddities for several replicates: 
SENSITIVE_1000_GFP(2:7,11,14);

% All replicates of zero concentration for Resistant 250 BF are collectively wrong for times 8 to 14.
RESISTANT_250_BF(:,1,8:14);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data cleaning by outlier removal and removing cells that seem to rise from the dead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SENSITIVE_500_BF(:,:,:)
%SENSITIVE_1000_BF(:,:,:)
%RESISTANT_250_BF(:,:,:)
%RESISTANT_500_BF(:,:,:)
%SENSITIVE_500_GFP(:,:,:)
%SENSITIVE_1000_GFP(:,:,:)

%% sensitive 500 GFP
%squeeze(SENSITIVE_500_GFP(1,1,:))
SENSITIVE_500_GFP(1,1,5:7) = NaN; % check, plus some extra
%squeeze(SENSITIVE_500_GFP(1,1,:))

%% sensitive 1000 GFP
% Single observations: 
SENSITIVE_1000_GFP(7,3:4,6) = NaN; % check
SENSITIVE_1000_GFP(7,5,7) = NaN; % check
SENSITIVE_1000_GFP(7,8:9,13) = NaN; % check for 8, not 9
SENSITIVE_1000_GFP(2:7,11,14) = NaN;

SENSITIVE_1000_GFP(7,10,15:17) = NaN; 
SENSITIVE_1000_GFP(7,11,15:16) = NaN;
SENSITIVE_1000_GFP(7,9,16) = NaN;
SENSITIVE_1000_GFP(7,9,18:20) = NaN;
SENSITIVE_1000_GFP(7,9,23) = NaN;
SENSITIVE_1000_GFP(7,11,23) = NaN;

% Replicate 7, conc 6 and 7, times 8 to 14, but we remove them from 8 and out.
%SENSITIVE_1000_GFP(:,6:7,8:25)
SENSITIVE_1000_GFP(7,6:7,8:25) = NaN;

% All replicates of zero concentration for Resistant 250 BF are collectively wrong for times 8 to 14.
% For inference we must censor entirely the zero concentration line.
RESISTANT_250_BF(:,1,8:14) = NaN;

%SENSITIVE_1000_GFP(:,:,:)
%RESISTANT_250_BF(:,:,:)

% Censor cells that seem to rise from the dead:
SENSITIVE_500_GFP(:,10:11,11:25) = 5;
SENSITIVE_500_GFP(:,9,14:25) = 5;
SENSITIVE_1000_GFP(:,10:11,12:25) = 1;

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

%SENSITIVE_500_BF(:,:,1)
%SENSITIVE_1000_BF(:,:,1)
%RESISTANT_250_BF(:,:,1)
%RESISTANT_500_BF(:,:,1)
%SENSITIVE_500_GFP(:,:,1)
%SENSITIVE_1000_GFP(:,:,1)
%RESISTANT_500_BF(2,3,14) %=3863

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Infer and plot estimated cell lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Plot hill functions
x=linspace(0,5,1000);
hill_sens_500_bf = hillfunc(x_final_sensitive_500_bf(2:5),x);
hill_sens_1000_bf = hillfunc(x_final_sensitive_1000_bf(2:5),x);
hill_res_250_bf = hillfunc(x_final_resistant_250_bf(2:5),x);
hill_res_500_bf = hillfunc(x_final_resistant_500_bf(2:5),x);
hill_sens_500_gfp = hillfunc(x_final_sensitive_500_gfp(2:5),x);
hill_sens_1000_gfp = hillfunc(x_final_sensitive_1000_gfp(2:5),x);

newcolors0 = [
    0 0 0.726562
    0 0 0.726562
    0.398437 0.398437 0.9765625
    0.398437 0.398437 0.9765625
    0.792968 0 0
    0.792968 0 0
];
figure
colororder(newcolors0)
hold on 
plot(x,hill_sens_500_bf, '--', 'linewidth', 2) %, '-o')
plot(x,hill_sens_1000_bf, '-', 'linewidth', 2) %, '-+')
plot(x,hill_sens_500_gfp, '--', 'linewidth', 2) %, '-^')
plot(x,hill_sens_1000_gfp, '-', 'linewidth', 2) %, '-p')
plot(x,hill_res_250_bf, '--', 'linewidth', 2) %, '->')
plot(x,hill_res_500_bf, '-', 'linewidth', 2) %, '-v')
ylabel('Hill function value')
xlabel('Concentration (x)')
title("Inferred hill functions")
legend('Sensitive 500 BF', 'Sensitive 1000 BF', 'Sensitive 500 GFP', 'Sensitive 1000 GFP', 'Resistant 250 BF', 'Resistant 500 BF')
saveas(gcf, "plots/dagim_monoclonal/inferred_cell_lines.png")

if PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting in time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tlabels = strings(1,N_t);
for t_index=1:N_t
    tlabels(N_t-t_index+1) = strcat("T=", int2str(Time(t_index)));
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
for t_index=N_t:-4:1
    %semilogx(Conc, MEAN_SENSITIVE_500_BF(:,t_index)', '-x')
    hold on
    errorbar(Conc, MEAN_SENSITIVE_500_BF(:,t_index)', STD_SENSITIVE_500_BF(:,t_index))
end
ylabel('Cell count')
xlabel('Concentration (M)')
ylim([0 inf])
legend(tlabels(1:4:N_t), 'location', 'northeast')
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
legend(tlabels(1:4:N_t), 'location', 'northeast')
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
legend(tlabels(1:4:N_t), 'location', 'northeast')
title("Resistant 250 BF (averaged)")
saveas(gcf, "plots/dagim_monoclonal/conc_axis_mono_BF_resistant_250.png")

figure
colororder(newcolors2)
hAx=axes;
hAx.XScale='log'
for t_index=N_t:-4:1
    %semilogx(Conc, MEAN_RESISTANT_500_BF(:,t_index)', '-x')
    hold on
    errorbar(Conc, MEAN_RESISTANT_500_BF(:,t_index)', STD_RESISTANT_500_BF(:,t_index))
end
ylabel('Cell count')
xlabel('Concentration (M)')
ylim([0 inf])
legend(tlabels(1:4:N_t), 'location', 'northeast')
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
legend(tlabels(1:4:N_t), 'location', 'northeast')
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
legend(tlabels(1:4:N_t), 'location', 'northeast')
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
    plot(Time, MEAN_RESISTANT_500_BF(c_index,:), '-x')
    %errorbar(Time, MEAN_RESISTANT_500_BF(c_index,:), STD_RESISTANT_500_BF(c_index, :))
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
    errorbar(Time, MEAN_RESISTANT_500_BF(c_index,:), STD_RESISTANT_500_BF(c_index, :))
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

end %(if PLOTTING loop)

