set(0,'DefaultFigureVisible','on')
% Plotting different hill functions to see the effect of delta_E
newcolors0 = [
    0 0 0.726562
    0 0 0.726562
    0.398437 0.398437 0.9765625
    0.398437 0.398437 0.9765625
    0.792968 0 0
    0.792968 0 0
];

%x = zeros(1,1000);
%x(2:1000) = logspace(-1,1,999);
max_N_c = 17;

% Setup 1
Conc1 = zeros(1,max_N_c);
Conc1(2:max_N_c) = logspace(-1,1,(max_N_c-1)); % Min C = 0, then logspace from 10-1 to 10^1

% Setup 2
Conc2 = zeros(1,max_N_c);
Conc2(2:max_N_c) = logspace(-1,1,(max_N_c-1));
Conc2(2:9) = [(Conc1(3)+Conc1(4))/2, Conc1(4), (Conc1(4)+Conc1(5))/2, Conc1(5), (Conc1(5)+Conc1(6))/2, Conc1(6), (Conc1(6)+Conc1(7))/2, Conc1(7)];
Conc2(10:17) = [(Conc1(8)+Conc1(9))/2, Conc1(9), (Conc1(9)+Conc1(10))/2, Conc1(10), (Conc1(10)+Conc1(11))/2, Conc1(11), (Conc1(11)+Conc1(12))/2, Conc1(12)];

%load MIN_BETA.mat
%load CELL_LINE_PARAMETERS.mat %resistant_cell_line high_dE medium_dE low_dE

% Plot hill functions
%y_hi = hillfunc(high_dE, x);
%y_med = hillfunc(medium_dE, x);
%y_med2 = hillfunc(medium_2_dE, x);
%y_lo = hillfunc(low_dE, x);
%y_resistant_cell_line = hillfunc(resistant_cell_line, x);
%y_827 = hillfunc(min_beta_E_827, x);
%y_1975 = hillfunc(min_beta_E_1975, x);
%figure
%h = axes;
%set(h,'xscale','log')
%ylim([0.98 1.005])
%for i=2:length(Conc)
%    line([Conc(i), Conc(i)], [0.98 1.005], 'Color', 'k');
%    hold on 
%end
%semilogx(x,y_resistant_cell_line, 'LineWidth', 2.5)
%%semilogx(x,y_med, 'LineWidth', 2.5)
%semilogx(x,y_med2, 'LineWidth', 2.5)
%ylabel('Hill function value')
%xlabel('Dose (x)')
%title("Cell lines, medium2 dE")
%legend('resistant cell line', 'sensitive cell line') % "medium dE", "medium 2 dE", "high dE") %, 'Low dE', 'Medium dE', 'High dE', 'H1975') %'High dE', 

%resistant_cell_line = min_beta_E_1975;
%sensitive_cell_line = min_beta_E_827;
%resistant_cell_line = [0.0275 0.9877 1.3144  9.1344];
%sensitive_cell_line = [0.0275 0.9877  0.3  9.1344];

% Old setup with med2: 
% Res:  [0.0275 0.9877 1.3144  9.1344]
% Sens: [0.0275 0.9877  0.3  9.1344]

% New ones that worked
% resistant_cell_line = [0.035, 0.9, 1.3, 9.1344];
% sensitive_cell_line = [0.035, 0.9, 0.3, 9.1344];

%resistant_cell_line = min_beta_E_1975;
%sensitive_cell_line = min_beta_E_827;

% Dagim's
%resistant_cell_line = [0.035, 0.9, 1.3, 9.1344];
%sensitive_cell_line = [0.035, 0.9, 0.3, 9.1344];

%resistant_cell_line = [base_rate_resistant, b, middle_E + dE/2, n] %2*(middle_E + dE/2), n]
%sensitive_cell_line = [base_rate_sensitive, b, middle_E - dE/2, n] %2*(middle_E - dE/2), n]

generation_number = 6
N_c = 17
stepsize_c = 2^4 / (N_c - 1);

% Plot growth rates
base_rate_resistant = 0.03;
base_rate_sensitive = 0.03;
b = 0.9;
middle_E = 0.857695898590894;
dE = 1;
n = 9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if generation_number == 5
    % Generation 5, Dagim's data
    base_rate_resistant = 0.03;
    base_rate_sensitive = 0.03;
    b_values = [0.3 0.5 0.7 0.9];
    outer_n_values = [1 3 5 7 9];
    outer_MixParams = 0.1:0.1:0.9;
    Noise_values = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500];
    E_sensitive = 1.0e-3;
    E_factors = [4 16 64 256];
elseif generation_number == 6
    % June 2021. Based on Dagim's data. Concentrations set loglinearly between lowest and highest
    base_rate_resistant = 0.03;
    base_rate_sensitive = 0.03;
    b_values = [0.3 0.5 0.7 0.9];
    outer_n_values = [0.1 0.3 0.5 0.7 0.9];
    outer_MixParams = 0.1:0.1:0.9;
    Noise_values = [0, 50, 100, 150, 200];
    E_sensitive = 1.0e-3;
    E_factors = [4 16 64 256];
    % Conc is set further down
end

figure
h = axes;
set(h,'xscale','log')
for slurm_array_task_id = 0:15 % MIX PARAM DOESN'T MATTER SO ONLY LOOK AT FIRST 16*5=80
    %if slurm_array_task_id == 3 || slurm_array_task_id == 6 || slurm_array_task_id == 9 || slurm_array_task_id == 12
    %if slurm_array_task_id == 3 || slurm_array_task_id == 7 || slurm_array_task_id == 11 || slurm_array_task_id == 15
    if slurm_array_task_id == 3 || slurm_array_task_id == 2 || slurm_array_task_id == 1 || slurm_array_task_id == 0
    slurm_array_task_id = slurm_array_task_id + 0*16 % Add multiples of 16 to move to next n value
    gridlength_b = length(b_values);
    %gridlength_E = length(dE_values);
    gridlength_E = length(E_factors);
    N_MixParams = length(outer_MixParams);
    N_n = length(outer_n_values);
    max_Noise = Noise_values(end);
    total_combinations = N_MixParams * N_n * gridlength_E * gridlength_b;

    % Find MixParam
    arr_size_MixParam = total_combinations / N_MixParams;
    index_MixParam = fix(slurm_array_task_id / arr_size_MixParam) + 1;
    MixParam = outer_MixParams(index_MixParam)
    effective_MixParam = min(MixParam, 1-MixParam);
    rest_MixParam = mod(slurm_array_task_id, arr_size_MixParam);

    % Find n value
    arr_size_Noise = arr_size_MixParam / N_n;
    index_n = fix(rest_MixParam / arr_size_Noise) + 1;
    n = outer_n_values(index_n) 
    rest_n = mod(rest_MixParam, arr_size_Noise);

    % Height and width
    b_index = fix(rest_n/gridlength_E) + 1; % whole integer division +1 is the row
    dE_index = mod(rest_n, gridlength_E) + 1; % The rest +1 is the column   
    % Set height by adjusting the b parameter
    b = b_values(b_index)
    % Set width by adjusting the E parameter of the resistant cell line
    %dE = dE_values(dE_index)
    %E_sensitive = middle_E + dE/2;
    %E_resistant = middle_E - dE/2;
    E_factor = E_factors(dE_index)
    E_resistant = E_sensitive * E_factor;
    sensitive_cell_line = [base_rate_sensitive, b, E_sensitive, n]
    resistant_cell_line = [base_rate_resistant, b, E_resistant, n]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if generation_number > 5
        Nested_conc = zeros(1,max_N_c);
        % The concentrations are logarithmically spaced between these values. 
        % The factors 0.05 * sensitive E was chosen from visual inspection of data plots
        % The factor 5 * resistant E was chosen from visual inspection of logarithmic plots of the data
        Nested_conc(2:max_N_c) = logspace(log10(0.05*sensitive_cell_line(3)),log10(5*resistant_cell_line(3)),(max_N_c-1)); % Min C = 0, then logspace from 10-1 to 10^1
        Conc = Nested_conc(1:stepsize_c:max_N_c);
    else
        Conc3 = find_the_right_concentrations(max_N_c, sensitive_cell_line, resistant_cell_line);
        Conc = Conc3(1:stepsize_c:max_N_c);
    end

    x = zeros(1,1000);
    x(2:1000) = logspace(log10(0.1*min(Conc(2:N_c))),log10(10^2*max(Conc(2:N_c))),999);
    y_resistant = ratefunc(resistant_cell_line, x);
    y_sensitive = ratefunc(sensitive_cell_line, x);

    semilogx(x,y_sensitive, '-b', 'LineWidth', 2.5)
    hold on
    semilogx(x,y_resistant, '-r', 'LineWidth', 2.5)
    %line(xlim, [0.035 0.035], 'Color', 'k' )
    %line(xlim, [0.035+log(0.9) 0.035+log(0.9)], 'Color', 'k' )
    %line(xlim, [0 0], 'Color', 'k' )
    %for i=2:length(Conc)
    %    line([Conc(i), Conc(i)], ylim, 'Color', 'k');
    %end
    %rrr = Conc(17);
    %line([rrr rrr], ylim, 'Color', 'r');
    %axis([10^(-6.5) 10^-4.5 -inf inf])
    %ylim([-0.08, 0.05])
    end
end
ylabel('Growth rate')
xlabel('Dose (x)')  
xlim([10^-5.6,10^2.3])
title(['n=', num2str(sensitive_cell_line(4))])
%title(['Growth rates' newline '             alpha          b             E          n' newline 'Sensitive   ', num2str(sensitive_cell_line), newline 'Resistant   ', num2str(resistant_cell_line)])
%legend('Sensitive', 'Resistant')
saveas(gcf, strcat("article/all_hills_reduced_n_", num2str(sensitive_cell_line(4)), ".png"))

% The total rate should become negative for high enough dose.
% This happens if (alpha + log(b) < 0) 
%high_dE(1) + log(high_dE(2));

%semilogx(x,y_827, 'LineWidth', 2)
%plot(x,y_hi)
%plot(x,y_lo)
%plot(x,y_1975)
%ylim([0.984,1.001]) 
