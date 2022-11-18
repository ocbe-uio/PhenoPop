
% load CELL_LINE_PARAMETERS.mat %fixed, high_dE meduim_dE low_dE
% delta_E = 'med';
% switch delta_E
% case 'hi'
%     sensitive_cell_line = high_dE % delta_E = 0.6
% case 'med'
%     sensitive_cell_line = medium_dE  % delta_E = 0.3
% case 'med2'
%     sensitive_cell_line = medium_2_dE  % delta_E = 0.3
% case 'low'
%     sensitive_cell_line = low_dE  % delta_E = 0.1
% end
% resistant_cell_line = resistant_dagim;
% sensitive_cell_line = sensitive_dagim;

% 'hi'
% sensitive_cell_line = [0.0275 0.983  0.01  9.1344]; % delta_E = 1.29
% 'med'
% sensitive_cell_line = [0.0275 0.983  0.7  9.1344];  % delta_E = 0.6
% low'
% sensitive_cell_line = [0.0275 0.983  1.2  9.1344];  % delta_E = 0.1

slurm_array_task_id = 262
mid_factor = 16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation 5, Dagim's data
base_rate_resistant = 0.03;
base_rate_sensitive = 0.03;
b_values = [0.3 0.5 0.7 0.9];
outer_n_values = [1 3 5 7 9];
outer_MixParams = 0.1:0.1:0.9;
Noise_values = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500];
E_sensitive = 1.0e-3;
E_factors = [4 16 64 256];

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
resistant_cell_line = [base_rate_resistant, 1, E_resistant, n]
middle_cell_line = [base_rate_resistant, b, mid_factor*E_sensitive, n] %    Conc(9) = Conc(5) * sqrt(Conc(12) / Conc(5));
MixParams = [0.3 0.3]; % Length k-1
MixParams = [0.5 0.5]; % Length k-1
MixParams = [1 0]; % Length k-1
MixParams = [0.5 0.5]; % Length k-1
MixParams = [0.5 0.0]; % Length k-1
Params = [sensitive_cell_line resistant_cell_line middle_cell_line];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_c = 17;
Time = [0 48 72]; %[0 12 24 36 48 60 72 84 96];
N_t = length(Time);
R = 1;
tlabels = strings(1,N_t);
for t_index=1:N_t
    tlabels(N_t-t_index+1) = strcat("T=", int2str(Time(t_index)));
end
stepsize_c = 2^4 / (N_c - 1);
%Conc3 = old_find_the_right_concentrations(N_c, sensitive_cell_line, resistant_cell_line);
Conc3 = find_the_right_concentrations(N_c, sensitive_cell_line, resistant_cell_line,Time);

max_N_c = N_c;
Conc4 = zeros(1,max_N_c);
Conc4(2:max_N_c) = logspace(log10(min(Conc3(2:max_N_c))),log10(max(Conc3(2:max_N_c))),(max_N_c-1)); % Min C = 0, then logspace from 10-1 to 10^1

Conc5 = zeros(1,max_N_c);
Conc5(2:max_N_c) = logspace(log10(0.05*E_sensitive),log10(E_resistant),(max_N_c-1)); % Min C = 0, then logspace from 10-1 to 10^1

Conc = Conc5(1:stepsize_c:max_N_c)
%temp_Nested_conc = zeros(1,max_N_c);
%temp_Nested_conc(2:max_N_c) = logspace(log10(min(Conc(2:N_c))),log10(max(Conc(2:N_c))),(max_N_c-1)); % Min C = 0, then logspace from 10-1 to 10^1
%Conc = temp_Nested_conc;

plot_N_c = 1000;
x = zeros(1,plot_N_c);
x(2:plot_N_c) = logspace(log10(0.1*min(Conc(2:N_c))),log10(10*max(Conc(2:N_c))),(plot_N_c-1));

%populations = generateDataEven(x,Time,R,MixParam,sensitive_cell_line,resistant_cell_line);
populations = generate_data_k_populations(x,Time,R,MixParams,Params);
max_cell_count = max(populations, [], "all");
% N_c = 17
% Conc = zeros(1,N_c)
% Conc(2:N_c) = logspace(-1,1,(N_c-1)) % Min C = 0, then logspace from 10-1 to 10^1

fig = figure %('Position',[800 800 400 300]);
movegui(fig,[1275 630]); % x y positions of bottom left corner
h = axes;
set(h,'xscale','log')
hold on
%for t_index=N_t:-1:1
%    semilogx(x, 100*populations(:,:,t_index)/max_cell_count, "linewidth",2)
%    hold on
%end
%for c_index=2:N_c
%    xline(Conc(c_index), "k")
%end
ylabel('Tissue response (% max)')
xlabel('Drug dose')
legend(tlabels)
%title(['Cell counts' newline '             alpha          b             E          n' newline 'Sensitive   ', num2str(sensitive_cell_line), newline 'Resistant   ', num2str(resistant_cell_line)])
saveas(gcf,[pwd strcat('/plots/difference/0presentation-cell_counts_at_some_times-slurmarray-', num2str(slurm_array_task_id), '.png')])

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

%% Plot data points and true viability curves
%% a) train data
%plot_N_c = 1000;
%increment = 1; % Step in time. Increase to plot less of the data points
%plot_Time = Time(1:increment:N_t); %[0 24 48 72 96]; %[0 12 24 26 48 60 72 84 96];
%plot_N_t = length(plot_Time);
%tlabels = strings(1,N_t);
%for t_index=1:plot_N_t
%    tlabels(plot_N_t-t_index+1) = strcat("T=", int2str(plot_Time(t_index)));
%end
%x = logspace(log10(0.1*min(Conc(2:N_c))),log10(10*max(Conc(2:N_c))),plot_N_c);
%plot_populations = generate_data_k_populations(x,plot_Time,1,MixParams,Params);
%
%fig = figure %('Position',[800 800 400 300]);
%colororder(repeated_colors);
%movegui(fig,[1275 630]); % x y positions of bottom left corner
%h = axes;
%%set(h,'xscale','log')
%for t_index=plot_N_t:-1:1
%    semilogx(x, plot_populations(:,:,t_index))
%    hold on
%    for r_index = 1:R
%        semilogx(Conc(2:N_c), DATA(r_index,2:N_c,increment*t_index - (increment-1)), '.', 'MarkerSize', 10,'HandleVisibility','off')
%    end
%end
%for c_index=2:N_c
%    xline(Conc(c_index), "k",'HandleVisibility','off')
%end
%ylabel('Cell count')
%xlabel('Drug dose')
%legend(tlabels)
%title(['Cell counts (Dots for 1 replicate)']) % newline '             alpha          b             E          n' newline 'Sensitive   ', num2str(sensitive_cell_line), newline 'Resistant   ', num2str(resistant_cell_line)])
%saveas(gcf, [pwd, '/plots/crossvalidation/logl-data-', num2str(generation_number), '-', num2str(i_rep), '-', savestring, '.png'])
