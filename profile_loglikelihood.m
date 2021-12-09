% Plot profile loglikelihood at selected solution

% Parameters
slurm_array_task_id = 330
savebool = true;
MixParam = 0.4;
resistant_dagim_5 = [0.035, 0.9, 1.3, 9.1344];
sensitive_dagim = [0.035, 0.9, 0.3, 9.1344];
resistant_cell_line = resistant_dagim_5;
sensitive_cell_line = sensitive_dagim;
Noise = 100;
R = 3;
N_c = 9;
N_t = 9;

% DATA generation    
outer_MixParams = 0.1:0.1:0.9;
outer_Noises = [0, 50, 100, 200, 500, 1000];
N_MixParams = length(outer_MixParams);
N_Noises = length(outer_Noises);
N_deltas = 1; %length(outer_deltas);
N_num_optims = 1; %length(outer_num_optims);

gridlength = 4;
N_c_values = 2.^(1:gridlength)+1; % [3,5,9,17];
N_t_values = 2.^(1:gridlength)+1; % [3,5,9,17];

total_permutations = N_MixParams * N_Noises * N_deltas * N_num_optims * gridlength * gridlength;

% MixParam
arr_size_MixParam = total_permutations / N_MixParams;
index_MixParam = fix(slurm_array_task_id / arr_size_MixParam) + 1;
MixParam = outer_MixParams(index_MixParam)
rest_MixParam = mod(slurm_array_task_id, arr_size_MixParam);

% Noise
arr_size_Noise = arr_size_MixParam / N_Noises;
index_Noise = fix(rest_MixParam / arr_size_Noise) + 1;
Noise = outer_Noises(index_Noise)
rest_Noise = mod(rest_MixParam, arr_size_Noise);

% N_c and N_t
n_c_index = fix(rest_Noise/gridlength) + 1; % whole integer division +1 is the row
n_t_index = mod(rest_Noise, gridlength) + 1; % The rest +1 is the column   
N_t = N_t_values(n_t_index)
N_c = N_c_values(n_c_index)

%Noise = 100
delta_E = 'med2'
num_optim = 200
N_repetitions = 1
seednr_1 = 43
effective_MixParam = min(MixParam, 1-MixParam);
Conc_num = 15

R_values = [3 6 12 15 24 30 50];
max_R = R_values(end);

% Generate data with noise
max_N_c = max(N_c_values);
max_N_t = max(N_t_values);

% Choose concentration locations
% Setup 1
Conc1 = zeros(1,max_N_c);
Conc1(2:max_N_c) = logspace(-1,1,(max_N_c-1)); % Min C = 0, then logspace from 10-1 to 10^1
% Setup 2
Conc2 = zeros(1,max_N_c);
Conc2(2:max_N_c) = logspace(-1,1,(max_N_c-1));
Conc2(2:9) = [(Conc1(3)+Conc1(4))/2, Conc1(4), (Conc1(4)+Conc1(5))/2, Conc1(5), (Conc1(5)+Conc1(6))/2, Conc1(6), (Conc1(6)+Conc1(7))/2, Conc1(7)];
Conc2(10:17) = [(Conc1(8)+Conc1(9))/2, Conc1(9), (Conc1(9)+Conc1(10))/2, Conc1(10), (Conc1(10)+Conc1(11))/2, Conc1(11), (Conc1(11)+Conc1(12))/2, Conc1(12)];
Nested_conc = Conc1;

resistant_cell_line = resistant_dagim_5;
sensitive_cell_line = sensitive_dagim;
sensitive_alpha = sensitive_cell_line(1);
resistant_alpha = resistant_cell_line(1);

Nested_time = linspace(0,96,max_N_t); % Min T = 0, max T = 96

NOISELESS_NESTED_DATA = generateDataEven(Nested_conc,Nested_time,max_R,MixParam,sensitive_cell_line,resistant_cell_line);
NESTED_DATA = zeros(max_R,max_N_c,max_N_t, N_repetitions);
for i_rep=1:N_repetitions
    rng(seednr_1 + i_rep - 1); % set seed
    NESTED_DATA(:,:,:,i_rep) = max(0, NOISELESS_NESTED_DATA + Noise*randn(max_R,max_N_c,max_N_t));
end

% The first 3 replicates are repeated over and over. 
NESTED_DATA(4:6,:,:,i_rep) = NESTED_DATA(1:3,:,:,i_rep);
NESTED_DATA(7:12,:,:,i_rep) = NESTED_DATA(1:6,:,:,i_rep);
NESTED_DATA(13:24,:,:,i_rep) = NESTED_DATA(1:12,:,:,i_rep);
NESTED_DATA(25:48,:,:,i_rep) = NESTED_DATA(1:24,:,:,i_rep);
NESTED_DATA(49:50,:,:,i_rep) = NESTED_DATA(1:2,:,:,i_rep);

% Define Time and Conc for this N_t and N_c
stepsize_t = 2^gridlength / (N_t - 1);
stepsize_c = 2^gridlength / (N_c - 1);
Conc = Nested_conc(1:stepsize_c:max_N_c)
Time = Nested_time(1:stepsize_t:max_N_t)

i_rep = 1;
DATA = NESTED_DATA(1:R, 1:stepsize_c:max_N_c, 1:stepsize_t:max_N_t, i_rep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_true = [MixParam sensitive_cell_line resistant_cell_line Noise]';
% Initial position for n = 6
x_final = [0.056744235947 0.097448309444 0.728734633501 1.757987561658 7.076051382591 0.079960460166 0.645561854508 2.073583436334 7.060310112730 1233.247662223828]';
% Some final position that we like
%x_final = [0.381105797002324, 0.036552392054100, 0.916135370903437, 0.275666078927534 sensitive_cell_line(4) resistant_cell_line Noise]';
% Some final position that is zero
%x_final = [0.00000, 0.053078279076776, 0.100621709151325, 4.980656679062155 sensitive_cell_line(4) resistant_cell_line Noise]';
% (True and final cannot be the same)

%f=@(x)vectorized_objective_function([x' sensitive_cell_line(4) resistant_cell_line Noise]',DATA,Conc,Time,R);
f=@(x)vectorized_objective_function(x,DATA,Conc,Time,R);

"f value at truth (which is not a local minimum):"
f(x_true)
%"f value at reasonable estimate"
%f(x_good)
"f value at estimate with p=0"
f(x_final)

% Plot likelihood profile for all parameters
paramtext = ["p", "alpha1", "b1", "E1", "n1", "alpha2", "b2", "E2", "n2", "sig"];
for iii=1:10;
    NNN = 100;
    inter_dist = abs(x_true(iii) - x_final(iii));
    test_values = min(x_true(iii), x_final(iii))-2*inter_dist:inter_dist/20:max(x_true(iii),x_final(iii))+2*inter_dist;
    % 1) Given true cell lines
    %f_values = zeros(NNN,1);
    %x_test = x_true;
    %for ii=1:length(test_values)
    %    x_test(iii) = test_values(ii);
    %    f_values(ii) = f(x_test);
    %end
    %figure
    %hold on
    %title(strcat("Given true cell lines, Mix ", num2str(MixParam), ", Noise ", int2str(Noise), ", Nc ", int2str(N_c), ", Nt ", int2str(N_t), ", R ", int2str(R)))
    %xline(x_true(iii), "-k")
    %xline(x_final(iii), "-r")
    %xlabel(paramtext(iii))
    %plot(test_values, f_values, '-ko')
    %legend("True", "ML Estimate", "loglikelihood")
    %hold off
    %if savebool
    %    saveas(gcf,[pwd strcat('/plots/profile-llhd/vec-given_true_cell_lines-seed-', int2str(seednr_1), '-num_optim-', int2str(num_optim), '-delta-E-', delta_E, '-MixParam-', num2str(MixParam), '-Noise-', int2str(Noise), '-Nc-', int2str(N_c), '-Nt_', int2str(N_t), '-R-', int2str(R), '-', int2str(iii), '.png')])
    %end

    % At estimated values
    f_values = zeros(NNN,1);
    x_test = x_final;
    for ii=1:length(test_values)
        x_test(iii) = test_values(ii);
        f_values(ii) = f(x_test);
    end
    if iii < 6
        figure('Position',[10 + (iii-1)*500 510 500 400])
    else 
        figure('Position',[10 + (iii-6)*500 10 500 400])
    end
    hold on
    title(strcat("With estimated cell lines, Mix ", num2str(MixParam), ", Noise ", int2str(Noise), ", Nc ", int2str(N_c), ", Nt ", int2str(N_t), ", R ", int2str(R)))
    xline(x_true(iii), "-k")
    xline(x_final(iii), "-r")
    xlabel(paramtext(iii))
    plot(test_values, f_values, '-ko')
    legend("Truth", "Iteration", "loglikelihood")
    hold off
    if savebool
        saveas(gcf,[pwd strcat('/plots/profile-llhd/vec-estimated_cell_lines-seed-', int2str(seednr_1), '-num_optim-', int2str(num_optim), '-delta-E-', delta_E, '-MixParam-', num2str(MixParam), '-Noise-', int2str(Noise), '-Nc-', int2str(N_c), '-Nt_', int2str(N_t), '-R-', int2str(R), '-', int2str(iii), '.png')])
    end
end

[x_true(1:9), x_final(1:9)]
