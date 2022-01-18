patno_array = ["MM210819", "MM130120", "MM19520", "MM3620", "MM1420"];
max_no_populations = 5 % The highest number of populations to try fitting to the data.


%%%%%%%%% Find max values %%%%%%%%%%
max_neg_ll = 0;
min_neg_ll = Inf;
max_difference = 0;

for patient_index = 1:5
patient_number = patno_array(patient_index)
for drugnumber = 1:4 %%%%%%%%%%%%%%%% parfor loop over drugs

load([strcat('./data/MM-patient_sample_data/negative_loglikelihood_values-', num2str(patient_number), '-drug-', num2str(drugnumber), '.mat')])
if max(negative_loglikelihood_values) > max_neg_ll
    max_neg_ll = max(negative_loglikelihood_values);
end
if min(negative_loglikelihood_values) < min_neg_ll
    min_neg_ll = min(negative_loglikelihood_values);
end

if max(negative_loglikelihood_values) - min(negative_loglikelihood_values) > max_difference
    max_difference = max(negative_loglikelihood_values) - min(negative_loglikelihood_values);
end

end 
end 

%%%%%%%%% Plot %%%%%%%%%%
for patient_index = 1:5
patient_number = patno_array(patient_index)
for drugnumber = 1:4 %%%%%%%%%%%%%%%% parfor loop over drugs

%load([strcat('./data/MM-patient_sample_data/negative_loglikelihood_values-', num2str(patient_number), '-drug-', num2str(drugnumber), '.mat')])
num_optim = 3000;
savestr = strcat(num2str(patient_number), '-drug-', num2str(drugnumber), '-num_optim-', num2str(num_optim));
load([strcat('./plots/MM-patient_sample_data/negLL/negative_loglikelihood_and_GR50-', savestr, '.mat')])

% Plot the negative loglikelihood values
fig = figure;
hold on
xlabel('Number of inferred populations')
ylabel('Negative loglikelihood')
%ylim([min_neg_ll max_neg_ll])
ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
title(strcat(patient_number, ' drug', {' '}, num2str(drugnumber)))
plot(1:max_no_populations, negative_loglikelihood_values, '.-k')
saveas(gcf, [pwd, '/plots/MM-patient_sample_data/negLL/negative_loglikelihood_values-', num2str(patient_number), '-drug-', num2str(drugnumber), '.png'])

end 
end 
