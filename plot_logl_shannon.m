case_array = ["TREAT E11", "TREAT E14", "TREAT E19", "TREAT E1975", "TREAT E41", "TREAT E827", "TREAT P11", "TREAT P14", "TREAT P19", "TREAT P1975", "TREAT P41", "TREAT P827"];
savestr_array = ["TREAT_E_11", "TREAT_E_14", "TREAT_E_19", "TREAT_E_1975", "TREAT_E_41", "TREAT_E_827", "TREAT_P_11", "TREAT_P_14", "TREAT_P_19", "TREAT_P_1975", "TREAT_P_41", "TREAT P827"];
max_no_populations = 5 % The highest number of populations to try fitting to the data.

%%%%%%%%% Find max values %%%%%%%%%%
max_neg_ll = 0;
min_neg_ll = Inf;
max_difference = 0;

for casenumber = 1:length(case_array) %%%%%%%%%%%%%%%% parfor loop over drugs
casename = case_array(casenumber);
savestr1 = savestr_array(casenumber);
savestr = strcat(savestr1, '-num_optim-', num2str(num_optim));

load([strcat('./data/shannon_data/negative_loglikelihood_values-', savestr, '.mat')])
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

%%%%%%%%% Plot %%%%%%%%%%
for casenumber = 1:length(case_array) %%%%%%%%%%%%%%%% parfor loop over drugs
casename = case_array(casenumber);
savestr1 = savestr_array(casenumber);
savestr = strcat(savestr1, '-num_optim-', num2str(num_optim));

load([strcat('./data/shannon_data/negative_loglikelihood_values-', savestr, '.mat')])

% Plot the negative loglikelihood values
fig = figure;
hold on
xlabel('Number of inferred populations')
ylabel('Negative loglikelihood')
%ylim([min_neg_ll max_neg_ll])
ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
title([strcat(newline, casename) strcat('Negative loglikelihood')])
plot(1:max_no_populations, negative_loglikelihood_values, '.-k')
saveas(gcf, [pwd, '/plots/shannon_data/negative_loglikelihood_values_scaled_', num2str(savestr), '.png'])
end 
