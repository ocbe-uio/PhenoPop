load([pwd '/data/elbow_just_loglikelihood-bic_values_ALL-46-gen-6-N_true_pop-3-num_optim-200-N_settings-16-Noise-300-Elim-1-plim-0-seed-46.mat'])

% ?
 %26.05

% Plot all test loglikelihood values
x = linspace(1, max_no_populations, max_no_populations);
fig = figure;
h = axes;
%set(h,'yscale','log')
xlabel('Number of inferred populations')
ylabel('Average test loglikelihood')
title(['Truth: ' num2str(N_true_populations) ' population(s)'])
hold on
for i_rep = 1:N_settings
    plot(x, average_loglikelihoods(i_rep, :), '.-k')
end
for i_rep = 1:N_settings
    differences_in_averages = average_loglikelihoods(i_rep, :) - min(average_loglikelihoods(i_rep, :));
    [min_value, min_pos] = min(differences_in_averages);
    plot(min_pos, average_loglikelihoods(i_rep, min_pos), '.', 'MarkerSize', 30)

    %standardized_differences = differences_in_averages / max(differences_in_averages);
    %plot(x, differences_in_averages, '.-', 'LineWidth', 1, 'MarkerSize', 10)
    %plot(x, standardized_differences, '.-', 'LineWidth', 1, 'MarkerSize', 10)
end
%saveas(gcf, [pwd, '/plots/crossvalidation/elbow-cv-', num2str(generation_number), '-Ntrue-', num2str(N_true_populations), '-N_settings-', num2str(N_settings), '-num_optim-', num2str(num_optim), '-', savestring, '.png'])
