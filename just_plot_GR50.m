    %%%%%%%%%%%%%%%%%%%% Plot errors in sensitive GR50 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% PLOT TYPE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% L1 error in p; Noise on x axis
    h1 = figure;
    colororder(newcolors4)
    set(h1,'Position',[300 1100 1000 800])
    
    for ii=0:(length(E_factors)*length(b_values)-1)
        for ncntr_index=1:length(all_ncntr_values)
            b_index = fix(ii/length(E_factors)) + 1; % whole integer division +1 is the row
            E_index = mod(ii, length(E_factors)) + 1; % The rest +1 is the column
            b = b_values(1 + length(b_values) - b_index);
            E_factor = E_factors(E_index);
            averages_gr50_sens = mean(all_log_L1_gr50_sens(:,:,b_index, E_index, ncntr_index), 1);
            averages_gr50_res = mean(all_log_L1_gr50_res(:,:,b_index, E_index, ncntr_index), 1);
    
            ax1 = subplot(length(E_factors),length(b_values),ii+1);
            set(gca, 'YScale', 'log')
            hold on
    
            % Plot absolute error for mixture parameter
            %yyaxis left
            ylim([0 1.1*gr50_error_max])
            colororder(newcolors4)
            plot(Noise_values, averages_gr50_sens, '-o', 'LineWidth', 1.5, 'MarkerSize', 2)
            %if ncntr_index == 1
            %    plot(Noise_values, averages_p, '-ro') %, 'LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',4)
            %elseif ncntr_index == 2
            %    plot(Noise_values, averages_p, '-bo') %, 'LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',4)
            %elseif ncntr_index == 3
            %    plot(Noise_values, averages_p, '-go') %, 'LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',4)
            %elseif ncntr_index == 4
            %    plot(Noise_values, averages_p, '-mo') %, 'LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',4)
            %end
    
            xlabel('Noise')
            if E_index == 1
                ylabel('Abs. error (GR50)')
            end
            title(strcat('E factor=', int2str(E_factor), ', b=', num2str(b))) %, ', dE=', num2str(dE)))
            %hold off
        end % loop over different (Nc,Nt,R)
        %legend('N_c=N_t=R=3', 'N_c=N_t=R=5', 'N_c=N_t=R=9', 'N_c=N_t=R=17')
    end % loop over subplots
    % add global legend
    Lgnd = legend('N_c=N_t=R=3', 'N_c=N_t=R=5', 'N_c=N_t=R=9', 'N_c=N_t=R=17');
    Lgnd.Position(1) = 0.8; % fraction of x direction
    Lgnd.Position(2) = 0.9; % fraction of y direction
    sgtitle(['L1 error in sensitive GR50. n=', int2str(n), '. Mixture of ', num2str(MixParam), ' sensitive cells.'], 'FontSize', 12, 'HorizontalAlignment', 'left') %Nc=', int2str(N_c), ', Nt=', int2str(N_t), ', R=', int2str(R), '.'
    saveas(gcf,[pwd strcat('/plots/difference/noise_x_axis-gen-', num2str(generation_number), '-0-GR1-sens-Nc-', int2str(N_c), '-Nt-', int2str(N_t), '-R-', int2str(R), '-num_optim-', int2str(num_optim), '-N-rep-', int2str(N_repetitions), '-seed-', int2str(seednr_1), '-MixParam-', num2str(MixParam), '-n-', num2str(n), '.png')])

    %%%%%%%%%%%%%%%%%%%% Plot errors in resistant GR50 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% PLOT TYPE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% L1 error in p; Noise on x axis
    h1 = figure;
    colororder(newcolors4)
    set(h1,'Position',[300 1100 1000 800])
    
    for ii=0:(length(E_factors)*length(b_values)-1)
        for ncntr_index=1:length(all_ncntr_values)
            b_index = fix(ii/length(E_factors)) + 1; % whole integer division +1 is the row
            E_index = mod(ii, length(E_factors)) + 1; % The rest +1 is the column
            b = b_values(1 + length(b_values) - b_index);
            E_factor = E_factors(E_index);
            averages_gr50_sens = mean(all_log_L1_gr50_sens(:,:,b_index, E_index, ncntr_index), 1);
            averages_gr50_res = mean(all_log_L1_gr50_res(:,:,b_index, E_index, ncntr_index), 1);
    
            ax1 = subplot(length(E_factors),length(b_values),ii+1);
            set(gca, 'YScale', 'log')
            hold on
    
            % Plot absolute error for mixture parameter
            %yyaxis left
            ylim([0 1.1*gr50_error_max])
            colororder(newcolors4)
            plot(Noise_values, averages_gr50_res, '-o', 'LineWidth', 1.5, 'MarkerSize', 2)
            %if ncntr_index == 1
            %    plot(Noise_values, averages_p, '-ro') %, 'LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',4)
            %elseif ncntr_index == 2
            %    plot(Noise_values, averages_p, '-bo') %, 'LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',4)
            %elseif ncntr_index == 3
            %    plot(Noise_values, averages_p, '-go') %, 'LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',4)
            %elseif ncntr_index == 4
            %    plot(Noise_values, averages_p, '-mo') %, 'LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',4)
            %end
    
            xlabel('Noise')
            if E_index == 1
                ylabel('Abs. error (GR50)')
            end
            title(strcat('E factor=', int2str(E_factor), ', b=', num2str(b))) %, ', dE=', num2str(dE)))
            %hold off
        end % loop over different (Nc,Nt,R)
        %legend('N_c=N_t=R=3', 'N_c=N_t=R=5', 'N_c=N_t=R=9', 'N_c=N_t=R=17')
    end % loop over subplots
    % add global legend
    Lgnd = legend('N_c=N_t=R=3', 'N_c=N_t=R=5', 'N_c=N_t=R=9', 'N_c=N_t=R=17');
    Lgnd.Position(1) = 0.8; % fraction of x direction
    Lgnd.Position(2) = 0.9; % fraction of y direction
    sgtitle(['L1 error in resistant GR50. n=', int2str(n), '. Mixture of ', num2str(MixParam), ' sensitive cells.'], 'FontSize', 12, 'HorizontalAlignment', 'left') %Nc=', int2str(N_c), ', Nt=', int2str(N_t), ', R=', int2str(R), '.'
    saveas(gcf,[pwd strcat('/plots/difference/noise_x_axis-gen-', num2str(generation_number), '-0-GR2-res-Nc-', int2str(N_c), '-Nt-', int2str(N_t), '-R-', int2str(R), '-num_optim-', int2str(num_optim), '-N-rep-', int2str(N_repetitions), '-seed-', int2str(seednr_1), '-MixParam-', num2str(MixParam), '-n-', num2str(n), '.png')])
