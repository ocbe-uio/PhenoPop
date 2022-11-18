function [Conc] = find_the_right_concentrations(max_N_c, sensitive_cell_line, resistant_cell_line, Time)
    % Find concentration locations based on desired growth rates for sensitive or resistant cells

    if max_N_c ~= 17
        error("Error. find_the_right_concentrations not built for Max_Nc ~= 17.")
    end

    alpha1 = sensitive_cell_line(1);
    b1 = sensitive_cell_line(2);
    E1 = sensitive_cell_line(3);
    n1 = sensitive_cell_line(4);
    alpha2 = resistant_cell_line(1);
    b2 = resistant_cell_line(2);
    E2 = resistant_cell_line(3); 
    n2 = resistant_cell_line(4);

    Conc = zeros(1, max_N_c);

    initial_guesses = repmat(E1, 1, max_N_c);
    initial_guesses(9:17) = E2;

    mid_time_index = (length(Time)+1)/2;
    mid_time = Time(mid_time_index);
                                       %Conc,Time,R,MixParam,sensitive,resistant
    hi_count_sensitive = generateDataEven(0,mid_time,1,1,[alpha1,b1,E1,n1],resistant_cell_line);
    lo_count_sensitive = generateDataEven(0,mid_time,1,1,[alpha1+log(b1),b1,E1,n1],resistant_cell_line);
    diff_sensitive = hi_count_sensitive - lo_count_sensitive;

    hi_count_resistant = generateDataEven(0,mid_time,1,0,sensitive_cell_line,[alpha2,b2,E2,n2]);
    lo_count_resistant = generateDataEven(0,mid_time,1,0,sensitive_cell_line,[alpha2+log(b2),b2,E2,n2]);
    diff_resistant = hi_count_resistant - lo_count_resistant;

    desired_counts = zeros(max_N_c);
    %desired_counts(2:8)   = hi_count_sensitive - diff_sensitive * [1/16 1/8 1/4 1/2 3/4 7/8 15/16];
    %desired_counts(10:17) = hi_count_resistant - diff_resistant * [1/16 1/8 1/4 1/2 3/4 7/8 15/16 255/256];
    desired_counts(2:8)   = hi_count_sensitive - diff_sensitive * [1/32 1/16 1/8 1/4 1/2 3/4 7/8];
    desired_counts(10:17) = hi_count_resistant - diff_resistant * [1/32 1/16 1/8 1/4 1/2 3/4 7/8 255/256];
    
    options_fsolve = optimoptions('fsolve','FiniteDifferenceType','central', 'Display', 'off');
    
    % Find conc s.t. at halftime, with single cell lines, we get the desired fraction of cell count
    for ii=2:8
        fun = @(x) generateDataEven(x,mid_time,1,1,sensitive_cell_line,resistant_cell_line) - desired_counts(ii);
        Conc(ii) = fsolve(fun, initial_guesses(ii), options_fsolve);
    end
    
    for ii=10:17
        fun = @(x) generateDataEven(x,mid_time,1,0,sensitive_cell_line,resistant_cell_line) - desired_counts(ii);
        Conc(ii) = fsolve(fun, initial_guesses(ii), options_fsolve);
    end

    % The middle Concentration location is directly between cell lines
    Conc(9) = Conc(5) * sqrt(Conc(12) / Conc(5)); %(E1 - E2)/2;

    %% We moved away from desired rates and instead did desired cell counts. Still found concentrations.
    % desired_rates = zeros(max_N_c);
    % desired_rates(2:8) = alpha1 + log(b1) * [1/16 1/8 1/4 1/2 3/4 7/8 15/16];
    % desired_rates(10:16) = alpha2 + log(b2) * [1/16 1/8 1/4 1/2 3/4 7/8 15/16];
    % desired_rates(17) = alpha2 + log(b2) * 255/256;
    % 
    % options_fsolve = optimoptions('fsolve','FiniteDifferenceType','central', 'Display', 'off');
    % 
    % for ii=2:8
    %     fun = @(x) ratefunc(sensitive_cell_line, x) - desired_rates(ii);
    %     Conc(ii) = fsolve(fun, initial_guesses(ii), options_fsolve);
    % end
    % 
    % for ii=9:17
    %     fun = @(x) ratefunc(resistant_cell_line, x) - desired_rates(ii);
    %     Conc(ii) = fsolve(fun, initial_guesses(ii), options_fsolve);
    % end
    % 
    % % The middle Concentration location is directly between cell lines
    % Conc(9) = Conc(5) * sqrt(Conc(12) / Conc(5)); %(E1 - E2)/2;
end
