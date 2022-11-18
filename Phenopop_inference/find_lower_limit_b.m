function [lower_bound_b] = find_lower_limit_b(remaining_cell_fraction, upper_limit_alpha, delta_t)
    lower_bound_b = exp(log(remaining_cell_fraction)./delta_t - upper_limit_alpha);
end 
