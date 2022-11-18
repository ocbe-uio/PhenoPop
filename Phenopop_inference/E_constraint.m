function [c,ceq] = E_constraint(x,lower_limit_E_ratio,no_populations)
    % Forces populations apart by requiring big enough minimum distance between log(E) values:
    E_values = x(3:5:5*no_populations-2);
    if length(E_values) == 1
        c = [];
        ceq = [];
    else
        % min_i,j {E_i/E_j} >= 2
        % equivalently: 
        % -(logE_i - logE_j)^2 + (log(2))^2 <= 0

        % shape = N observations by P=1 dimensions
        logE = log(E_values);
        alldist = pdist(logE);
        mindist = min(alldist,[],'all');
        
        c = -mindist^2 + (log(lower_limit_E_ratio))^2;
        ceq = [];
    end    
end 
