% Function returning upper and lower bounds
function [lb, ub, A_inequality, b_inequality, nonlcon]=optimization_constraints_real_data(no_populations, INFER_SIGMA, USE_TWO_NOISE_LEVELS, lower_limit_mixparam, lower_limit_E_ratio, lb_all, ub_all)
    % xx = [alpha1 b1 E1 n1 p1, alpha2 b2 E2 n2 p2 ... alpha_k b_k E_k n_k sig]
    %alpha_m = Params(5*m-4);
    %b_m = Params(5*m-3);
    %E_m = Params(5*m-2);
    %n_m = Params(5*m-1);
    %p_m = Params(5*m); %except p_k which is 1 minus sum(p_m)
    %sig=Params(end);
    if INFER_SIGMA
        lb=zeros(5*no_populations,1);
        lb(length(lb))=1e-6; %To avoid Inf or NaN in log() terms we require positive Noise
        if USE_TWO_NOISE_LEVELS
            lb=zeros(5*no_populations + 1,1);
            lb(length(lb))     = lb_all(5); % 1e-6; %To avoid Inf or NaN in log() terms we require positive Noise
            lb(length(lb) - 1) = lb_all(5); % 1e-6; %Also for sigH
        end
    else
        lb=zeros(5*no_populations - 1,1);
        % Drop Noise parameter
    end
    lb(1:5:5*no_populations-4)= lb_all(1); %-0.1; % alpha
    lb(2:5:5*no_populations-3)= lb_all(2); %0; % b
    lb(3:5:5*no_populations-2)= lb_all(3); %-6; %For log scale sampling in E
    lb(4:5:5*no_populations-1)= lb_all(4); % n
    lb(5:5:5*no_populations-5)= lower_limit_mixparam; % Mixture parameter required to be more than zero
    
    ub=1000*ones(length(lb),1);
    ub(1:5:5*no_populations-4)= ub_all(1); %0.1; % alpha
    ub(2:5:5*no_populations-3)= ub_all(2); %1; % b
    ub(3:5:5*no_populations-2)= ub_all(3); %4; % E
    ub(4:5:5*no_populations-1)= ub_all(4); %100; % n
    ub(5:5:5*no_populations-5)= 1; % Mixture parameter
    if INFER_SIGMA
        ub(length(ub))= ub_all(5); % 5000; % Noise
        if USE_TWO_NOISE_LEVELS
            ub(length(ub) - 1)= ub_all(5); % 5000; % Also for sigH
        end
    end

    % The k-1 mixture parameters must sum to less than 1. A*x < b where A picks out mixture params
    if INFER_SIGMA
        A_inequality = zeros(1,5*no_populations);
        if USE_TWO_NOISE_LEVELS
            A_inequality = zeros(1,5*no_populations + 1); % Also for sigH
        end
    else
        A_inequality = zeros(1,5*no_populations - 1); % Drop Noise parameter
    end
    % Making mixture parameters sum to 1
    A_inequality(5:5:5*no_populations-5) = 1;
    b_inequality = 1-lower_limit_mixparam;

    nonlcon = @(x)E_constraint(x,lower_limit_E_ratio,no_populations);

end
