function [L, gradL] = vectorized_objective_function(Params,DATA,concvec,timevec,NR)
    %Params are model parameters. DATA is cell counts at each dependent
    %variable condition. size(DATA) = (NR,NC,NT). 
    %concvec is list of concentration considered.
    %timevec is list of times, and NR is number of replicates at each condition.
    %To predict, we multiply DATA(k,j,1) by {p*popfunc(Param1,concvec(j),timevec(i)) + (1-p)*popfunc(Param2,concvec(j),timevec(i))}
    %Predictions are the same for all replicates, but residuals differ between replicates.
    %Params' % print at every evaluation
    p = Params(1); %mixture parameter
    alpha1 = Params(2);
    b1 = Params(3);
    E1 = Params(4);
    n1 = Params(5);
    alpha2 = Params(6);
    b2 = Params(7);
    E2 = Params(8); 
    n2 = Params(9);
    sig=Params(10);
    Param1=Params(2:5); %alpha,b,E,n pop.1
    Param2=Params(6:9); %alpha,b,E,n pop.2
    if length(Params) ~= 10
        error('Wrong number of parameters. Use one noise level.')
    end
    NC=length(concvec);
    NT=length(timevec);

    % Initial cell counts to be used in prediction
    time_zero_DATA = DATA(:,:,1); %size (NR,NC,1)
    rep_time_zero_DATA = repmat(time_zero_DATA,1,1,NT); %size (NR,NC,NT)

    % Prediction at time T by exponential growth with rate affected by hill function of concentration, per cell line 
    pred_multiplier = p*vectorized_popfunc(Param1,concvec,timevec) ...
                        + (1-p)*vectorized_popfunc(Param2,concvec,timevec); % size (1,NC,NT)
    repeated_multiplier = repmat(pred_multiplier,NR,1,1); % size (NR,NC,NT)

    pred_data = rep_time_zero_DATA .* repeated_multiplier; % size (NR,NC,NT)
    %pred_data = max(0, pred_data); % Should we rectify negative values?
    
    resid = DATA(:,:,2:NT) - pred_data(:,:,2:NT); % Residuals at time zero do not enter in the calculation (and are zero anyway)
    sum_squared_resid = sum(resid.^2, 'all', 'omitnan'); % size 1
    resid_term = sum_squared_resid ./ (2*sig^2);

    % L is the negative normal loglikelihood
    L = resid_term + NC*(NT-1)*NR*log(sig);
    %L = resid_term/(NC*(NT-1)*NR) + log(sig);

    % lsqnonlin: 
    %L = reshape(resid, [NC*(NT-1)*NR,1,1]);

    if nargout > 1
        % Analytical gradient
        % (The reason in all this voodoo can be found in a document titled "Analytical gradient")
        % (This document will be kindly provided by Even Moa Myklebust upon request)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % p1, mixture parameter     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        exp_ast_1 = exp(alpha1*timevec)'; % size 1 x Nt
        exp_ast_2 = exp(alpha2*timevec)';
        % The bracket term to the power of t
        upper_term_1 = E1.^n1 + b1*concvec.^n1; % size 1 x Nc, equal to concvec
        lower_term_1 = E1.^n1 + concvec.^n1;
        E_fraction_1 = (upper_term_1 ./ lower_term_1);
        repeated_E_fraction_1 = repmat(E_fraction_1, NT, 1,1);
        Ebf_bracket_to_t_1 = repeated_E_fraction_1 .^ (timevec'); % size Nt by Nc
        % t increases down, d increases to the right. 
        upper_term_2 = E2.^n2 + b2*concvec.^n2; % size 1 x Nc, equal to concvec
        lower_term_2 = E2.^n2 + concvec.^n2;
        E_fraction_2 = (upper_term_2 ./ lower_term_2);
        repeated_E_fraction_2 = repmat(E_fraction_2, NT, 1,1); 
        Ebf_bracket_to_t_2 = repeated_E_fraction_2 .^ (timevec'); % size Nt by Nc

        exp_times_bracket_t_1 = exp_ast_1 .* Ebf_bracket_to_t_1;
        exp_times_bracket_t_2 = exp_ast_2 .* Ebf_bracket_to_t_2;
        dell_f_dell_p = exp_times_bracket_t_1 - exp_times_bracket_t_2; % size Nt by Nc, this is nice to print
        % Interpretation of output: 
        % This is the effect on the f function when increasing p. Can't change it for time zero. 
        % For the lowest concentration (first column), the effect is large at high time points.
        % For the highest concentration at highest time, the effect is small as there are hardly any cells there anyway.

        % General steps
        repeated_dell_f_dell_p = repmat(dell_f_dell_p, 1,1,NR);
        reshape_dell_f_dell_p = permute(repeated_dell_f_dell_p, [3,2,1]); % Transpose and copy to get the golden shape
        resid_times_dell_f_dell_p = resid .* reshape_dell_f_dell_p(:,:,2:NT);
        t_sum_prod_p = sum(resid_times_dell_f_dell_p,3); % shape R by Nc
        product_X0_times_t_sum_prod_p = DATA(:,:,1) .* t_sum_prod_p;
        sum_resid_prod_p = sum(product_X0_times_t_sum_prod_p, 'all');
        partial_wrt_p = sum_resid_prod_p / (sig^2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % alpha1, alpha2            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dell_f_dell_a1 = p * timevec' .* exp_times_bracket_t_1; % size Nt by Nc
        dell_f_dell_a2 = (1-p) * timevec' .* exp_times_bracket_t_2; % size Nt by Nc
        
        % General steps
        repeated_dell_f_dell_a1 = repmat(dell_f_dell_a1, 1,1,NR);
        reshape_dell_f_dell_a1 = permute(repeated_dell_f_dell_a1, [3,2,1]); % Transpose and copy to get the golden shape
        resid_times_dell_f_dell_a1 = resid .* reshape_dell_f_dell_a1(:,:,2:NT);
        t_sum_prod_a1 = sum(resid_times_dell_f_dell_a1,3); % shape R by Nc
        product_X0_times_t_sum_prod_a1 = DATA(:,:,1) .* t_sum_prod_a1;
        sum_resid_prod_a1 = sum(product_X0_times_t_sum_prod_a1, 'all');
        partial_wrt_a1 = sum_resid_prod_a1 / (sig^2);
        
        repeated_dell_f_dell_a2 = repmat(dell_f_dell_a2, 1,1,NR);
        reshape_dell_f_dell_a2 = permute(repeated_dell_f_dell_a2, [3,2,1]); % Transpose and copy to get the golden shape
        resid_times_dell_f_dell_a2 = resid .* reshape_dell_f_dell_a2(:,:,2:NT);
        t_sum_prod_a2 = sum(resid_times_dell_f_dell_a2,3); % shape R by Nc
        product_X0_times_t_sum_prod_a2 = DATA(:,:,1) .* t_sum_prod_a2;
        sum_resid_prod_a2 = sum(product_X0_times_t_sum_prod_a2, 'all');
        partial_wrt_a2 = sum_resid_prod_a2 / (sig^2);
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % b1, b2                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ebf_bracket_to_t_minus_one_1 = repeated_E_fraction_1 .^ (timevec' - 1); % size Nt by Nc
        pt_exp_bracket_t_minus_one_1 = p * timevec' .* exp_ast_1 .* Ebf_bracket_to_t_minus_one_1; % size Nt by Nc
        d_over_E_plus_d_1 = (concvec.^n1 ./ lower_term_1);
        dell_f_dell_b1 = pt_exp_bracket_t_minus_one_1 .* d_over_E_plus_d_1;

        Ebf_bracket_to_t_minus_one_2 = repeated_E_fraction_2 .^ (timevec' - 1); % size Nt by Nc
        pt_exp_bracket_t_minus_one_2 = (1-p) * timevec' .* exp_ast_2 .* Ebf_bracket_to_t_minus_one_2; % size Nt by Nc
        d_over_E_plus_d_2 = (concvec.^n2 ./ lower_term_2);
        dell_f_dell_b2 = pt_exp_bracket_t_minus_one_2 .* d_over_E_plus_d_2;

        % General steps
        repeated_dell_f_dell_b1 = repmat(dell_f_dell_b1, 1,1,NR);
        reshape_dell_f_dell_b1 = permute(repeated_dell_f_dell_b1, [3,2,1]); % Transpose and copy to get the golden shape
        resid_times_dell_f_dell_b1 = resid .* reshape_dell_f_dell_b1(:,:,2:NT);
        t_sum_prod_b1 = sum(resid_times_dell_f_dell_b1,3); % shape R by Nc
        product_X0_times_t_sum_prod_b1 = DATA(:,:,1) .* t_sum_prod_b1;
        sum_resid_prod_b1 = sum(product_X0_times_t_sum_prod_b1, 'all');
        partial_wrt_b1 = sum_resid_prod_b1 / (sig^2);
        
        repeated_dell_f_dell_b2 = repmat(dell_f_dell_b2, 1,1,NR);
        reshape_dell_f_dell_b2 = permute(repeated_dell_f_dell_b2, [3,2,1]); % Transpose and copy to get the golden shape
        resid_times_dell_f_dell_b2 = resid .* reshape_dell_f_dell_b2(:,:,2:NT);
        t_sum_prod_b2 = sum(resid_times_dell_f_dell_b2,3); % shape R by Nc
        product_X0_times_t_sum_prod_b2 = DATA(:,:,1) .* t_sum_prod_b2;
        sum_resid_prod_b2 = sum(product_X0_times_t_sum_prod_b2, 'all');
        partial_wrt_b2 = sum_resid_prod_b2 / (sig^2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % E1, E2                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        upper_E1 = n1 .* E1.^(n1 - 1) .* (concvec .^n1) .* (1 - b1);
        frac_n_E_d_one_min_b_1 = upper_E1 ./ (lower_term_1 .^2);
        dell_f_dell_E1 = pt_exp_bracket_t_minus_one_1 .* frac_n_E_d_one_min_b_1;

        upper_E2 = n2 .* E2.^(n2 - 1) .* (concvec .^n2) .* (1 - b2);
        frac_n_E_d_one_min_b_2 = upper_E2 ./ (lower_term_2 .^2);
        dell_f_dell_E2 = pt_exp_bracket_t_minus_one_2 .* frac_n_E_d_one_min_b_2;

        % General steps
        repeated_dell_f_dell_E1 = repmat(dell_f_dell_E1, 1,1,NR);
        reshape_dell_f_dell_E1 = permute(repeated_dell_f_dell_E1, [3,2,1]); % Transpose and copy to get the golden shape
        resid_times_dell_f_dell_E1 = resid .* reshape_dell_f_dell_E1(:,:,2:NT);
        t_sum_prod_E1 = sum(resid_times_dell_f_dell_E1,3); % shape R by Nc
        product_X0_times_t_sum_prod_E1 = DATA(:,:,1) .* t_sum_prod_E1;
        sum_resid_prod_E1 = sum(product_X0_times_t_sum_prod_E1, 'all');
        partial_wrt_E1 = sum_resid_prod_E1 / (sig^2);

        repeated_dell_f_dell_E2 = repmat(dell_f_dell_E2, 1,1,NR);
        reshape_dell_f_dell_E2 = permute(repeated_dell_f_dell_E2, [3,2,1]); % Transpose and copy to get the golden shape
        resid_times_dell_f_dell_E2 = resid .* reshape_dell_f_dell_E2(:,:,2:NT);
        t_sum_prod_E2 = sum(resid_times_dell_f_dell_E2,3); % shape R by Nc
        product_X0_times_t_sum_prod_E2 = DATA(:,:,1) .* t_sum_prod_E2;
        sum_resid_prod_E2 = sum(product_X0_times_t_sum_prod_E2, 'all');
        partial_wrt_E2 = sum_resid_prod_E2 / (sig^2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % n1, n2                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For concentration = 0, the derivative is 0.
        % The first concentration must always be 0.
        upper_n1_nonzero = E1.^n1 .* (b1 - 1) .* concvec(2:NT).^n1 .* log( concvec(2:NT) ./ E1);
        upper_n1 = [0 upper_n1_nonzero]; % shape 1 x Nc
        frac_E_b_minus_one_d_log_dE_1 = upper_n1 ./ (lower_term_1 .^2);
        dell_f_dell_n1 = pt_exp_bracket_t_minus_one_1 .* frac_E_b_minus_one_d_log_dE_1;

        upper_n2_nonzero = E2.^n2 .* (b2 - 1) .* concvec(2:NT).^n2 .* log( concvec(2:NT) ./ E2);
        upper_n2 = [0 upper_n2_nonzero]; % shape 1 x Nc
        frac_E_b_minus_one_d_log_dE_2 = upper_n2 ./ (lower_term_2 .^2);
        dell_f_dell_n2 = pt_exp_bracket_t_minus_one_2 .* frac_E_b_minus_one_d_log_dE_2;
        
        % General steps
        repeated_dell_f_dell_n1 = repmat(dell_f_dell_n1, 1,1,NR);
        reshape_dell_f_dell_n1 = permute(repeated_dell_f_dell_n1, [3,2,1]); % Transpose and copy to get the golden shape
        resid_times_dell_f_dell_n1 = resid .* reshape_dell_f_dell_n1(:,:,2:NT);
        t_sum_prod_n1 = sum(resid_times_dell_f_dell_n1,3); % shape R by Nc
        product_X0_times_t_sum_prod_n1 = DATA(:,:,1) .* t_sum_prod_n1;
        sum_resid_prod_n1 = sum(product_X0_times_t_sum_prod_n1, 'all');
        partial_wrt_n1 = sum_resid_prod_n1 / (sig^2);

        repeated_dell_f_dell_n2 = repmat(dell_f_dell_n2, 1,1,NR);
        reshape_dell_f_dell_n2 = permute(repeated_dell_f_dell_n2, [3,2,1]); % Transpose and copy to get the golden shape
        resid_times_dell_f_dell_n2 = resid .* reshape_dell_f_dell_n2(:,:,2:NT);
        t_sum_prod_n2 = sum(resid_times_dell_f_dell_n2,3); % shape R by Nc
        product_X0_times_t_sum_prod_n2 = DATA(:,:,1) .* t_sum_prod_n2;
        sum_resid_prod_n2 = sum(product_X0_times_t_sum_prod_n2, 'all');
        partial_wrt_n2 = sum_resid_prod_n2 / (sig^2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sig                       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        partial_wrt_sig = sum_squared_resid / (sig^3) - NC*(NT-1)*NR / sig;

        analytical_gradient = [partial_wrt_p, partial_wrt_a1, partial_wrt_b1, partial_wrt_E1, partial_wrt_n1, partial_wrt_a2, partial_wrt_b2, partial_wrt_E2, partial_wrt_n2, partial_wrt_sig];
        gradL = - analytical_gradient;
        %gradL = - analytical_gradient / (NC*(NT-1)*NR);
    end
end

