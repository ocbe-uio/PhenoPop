function [L] = vectorized_objective_function_two_noise_levels(Params,DATA,concvec,timevec,NR,ConcT,TimeT)
    %Params are model parameters. DATA is cell counts at each dependent.
    %variable condition. size(DATA) = (NR,NC,NT). 
    %concvec is list of concentration considered.
    %timevec is list of times, and NR is number of replicates at each condition.
    %ConcT is a threshold in concentration.
    %TimeT is a threshold in time.
    %How it works: If the concentration is BELOW the threshold and the time is ABOVE the threshold, then the high noise level (sigH) is used.
    %The idea is that with low concentration the cell count will be large after enough time has passed.
    %To predict, we multiply DATA(k,j,1) by {p*popfunc(Param1,concvec(j),timevec(i)) + (1-p)*popfunc(Param2,concvec(j),timevec(i))}.
    %Predictions are the same for all replicates, but residuals differ between replicates.
    %Params' %uncomment this to print position at every evaluation (does not work because also used in numerical gradient).
    p = Params(1); %mixture parameter
    Param1=Params(2:5); %alpha,b,E,n pop.1
    Param2=Params(6:9); %alpha,b,E,n pop.2
    sigH=Params(10);
    sigL=Params(11);
    if length(Params) ~11
        error('Wrong number of parameters. Use two noise levels.')
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

    % Count the number of High and Low noise levels
    Hi_Lo_indicators = ones(NR,NC,NT-1);
    N_hi = sum(Hi_Lo_indicators(:,concvec <= ConcT,timevec(2:NT) >= TimeT), 'all');
    N_low = NR*NC*(NT-1) - N_hi;
    
    % Make arrays to divide squared residuals elementwise by
    sig_Hi_Lo_array = sigL*ones(NR,NC,NT-1);
    sig_Hi_Lo_array(:,concvec <= ConcT,timevec(2:NT) >= TimeT) = sigH;

    resid_terms = resid.^2 ./ (2 .* sig_Hi_Lo_array.^2);
    resid_sum = sum(resid_terms, 'all', 'omitnan'); % size 1

    % L is the negative normal loglikelihood
    L = resid_sum + (N_hi/2)*log(2*pi*sigH^2) + (N_low/2)*log(2*pi*sigL^2);
end
