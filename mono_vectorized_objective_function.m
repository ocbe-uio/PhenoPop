function L=mono_vectorized_objective_function(Params,DATA,concvec,timevec,NR)
    %Params are model parameters. DATA is cell counts at each dependent
    %variable condition. size(DATA) = (NR,NC,NT). 
    %ConcVector is list of concentration considered.
    %TimeVector is list of times, and NR is number of replicates at each condition.
    %To predict, we multiply DATA(k,j,1) by {p*popfunc(Param1,concvec(j),timevec(i)) + (1-p)*popfunc(Param2,concvec(j),timevec(i))}
    %Predictions are the same for all replicates, but residuals differ between replicates.
%%%%    p=1 %%%%% Params(1); % Mixture parameter
    Param1=Params(1:4); %alpha,b,E,n pop.1
    sig=Params(5);
    if length(Params) > 5
        error('There are too many parameters being passed to ObjectiveFunction!')
    end
    NC=length(concvec);
    NT=length(timevec);

    % Initial cell counts to be used in prediction
    time_zero_DATA = DATA(:,:,1); %size (NR,NC,1)
    rep_time_zero_DATA = repmat(time_zero_DATA,1,1,NT); %size (NR,NC,NT)

    % Prediction at time T by exponential growth with rate affected by hill function of concentration, per cell line 
    pred_multiplier = 1.0 * vectorized_popfunc(Param1,concvec,timevec); % size (1,NC,NT)
    repeated_multiplier = repmat(pred_multiplier,NR,1,1); % size (NR,NC,NT)

    pred_data = rep_time_zero_DATA .* repeated_multiplier; % size (NR,NC,NT)
    %pred_data = max(0, pred_data);
    
    resid = DATA(:,:,2:NT) - pred_data(:,:,2:NT); % Cell counts at time zero do not enter in the residual calculation
    %DATA(:,:,2:NT)
    %pred_data(:,:,2:NT)
    sum_squared_resid = sum(resid.^2, 'all', 'omitnan'); % size 1
    resid_term = sum_squared_resid ./ (2*sig^2);
    L = resid_term + NC*NT*NR*log(sig);
    end
    
    