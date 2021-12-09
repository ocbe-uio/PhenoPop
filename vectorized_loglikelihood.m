%   ACOV = MLECOV(PARAMS, DATA, 'nloglf',NLOGLF) allows you to define a
%   distribution by its log-likelihood function.  NLOGLF is a function
%   handle, specified using @, that accepts the four input arguments
%      PARAMS - a vector of distribution parameter values
%      DATA   - a vector of data
%      CENS   - a boolean vector of censoring values
%      FREQ   - a vector of integer data frequencies
%   NLOGLF must accept all four arguments even if you do not supply the
%   'censoring' or 'frequency' name/value pairs (see below).  However,
%   NLOGLF can safely ignore its CENS and FREQ arguments in that case.
%   NLOGLF returns a scalar negative log-likelihood value, and optionally,
%   the negative log-likelihood gradient vector (see the 'options'
%   name/value pair below).

function vector_f = vectorized_loglikelihood(Params,data_vector,cens,freq,Conc,Time,R)
    DATA = reshape(data_vector, [R length(Conc) length(Time)]);
    vector_f = ObjectiveFunctionEven(Params,DATA,Conc,Time,R);
end
