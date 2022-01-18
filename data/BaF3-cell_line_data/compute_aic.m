function aic = compute_aic(fval, NP)
%data must be in (rep, conc, time) dimensions. compute the Akaike
%Information Criterion of the negative log likelihood
%2*NP + log(fval)
%Here NP is the number of parameters inferred. 

aic = 2*NP + 2*fval;

