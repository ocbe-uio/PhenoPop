function [bic, size_data, num_obs] = compute_bic(data, fval, NP)
%This function computes the Bayes Information Criterion. It takes the
%number of observations in argument data,  the negative log likelihood of a
%model and the number of inferred parameters. 

% returns the bic:
%  bic(1) = total bic,
%  bic(2) = NumParams*log(num_obs),
%  bic(3) = 2*fval.
% the size of the data set, and the number of removed points. 

[NR, NC, NT] = size(data);
size_data=NR*NC*NT;
num_obs = size_data - get_number_discarded_observations(data); %throw out replicates whose initial timepoint is NaN, so do not get counted in LL

bic(1) = NP*log(num_obs) + 2*fval;
bic(2) = NP*log(num_obs);
bic(3) = 2*fval;




