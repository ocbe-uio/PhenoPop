case_array = ["TREAT E11", "TREAT E14", "TREAT E19", "TREAT E1975", "TREAT E41", "TREAT E827", "TREAT P11", "TREAT P14", "TREAT P19", "TREAT P1975", "TREAT P41", "TREAT P827"];
savestr_array = ["TREAT_E_11", "TREAT_E_14", "TREAT_E_19", "TREAT_E_1975", "TREAT_E_41", "TREAT_E_827", "TREAT_P_11", "TREAT_P_14", "TREAT_P_19", "TREAT_P_1975", "TREAT_P_41", "TREAT P827"];
max_no_populations = 5 % The highest number of populations to try fitting to the data.
num_optim = 3000 % Number of starting points in maximum likelihood optimization.
INFER_SIGMA = true % If true, sigma is estimated by MLE. If false, the true value is given.
USE_TWO_NOISE_LEVELS = true
gr50values = zeros(length(case_array), max_no_populations, max_no_populations);

%%%%%%%%% Find max values %%%%%%%%%%
max_neg_ll = 0;
min_neg_ll = Inf;
max_difference = 0;

for casenumber = 1:length(case_array) %%%%%%%%%%%%%%%% parfor loop over drugs
casename = case_array(casenumber);
savestr1 = savestr_array(casenumber);
savestr = strcat(savestr1, '-num_optim-', num2str(num_optim));

load([strcat('./data/shannon_data/negative_loglikelihood_values-', savestr, '.mat')])

Time = [0 24 48 72];
Conc = [0 0.1 1 10];

for no_populations = 1:max_no_populations % how many populations we infer
for jj = 1:no_populations % which of those population we find the gr50 for

%no_populations
%jj
%if no_populations ~= 2 || (jj ~= 2 && jj ~= 5)
jj_params = squeeze(x_finals_temp(1,no_populations,5*jj-4:5*jj-1))';
gr50values(casenumber, no_populations, jj) = find_gr50(jj_params, Conc);
%end 

end
end

casename
gr50values(casenumber,:,:)
end 
