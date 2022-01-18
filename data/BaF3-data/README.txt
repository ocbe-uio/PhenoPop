README - 5-20-21

Goal: The files complied in this folder are used in the mixtures project, where the goal is to infer a model to describe the data - MONOCLONAL_DATA and Mixure_Data - with the correct number of subpopulation, ie a one and two population model respectively. 

The files in this folder are authored by Jasmine Noory. Many codes are adapted scripts authored by Even Moa as well. 


Configuration: This folder is self-contained. As long as files are run within this directory, everything should run properly.



Description of Data:
The main data we use for inference is MONOCLONAL_DATA.mat and Mixture_Data.mat. These  files each contain 4 variables saving results from Dagim’s experiments. 

The monoclonal data variables track raw cell counts at 3 hr increments between times 9 and 48 hours for 11 experimental concentrations with 7 biological replicates.
The datasets are named after the intended initial cell count and cell type. 
	
RESISTANT_250_BF,	     RESISTANT_500_BF,       SENSITIVE_500_B,       SENSITIVE_1000_BF

The mixture data set is the same experimental setup, with 14 biological replicates. The datasets are named after the intended mixture ratio, parts sensitive to parts resistant

BF_11		BF_21		BF_12		BF_41

The generation and cleaning of these data sets can be found in Create_and_Clean_Datasets


Description of inference script
The main inference script is Inference_Script. We infer governing growth parameters [alpha, b, E, n, p] (p mixture) for every subpopulation and two noise levels using MLE.
We minimize the negative log likelihood computed with the function vectorized_objective_function_two_noise_levels_k_subpopulations.m.

Results are saved in the .mat files Results_Monoclonal_max4pop_infer and Results_Mixture_max4pop_infer. The saved variables are 
xfinals : The inferred models in array of size (model params, number of inferred models, number of datasets) The convention for 
the order of inferred parameters is 
[alpha1, b1, E1, n1, p1, alpha2, b2, E2, n2, p2, ... , alphak, bk, Ek, nk, sigH, sigL].
fvals : The fvals of each inferred model. Array of size(1, num of inferred models, number of datasets) 
Model score arrays: contain information criteria scores based on fvals and number of parameters inferred. 
Bic_array
bic_squared_num_params_array
 aic arrays.

Description of Plotting Files
There are files to plot the best_inferred_rates, plot elbow plots using the information criterion scores and negative log likelihood,
and plot_population_predictions. These files produce plots in loops with the convention that the datasets are enumerated. 

To run these files, several booleans must be defined initially to determine which datasets we work with.

Supporting Functions
compute_aic
compute_bic
get_fitted_rates
get_number_discarded_observations
ratefunc
vectorized_inference_k_subpopulations_two_noise_levels
vectorized_popfunc
plot_linear_timecourse

