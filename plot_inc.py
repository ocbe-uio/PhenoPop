from moodboard import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from matplotlib.patches import Patch
import scipy.stats

# Plotting results from increasing values, with errorbars
base_case = 3
N_repetitions = 30
num_optim = 1000
seednr_1 = 43
generation_number = 5
Noise = 100
base_case = 3
MixParams = [0.4]
base_rate_sensitive = 0.03
base_rate_middle = 0.03
base_rate_resistant = 0.03
b1 = 0.3
b2 = 0.4
b3 = 0.5
E_sensitive = 1.0e-4
E_middle = 1.0e-2
E_resistant = 1.0e-1
n = 3
#sensitive_cell_line = [base_rate_sensitive, b1, E_sensitive, n]
#middle_cell_line = [base_rate_middle, b2, E_middle, n]
#resistant_cell_line = [base_rate_resistant, b3, E_resistant, n]
#Params = [sensitive_cell_line, resistant_cell_line]
#gridlength = 4
#N_c_values = 2.^(1:gridlength)+1 %[3 5 9 17]
#N_t_values = 2.^(1:gridlength)+1 %[3 5 9 17]
#R_values = 2.^(1:gridlength)+1 %[3 5 9 17]
#max_N_c = max(N_c_values)
#max_N_t = max(N_t_values)
#max_R = max(R_values)

# Load data
RCT_array = ['R', 'C', 'T']
for ii in range(len(RCT_array)):
    change_value = RCT_array[ii]
    filename = './data/difference/inc-' + str(change_value) + '-gen-' + str(generation_number) + 'base_case' + str(base_case) + '-num_optim-' + str(num_optim) + '-N-rep-' + str(N_repetitions) + '-Noise-' + str(Noise) + '-seed-' + str(seednr_1) + '-MixParams-' + str(MixParams[0]) + '-n-' + str(n) + '-b1-' + str(b1) + '-b2-' + str(b2) + '-b3-' + str(b3) + '-Esens-' + str(E_sensitive) + '-Eres-' + str(E_resistant) + '.mat'
    mat_contents = scipy.io.loadmat(filename)
    #print(mat_contents)
    if change_value == "R":
        mean_errors_mixture_fraction_r = mat_contents['mean_errors_mixture_fraction'][0]
        mean_logmaxratio_GR50_1_r = mat_contents['mean_logmaxratio_GR50_1'][0]
        mean_logmaxratio_GR50_2_r = mat_contents['mean_logmaxratio_GR50_2'][0]
        mean_errors_alpha_1_r = mat_contents['mean_errors_alpha_1'][0]
        mean_errors_alpha_2_r = mat_contents['mean_errors_alpha_2'][0]
        std_errors_mixture_fraction_r = mat_contents['std_errors_mixture_fraction'][0]
        std_logmaxratio_GR50_1_r = mat_contents['std_logmaxratio_GR50_1'][0]
        std_logmaxratio_GR50_2_r = mat_contents['std_logmaxratio_GR50_2'][0]
        std_errors_alpha_1_r = mat_contents['std_errors_alpha_1'][0]
        std_errors_alpha_2_r = mat_contents['std_errors_alpha_2'][0]
    elif change_value == "C":
        mean_errors_mixture_fraction_nc = mat_contents['mean_errors_mixture_fraction'][0]
        mean_logmaxratio_GR50_1_nc = mat_contents['mean_logmaxratio_GR50_1'][0]
        mean_logmaxratio_GR50_2_nc = mat_contents['mean_logmaxratio_GR50_2'][0]
        mean_errors_alpha_1_nc = mat_contents['mean_errors_alpha_1'][0]
        mean_errors_alpha_2_nc = mat_contents['mean_errors_alpha_2'][0]        
        std_errors_mixture_fraction_nc = mat_contents['std_errors_mixture_fraction'][0]
        std_logmaxratio_GR50_1_nc = mat_contents['std_logmaxratio_GR50_1'][0]
        std_logmaxratio_GR50_2_nc = mat_contents['std_logmaxratio_GR50_2'][0]
        std_errors_alpha_1_nc = mat_contents['std_errors_alpha_1'][0]
        std_errors_alpha_2_nc = mat_contents['std_errors_alpha_2'][0]
    elif change_value == "T":
        mean_errors_mixture_fraction_nt = mat_contents['mean_errors_mixture_fraction'][0]
        mean_logmaxratio_GR50_1_nt = mat_contents['mean_logmaxratio_GR50_1'][0]
        mean_logmaxratio_GR50_2_nt = mat_contents['mean_logmaxratio_GR50_2'][0]
        mean_errors_alpha_1_nt = mat_contents['mean_errors_alpha_1'][0]
        mean_errors_alpha_2_nt = mat_contents['mean_errors_alpha_2'][0]
        std_errors_mixture_fraction_nt = mat_contents['std_errors_mixture_fraction'][0]
        std_logmaxratio_GR50_1_nt = mat_contents['std_logmaxratio_GR50_1'][0]
        std_logmaxratio_GR50_2_nt = mat_contents['std_logmaxratio_GR50_2'][0]
        std_errors_alpha_1_nt = mat_contents['std_errors_alpha_1'][0]
        std_errors_alpha_2_nt = mat_contents['std_errors_alpha_2'][0]

# Test for significance
# https://se.mathworks.com/help/stats/anova1.html h = ttest2(x,y,'alpha',0.05/27)

dof = N_repetitions - 1 
# For 30 repetitions there are 29 degrees of freedom. 
# Two-sided test with alpha 0.05: Critical value of 2.045. Apply to results below
#critical_value = scipy.stats.t.ppf(q=1-.05/2,df=29)
#print(critical_value)

#find T critical value #https://www.khanacademy.org/math/ap-statistics/xfb5d8e68:inference-quantitative-means/two-sample-t-test-means/v/two-sample-t-test-for-difference-of-means
bonferroni_corrected_significance_level = 0.05/27
critical_value = scipy.stats.t.ppf(q=1-bonferroni_corrected_significance_level/2,df=29)
print(critical_value)

# Mixture fraction
# R vs Nc
joint_std_estimate = np.sqrt(std_errors_mixture_fraction_r**2/N_repetitions + std_errors_mixture_fraction_nc**2/N_repetitions)
t_mix_r_nc = (mean_errors_mixture_fraction_r - mean_errors_mixture_fraction_nc) / joint_std_estimate
print("\nt_mix_r_nc: ", t_mix_r_nc)
print("Significant?: ", abs(t_mix_r_nc) > critical_value)
# R vs Nt
joint_std_estimate = np.sqrt(std_errors_mixture_fraction_r**2/N_repetitions + std_errors_mixture_fraction_nt**2/N_repetitions)
t_mix_r_nt = (mean_errors_mixture_fraction_r - mean_errors_mixture_fraction_nt) / joint_std_estimate
print("\nt_mix_r_nt: ", t_mix_r_nt)
print("Significant?: ", abs(t_mix_r_nt) > critical_value)
# Nc vs Nt
joint_std_estimate = np.sqrt(std_errors_mixture_fraction_nc**2/N_repetitions + std_errors_mixture_fraction_nt**2/N_repetitions)
t_mix_nc_nt = (mean_errors_mixture_fraction_nc - mean_errors_mixture_fraction_nt) / joint_std_estimate
print("\nt_mix_nc_nt: ", t_mix_nc_nt)
print("Significant?: ", abs(t_mix_nc_nt) > critical_value)
# GR50 clone 1
# R vs Nc
joint_std_estimate = np.sqrt(std_logmaxratio_GR50_1_r**2/N_repetitions + std_logmaxratio_GR50_1_nc**2/N_repetitions)
t_GR50_1_r_nc = (mean_logmaxratio_GR50_1_r - mean_logmaxratio_GR50_1_nc) / joint_std_estimate
print("\nt_GR50_1_r_nc: ", t_GR50_1_r_nc)
print("Significant?: ", abs(t_GR50_1_r_nc) > critical_value)
# R vs Nt
joint_std_estimate = np.sqrt(std_logmaxratio_GR50_1_r**2/N_repetitions + std_logmaxratio_GR50_1_nt**2/N_repetitions)
t_GR50_1_r_nt = (mean_logmaxratio_GR50_1_r - mean_logmaxratio_GR50_1_nt) / joint_std_estimate
print("\nt_GR50_1_r_nt: ", t_GR50_1_r_nt)
print("Significant?: ", abs(t_GR50_1_r_nt) > critical_value)
# Nc vs Nt
joint_std_estimate = np.sqrt(std_logmaxratio_GR50_1_nc**2/N_repetitions + std_logmaxratio_GR50_1_nt**2/N_repetitions)
t_GR50_1_nc_nt = (mean_logmaxratio_GR50_1_nc - mean_logmaxratio_GR50_1_nt) / joint_std_estimate
print("\nt_GR50_1_nc_nt: ", t_GR50_1_nc_nt)
print("Significant?: ", abs(t_GR50_1_nc_nt) > critical_value)
# GR50 clone 1
# R vs Nc
joint_std_estimate = np.sqrt(std_logmaxratio_GR50_2_r**2/N_repetitions + std_logmaxratio_GR50_2_nc**2/N_repetitions)
t_GR50_2_r_nc = (mean_logmaxratio_GR50_2_r - mean_logmaxratio_GR50_2_nc) / joint_std_estimate
print("\nt_GR50_2_r_nc: ", t_GR50_2_r_nc)
print("Significant?: ", abs(t_GR50_2_r_nc) > critical_value)
# R vs Nt
joint_std_estimate = np.sqrt(std_logmaxratio_GR50_2_r**2/N_repetitions + std_logmaxratio_GR50_2_nt**2/N_repetitions)
t_GR50_2_r_nt = (mean_logmaxratio_GR50_2_r - mean_logmaxratio_GR50_2_nt) / joint_std_estimate
print("\nt_GR50_2_r_nt: ", t_GR50_2_r_nt)
print("Significant?: ", abs(t_GR50_2_r_nt) > critical_value)
# Nc vs Nt
joint_std_estimate = np.sqrt(std_logmaxratio_GR50_2_nc**2/N_repetitions + std_logmaxratio_GR50_2_nt**2/N_repetitions)
t_GR50_2_nc_nt = (mean_logmaxratio_GR50_2_nc - mean_logmaxratio_GR50_2_nt) / joint_std_estimate
print("\nt_GR50_2_nc_nt: ", t_GR50_2_nc_nt)
print("Significant?: ", abs(t_GR50_2_nc_nt) > critical_value)

#mean_errors_mixture_fraction_r =    np.array([1,1,1,1])*1
#mean_errors_mixture_fraction_nc =   np.array([1,1,1,1])*2
#mean_errors_mixture_fraction_nt =   np.array([1,1,1,1])*3
#std_errors_mixture_fraction_r =     np.array([1,1,1,1])*0.1
#std_errors_mixture_fraction_nc =    np.array([1,1,1,1])*0.2
#std_errors_mixture_fraction_nt =    np.array([1,1,1,1])*0.3

#def make_inc_plot(base_case, change_value, values_p, values_gr50_1, values_gr50_2, std_p, std_GR50_1, std_GR50_2, title, savename):
#    # 95 % confidence intervals calculated using the t-distribution with 29 degrees of freedom. 
#    fig, ax1 = plt.subplots()
#    ax2 = ax1.twinx()
#    # GR50s
#    ax2.errorbar([3,5,9,17], values_gr50_1, yerr=2.045*std_GR50_1/np.sqrt(N_repetitions), label=r"GR50 clone $\rm 1$", color=paul_tol[2], linewidth=2, capsize=6, zorder=1)
#    ax2.errorbar([3,5,9,17], values_gr50_2, yerr=2.045*std_GR50_2/np.sqrt(N_repetitions), label=r"GR50 clone $\rm 2$", color=paul_tol[2], linewidth=2, linestyle="--", capsize=6, zorder=1)
#    ax2.set_ylabel("Max log ratio GR50", color=paul_tol[2])
#    # Mixture fraction
#    ax1.errorbar([3,5,9,17], values_p, yerr=2.045*std_p/np.sqrt(N_repetitions), label="Mixture fraction", color=paul_tol[1], linewidth=2, capsize=6, zorder=30)
#    ax1.set_ylabel("Absolute error in mixture fraction", color=paul_tol[1])
#    ax1.set_xlabel(change_value)
#    ax1.set_xticks([3,5,9,17])
#    #plt.vlines(base_case, ymin=0, ymax=1.05*max(max(values_r), max(values_nc), max(values_nt)), color='k', linewidth=2)
#    plt.title(title)
#    plt.tight_layout()
#    # ask matplotlib for the plotted objects and their labels
#    lines, labels = ax1.get_legend_handles_labels()
#    lines2, labels2 = ax2.get_legend_handles_labels()
#    ax2.legend(lines + lines2, labels + labels2, loc='upper right', frameon=False)
#    plt.savefig(savename)
#    #plt.show()
#
## Fig 1 increasing R
#change_value = "Replicates"
#values_p = mean_errors_mixture_fraction_r
#values_gr50_1 = mean_logmaxratio_GR50_1_r
#values_gr50_2 = mean_logmaxratio_GR50_2_r
#std_p = std_errors_mixture_fraction_r
#std_GR50_1 = std_logmaxratio_GR50_1_r
#std_GR50_2 = std_logmaxratio_GR50_2_r
#title = r"$N_c = 3, N_t = 3$"
#savename = "./article/incplots/inc-1-r.png"
#make_inc_plot(base_case, change_value, values_p, values_gr50_1, values_gr50_2, std_p, std_GR50_1, std_GR50_2, title, savename)
#
## Fig 2 increasing Nc
#change_value = "Observed concentrations"
#values_p = mean_errors_mixture_fraction_nc
#values_gr50_1 = mean_logmaxratio_GR50_1_nc
#values_gr50_2 = mean_logmaxratio_GR50_2_nc
#std_p = std_errors_mixture_fraction_nc
#std_GR50_1 = std_logmaxratio_GR50_1_nc
#std_GR50_2 = std_logmaxratio_GR50_2_nc
#title = r"$N_t = 3, R = 3$"
#savename = "./article/incplots/inc-1-nc.png"
#make_inc_plot(base_case, change_value, values_p, values_gr50_1, values_gr50_2, std_p, std_GR50_1, std_GR50_2, title, savename)
#
## Fig 3 increasing Nt
#change_value = "Observed timepoints"
#values_p = mean_errors_mixture_fraction_nt
#values_gr50_1 = mean_logmaxratio_GR50_1_nt
#values_gr50_2 = mean_logmaxratio_GR50_2_nt
#std_p = std_errors_mixture_fraction_nt
#std_GR50_1 = std_logmaxratio_GR50_1_nt
#std_GR50_2 = std_logmaxratio_GR50_2_nt
#title = r"$N_c = 3, R = 3$"
#savename = "./article/incplots/inc-1-nt.png"
#make_inc_plot(base_case, change_value, values_p, values_gr50_1, values_gr50_2, std_p, std_GR50_1, std_GR50_2, title, savename)


# c0c0c0 # 868686
def make_inc_plot(base_case, values_r, values_nc, values_nt, std_r, std_nc, std_nt, ylabel, savename):
    # 95 % confidence intervals calculated using the t-distribution with 29 degrees of freedom. 
    fig = plt.figure()
    plt.errorbar([3,5,9,17], values_r, yerr=2.045*std_r/np.sqrt(N_repetitions), color=paul_tol[5], capsize=6) #, label=r"$\rm Changing \enspace R$"
    plt.errorbar([3,5,9,17], values_nt, yerr=2.045*std_nt/np.sqrt(N_repetitions), color=paul_tol[4], capsize=6) #, label=r"$\rm Changing \enspace N_t$"
    plt.errorbar([3,5,9,17], values_nc, yerr=2.045*std_nc/np.sqrt(N_repetitions), color=paul_tol[0], capsize=6) #, label=r"$\rm Changing \enspace N_c$"
    #plt.vlines(base_case, ymin=0, ymax=1.05*max(max(values_r), max(values_nc), max(values_nt)), color='k', linewidth=2)
    #plt.xlabel(r"$\rm R, N_c, N_t $")
    plt.ylabel(ylabel)
    legend_elements = []
    legend_elements = legend_elements + [Patch(label=r"$\rm N_c=3$, $\rm N_t=3$. R value follows x axis", color=paul_tol[5])]
    legend_elements = legend_elements + [Patch(label=r"$\rm N_c=3$, $\rm R=3$. $\rm N_t$ value follows x axis", color=paul_tol[4])]
    legend_elements = legend_elements + [Patch(label=r"$\rm N_t=3$, $\rm R=3$. $\rm N_c$ value follows x axis", color=paul_tol[0])]
    plt.legend(handles=legend_elements, bbox_to_anchor=(0.2,-0.04,1,-0.04), loc='upper left', ncol=1, mode="expand", borderaxespad=0.0, frameon=False) # bbox:(x y width height) # upper left means upper left position is fixed. 
    #plt.legend(bbox_to_anchor=(0,-0.1,1,-0.1), loc='upper left', borderaxespad=0.0, frameon=False) # bbox:(x y width height) # upper left means upper left position is fixed. 
    plt.ylim(ymin=0)
    plt.xticks([3,5,9,17])
    plt.tight_layout()
    plt.savefig(savename)
    #plt.show()

# Fig 1: Error in mixture fraction (population 1)
values_r = mean_errors_mixture_fraction_r
values_nc = mean_errors_mixture_fraction_nc
values_nt = mean_errors_mixture_fraction_nt
std_r = std_errors_mixture_fraction_r
std_nc = std_errors_mixture_fraction_nc
std_nt = std_errors_mixture_fraction_nt
title = "Mixture fraction"
ylabel = "Absolute error in mixture fraction"
savename = "./article/incplots/inc-1-mixture-fraction.png"
make_inc_plot(base_case, values_r, values_nc, values_nt, std_r, std_nc, std_nt, ylabel, savename)

# Fig 2 GR50 clone 1
title = "GR50 clone 1"
values_r = mean_logmaxratio_GR50_1_r
values_nc = mean_logmaxratio_GR50_1_nc
values_nt = mean_logmaxratio_GR50_1_nt
std_r = std_logmaxratio_GR50_1_r
std_nc = std_logmaxratio_GR50_1_nc
std_nt = std_logmaxratio_GR50_1_nt
ylabel = "Max GR50 log ratio clone 1"
savename = "./article/incplots/inc-2-GR50-clone-1.png"
make_inc_plot(base_case, values_r, values_nc, values_nt, std_r, std_nc, std_nt, ylabel, savename)

# Fig 3 GR50 clone 2
title = "GR50 clone 2"
values_r = mean_logmaxratio_GR50_2_r
values_nc = mean_logmaxratio_GR50_2_nc
values_nt = mean_logmaxratio_GR50_2_nt
std_r = std_logmaxratio_GR50_2_r
std_nc = std_logmaxratio_GR50_2_nc
std_nt = std_logmaxratio_GR50_2_nt
ylabel = "Max GR50 log ratio clone 2"
savename = "./article/incplots/inc-3-GR50-clone-2.png"
make_inc_plot(base_case, values_r, values_nc, values_nt, std_r, std_nc, std_nt, ylabel, savename)

## Fig 4 Unmodified growth rate clone 1
#title = "Unmodified growth rate clone 1"
#values_r = mean_errors_alpha_1_r
#values_nc = mean_errors_alpha_1_nc
#values_nt = mean_errors_alpha_1_nt
#std_r = std_errors_alpha_1_r
#std_nc = std_errors_alpha_1_nc
#std_nt = std_errors_alpha_1_nt
#ylabel = "Absolute error in GR"
#savename = "./article/incplots/inc-4-GR-0-clone-1.png"
#make_inc_plot(base_case, values_r, values_nc, values_nt, std_r, std_nc, std_nt, ylabel, savename)
#
## Fig 5 Unmodified growth rate clone 2
#title = "Unmodified growth rate clone 2"
#values_r = mean_errors_alpha_2_r
#values_nc = mean_errors_alpha_2_nc
#values_nt = mean_errors_alpha_2_nt
#std_r = std_errors_alpha_2_r
#std_nc = std_errors_alpha_2_nc
#std_nt = std_errors_alpha_2_nt
#ylabel = "Absolute error in GR"
#savename = "./article/incplots/inc-5-GR-0-clone-2.png"
#make_inc_plot(base_case, values_r, values_nc, values_nt, std_r, std_nc, std_nt, ylabel, savename)

