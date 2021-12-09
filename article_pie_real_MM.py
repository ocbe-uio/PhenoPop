from make_pie_plot import *
#from make_patient_pie_plot import *
from moodboard import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# Pie charts, where the slices will be ordered and plotted counter-clockwise:
#plt.pie(x, explode=None, labels=None, colors=None, autopct=None0 pctdistance=0.6, 
#    shadow=False, labeldistance=1.1, startangle=0,, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'} radius=1, counterclock=True, 
#    wedgeprops=None, textprops=None, center=(0, 0), frame=False, rotatelabels=False, 
#    *, normalize=None, data=None)[source]
#ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
mixturecolors = paul_tol_darker #cartocolors_prism #cartocolors_3_5
mixturecolors_light = paul_tol_lighter #cartocolors_prism #cartocolors_3_5

patient_names = ["MM1420", "MM19520", "MM3620", "MM210819", "MM130120"] # Patients 1,2,3,4,5 with ids 0,1,2,3,4
colordict = {
"MM1420" : [mixturecolors[3], mixturecolors_light[3]], #desaturate(mixturecolors[3], 0.7)], # np.array([1, 0.6])*mixturecolors[4],
"MM130120" : [mixturecolors[4], mixturecolors_light[4]], #desaturate(mixturecolors[5], 0.7)], # np.array([1, 0.6])*mixturecolors[6],
"MM210819" : [mixturecolors[5], mixturecolors_light[5]], #desaturate(mixturecolors[4], 0.7)], # np.array([1, 0.6])*mixturecolors[5],
"MM19520" : [mixturecolors[6], mixturecolors_light[6]], #desaturate(mixturecolors[6], 0.7)], # np.array([1, 0.6])*mixturecolors[7],
"MM3620" : [mixturecolors[7], mixturecolors_light[7]], #desaturate(mixturecolors[7], 0.7)], # np.array([1, 0.6])*mixturecolors[8],
}
labels_all_patients = [
["Patient MM1420 clone 1", "Patient MM1420 clone 2"], 
["Patient MM19520 clone 1", "Patient MM19520 clone 2"], 
["Patient MM3620 clone 1", "Patient MM3620 clone 2"],
["Patient MM210819 clone 1", "Patient MM210819 clone 2"], 
["Patient MM130120 clone 1", "Patient MM130120 clone 2"], 
]
# Conventions: 
# Conc[0] is really 0 for all cases, but set Conc[0] to a custom value wher 0 fits nicely instead. 
# Concentrations are originally in nanomolar
# Then the tick text at that position is changed to 0
# MM130120 alias MM720

## define concentration (from Elise's code: data_plots.py)                       My understanding: (which is then wrong)
#concentration_drug1 = [0, 1, 5, 10, 50, 100, 500, 1000, 10000] #ixa/veneto      ixa
#concentration_drug2 = [0, 0.1, 1, 5, 10, 100, 500, 1000, 10000] #thalido/melf   thali
#concentration_drug3 = [0, 0.1, 0.5, 1, 5, 10, 50, 100, 1000] #seli              seli / melflufen
#concentration_drug4 = [0, 0.1, 1, 10, 50, 100, 500, 1000, 10000] #dexa          dexa / venetoclax
Conc_venetoclax =    [0, 0.1, 1, 10, 50, 100, 500, 1000, 10000] # checked
Conc_melflufen =     [0, 0.1, 0.5, 1, 5, 10, 50, 100, 1000]     # checked
Conc_ixazomib =      [0, 1, 5, 10, 50, 100, 500, 1000, 10000]   # checked
Conc_thalidomide =   [0, 0.1, 1, 5, 10, 100, 500, 1000, 10000]  # checked
Conc_selinexor =     [0, 0.1, 0.5, 1, 5, 10, 50, 100, 1000]     # checked
Conc_dexamethasone = [0, 0.1, 1, 10, 50, 100, 500, 1000, 10000] # checked

######################################################
# Fig 3C real-MM data
######################################################
Conc_1_210819 = [0, 1, 5, 10, 50, 100, 500, 1000, 10000]
Conc_2_210819 = [0, 0.1, 1, 5, 10, 100, 500, 1000, 10000]
Conc_3_210819 = [0, 0.1, 0.5, 1, 5, 10, 50, 100, 1000]
Conc_4_210819 = [0, 0.1, 1, 10, 50, 100, 500, 1000, 10000]

Conc_1_130120 = [0, 1, 5, 10, 50, 100, 500, 1000, 10000]
Conc_2_130120 = [0, 0.1, 1, 5, 10, 100, 500, 1000, 10000]
Conc_3_130120 = [0, 0.1, 0.5, 1, 5, 10, 50, 100, 1000]
Conc_4_130120 = [0, 0.1, 1, 10, 50, 100, 500, 1000, 10000]

Conc_1_19520 =  [0, 0.1, 1, 10, 50, 100, 500, 1000, 10000]
Conc_2_19520 =  [0, 0.1, 0.5, 1, 5, 10, 50, 100, 1000]
Conc_3_19520 =  [0, 0.1, 0.5, 1, 5, 10, 50, 100, 1000]
Conc_4_19520 =  [0, 0.1, 1, 10, 50, 100, 500, 1000, 10000]

Conc_1_3620 =   [0, 0.1, 1, 10, 50, 100, 500, 1000, 10000]
Conc_2_3620 =   [0, 0.1, 0.5, 1, 5, 10, 50, 100, 1000]
Conc_3_3620 =   [0, 0.1, 0.5, 1, 5, 10, 50, 100, 1000]
Conc_4_3620 =   [0, 0.1, 1, 10, 50, 100, 500, 1000, 10000]

Conc_1_1420 =   [0, 0.1, 1, 10, 50, 100, 500, 1000, 10000]
Conc_2_1420 =   [0, 0.1, 0.5, 1, 5, 10, 50, 100, 1000]
Conc_3_1420 =   [0, 0.1, 0.5, 1, 5, 10, 50, 100, 1000]
Conc_4_1420 =   [0, 0.1, 1, 10, 50, 100, 500, 1000, 10000]


# Patient MM1420
true_mixtures = [0]
#true_gr50s = [0]
true_gr50_for_plotting = [0]
#max_N_populations = len(true_gr50_for_plotting)
###
# Drugs 1,2,3,4
N_inferred_populations = [2,2,2,2]
estimated_mixtures_MM1420 = [
[0.42, 0.58],
[0.62, 0.38],
[0.28, 0.72],
[0.13, 0.87]
]
estimated_gr50s_MM1420 = [
[499.88, np.Inf], #566.25],
[971.41, np.Inf], #39.17],
[9.24, 246.10], #[246.10, np.inf], #9.24],
[98.59, np.Inf], #1015.61]
]
N_cases = len(estimated_gr50s_MM1420)
gr50colors_inferred = [
mixturecolors,
mixturecolors,
mixturecolors,
mixturecolors
]
vlinerange = [0.7, 4.3]
title = 'Patient MM1420\nMixture fractions                                               GR50 values                       '
savename = "./article/Fig 3 GR50 pie charts/MM1420.png"
legendcolors = mixturecolors[0:2]
#make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures_MM1420, estimated_gr50s_MM1420, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, KNOW_TRUTH=False, legendcolors=legendcolors)


# Patient MM3620
true_mixtures = [0]
#true_gr50s = [0]
true_gr50_for_plotting = [0]
#max_N_populations = len(true_gr50_for_plotting)
###
# Drugs 1,2,3,4
N_inferred_populations = [2,1,1,1]
estimated_mixtures_MM3620 = [
[0.77, 0.23],
[1],
[1],
[1]
]
estimated_gr50s_MM3620 = [
[106.56, np.Inf], #[106.56, 164.14], #[106.56, np.Inf], # 164.14],
[660.89],
[980.36],
[4.20]
]
N_cases = len(estimated_gr50s_MM3620)
gr50colors_inferred = [
mixturecolors,
mixturecolors,
mixturecolors,
mixturecolors
]
vlinerange = [0.7, 4.3]
title = 'Patient MM3620\nMixture fractions                                               GR50 values                       '
savename = "./article/Fig 3 GR50 pie charts/MM3620.png"
legendcolors = mixturecolors[0:2]
#make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures_MM3620, estimated_gr50s_MM3620, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, KNOW_TRUTH=False, legendcolors=legendcolors)

# Patient MM19520
true_mixtures = [0]
#true_gr50s = [0]
true_gr50_for_plotting = [0]
#max_N_populations = len(true_gr50_for_plotting)
###
# Drugs 1,2,3,4
N_inferred_populations = [2,1,1,1]
estimated_mixtures_MM19520 = [
[0.5, 0.5],
[1],
[1],
[1]
]
estimated_gr50s_MM19520 = [
[1177.21, np.Inf], #44.25],
[624.59],
[666.75],
[Conc_dexamethasone[1]*0.1] #3.28e-6]
]
N_cases = len(estimated_gr50s_MM19520)
gr50colors_inferred = [
mixturecolors,
mixturecolors,
mixturecolors,
mixturecolors
]
vlinerange = [0.7, 4.3]
title = 'Patient MM19520\nMixture fractions                                               GR50 values                       '
savename = "./article/Fig 3 GR50 pie charts/MM19520.png"
legendcolors = mixturecolors[0:2]
#make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures_MM19520, estimated_gr50s_MM19520, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, KNOW_TRUTH=False, legendcolors=legendcolors)
# Patient MM130120
true_mixtures = [0]
#true_gr50s = [0]
true_gr50_for_plotting = [0]
#max_N_populations = len(true_gr50_for_plotting)
###
# Drugs 1,2,3,4
N_inferred_populations = [1,1,1,1]
estimated_mixtures_MM130120 = [
[1],
[1],
[1],
[1]
]
estimated_gr50s_MM130120 = [
[52.13],
[9645.55],
[147.96],
[992.59]
]
N_cases = len(estimated_gr50s_MM130120)
gr50colors_inferred = [
mixturecolors,
mixturecolors,
mixturecolors,
mixturecolors
]
vlinerange = [0.7, 4.3]
title = 'Patient MM130120\nMixture fractions                                               GR50 values                       '
savename = "./article/Fig 3 GR50 pie charts/MM130120.png"
legendcolors = mixturecolors[0:2]
#make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures_MM130120, estimated_gr50s_MM130120, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, KNOW_TRUTH=False, legendcolors=legendcolors)

# Patient MM210819

true_mixtures = [0]
#true_gr50s = [0]
true_gr50_for_plotting = [0]
#max_N_populations = len(true_gr50_for_plotting)
###
# Drugs 1,2,3,4
N_inferred_populations = [1,1,1,1]
estimated_mixtures_MM210819 = [
[1],
[1],
[1],
[1]
]
estimated_gr50s_MM210819 = [
[56.29],
[726.85],
[129.35],
[8371.74]
]
N_cases = len(estimated_gr50s_MM210819)
gr50colors_inferred = [
mixturecolors,
mixturecolors,
mixturecolors,
mixturecolors
]
vlinerange = [0.7, 4.3]
title = 'Patient MM210819\nMixture fractions                                               GR50 values                       '
savename = "./article/Fig 3 GR50 pie charts/MM210819.png"
legendcolors = mixturecolors[0:2]
#make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures_MM210819, estimated_gr50s_MM210819, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, KNOW_TRUTH=False, legendcolors=legendcolors)

##########################################
# One for each drug
##########################################
# Venetoclax
Conc = [0.01] + Conc_venetoclax[1:]
# Patients with this drug according to final_results.ods
# 1420 nr 0
# 19520 nr 0
# 3620 nr 0
labels_patients = [labels_all_patients[idx] for idx in [0,1,2]] # [0,3,4]]
true_mixtures = []
true_gr50_for_plotting = []
truegr50regionlimits = []
estimated_mixtures = [
estimated_mixtures_MM1420[0],
estimated_mixtures_MM19520[0],
estimated_mixtures_MM3620[0],
]
estimated_gr50s = [
estimated_gr50s_MM1420[0],
estimated_gr50s_MM19520[0],
estimated_gr50s_MM3620[0],
]
estimated_gr50_conc_indices = [[0,0],[0,0],[0,0]]
for ii in range(len(estimated_gr50s)):
    for jj in range(len(estimated_gr50s[ii])):
        elem = estimated_gr50s[ii][jj]
        try:
            estimated_gr50_conc_indices[ii][jj] = next(idx for idx, value in enumerate(Conc) if value > elem)
        except:
            estimated_gr50_conc_indices[ii][jj] = -1
# Currently uses last interval if GR50 greater than max_conc.
bb = estimated_gr50_conc_indices
estimatedgr50regionlimits = [
[[Conc[bb[0][0]-1], Conc[bb[0][0]]], [Conc[bb[0][1]-1], Conc[bb[0][1]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]], [Conc[bb[1][1]-1], Conc[bb[1][1]]]],
[[Conc[bb[2][0]-1], Conc[bb[2][0]]], [Conc[bb[2][1]-1], Conc[bb[2][1]]]],
]
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
gr50colors_inferred = [colordict["MM1420"][0:N_inferred_populations[0]], colordict["MM19520"][0:N_inferred_populations[1]], colordict["MM3620"][0:N_inferred_populations[2]]]
vlinerange = [0.3, N_cases - 0.3]
title = 'Mixture fractions                      GR50 values                                                         '#
concticks=[Conc[0], 0.1, 1, 10**1, 10**2, 10**3, 10**4] # position values
savename = "./article/Fig 3 GR50 pie charts/MM-Venetoclax.png"
xlabel=r"$\rm Venetoclax \enspace concentration \enspace ( \mu M$)"
pctdistance=1.7
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [10,6]
#bbox_to_anchor = (1.01, N_cases*12/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.16, -0.25, 1, -0.25) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
logaxis=True
titlepadding = 0.0
wspace = 0.3
ylim = [-1,N_cases+1]
legendcolors = mixturecolors[0:2]
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc=Conc, KNOW_TRUTH=False, concticks=concticks, pctdistance=pctdistance, textprops=textprops, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, figsize=figsize, xlabel=xlabel, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, logaxis=logaxis, titlepadding=titlepadding, wspace=wspace, labels_patients=labels_patients, ylim=ylim, legendcolors=legendcolors)

# Melflufen
Conc = [0.01] + Conc_melflufen[1:]
# Patients with this drug according to final_results.ods
# 1420 nr 1
# 19520 nr 1
# 3620 nr 1
labels_patients = [labels_all_patients[idx] for idx in [0,1,2]] # [0,3,4]]
true_mixtures = []
true_gr50_for_plotting = []
truegr50regionlimits = []
estimated_mixtures = [
estimated_mixtures_MM1420[1],
estimated_mixtures_MM19520[1],
estimated_mixtures_MM3620[1],
]
estimated_gr50s = [
estimated_gr50s_MM1420[1],
estimated_gr50s_MM19520[1],
estimated_gr50s_MM3620[1],
]
estimated_gr50_conc_indices = [[0,0],[0,0],[0,0]]
for ii in range(len(estimated_gr50s)):
    for jj in range(len(estimated_gr50s[ii])):
        elem = estimated_gr50s[ii][jj]
        try:
            estimated_gr50_conc_indices[ii][jj] = next(idx for idx, value in enumerate(Conc) if value > elem)
        except:
            estimated_gr50_conc_indices[ii][jj] = -1
# Currently uses last interval if GR50 greater than max_conc.
bb = estimated_gr50_conc_indices
estimatedgr50regionlimits = [
[[Conc[bb[0][0]-1], Conc[bb[0][0]]], [Conc[bb[0][1]-1], Conc[bb[0][1]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]], [Conc[bb[1][1]-1], Conc[bb[1][1]]]],
[[Conc[bb[2][0]-1], Conc[bb[2][0]]], [Conc[bb[2][1]-1], Conc[bb[2][1]]]],
]
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
gr50colors_inferred = [colordict["MM1420"][0:N_inferred_populations[0]], colordict["MM19520"][0:N_inferred_populations[1]], colordict["MM3620"][0:N_inferred_populations[2]]]
vlinerange = [0.8, N_cases - 0.8]
title = '   Mixture fractions                      GR50 values                                                         '#
concticks=[Conc[0], 0.1, 1, 10**1, 10**2, 10**3] # position values
savename = "./article/Fig 3 GR50 pie charts/MM-Melflufen.png"
xlabel=r"$\rm Melflufen \enspace concentration \enspace ( \mu M$)"
pctdistance=1.7
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [10,6]
#bbox_to_anchor = (1.01, N_cases*12/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.16, -0.25, 1, -0.25) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
logaxis=True
titlepadding = 10.0
wspace = 0.3
ylim = [0,N_cases]
legendcolors = mixturecolors[0:2]
arrow_multiplier = 2.3
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc=Conc, KNOW_TRUTH=False, concticks=concticks, pctdistance=pctdistance, textprops=textprops, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, figsize=figsize, xlabel=xlabel, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, logaxis=logaxis, titlepadding=titlepadding, wspace=wspace, labels_patients=labels_patients, ylim=ylim, legendcolors=legendcolors, arrow_multiplier=arrow_multiplier)

# Ixazomib
Conc = [0.1] + Conc_ixazomib[1:]
# Patients with this drug according to final_results.ods
# MM210819 nr 0
# MM130120 alias MM720 nr 0
labels_patients = [labels_all_patients[idx] for idx in [3,4]] # [1,2]]
true_mixtures = []
true_gr50_for_plotting = []
truegr50regionlimits = []
estimated_mixtures = [
estimated_mixtures_MM210819[0],
estimated_mixtures_MM130120[0],
]
estimated_gr50s = [
estimated_gr50s_MM210819[0],
estimated_gr50s_MM130120[0],
]
estimated_gr50_conc_indices = [[0,0],[0,0]]
for ii in range(len(estimated_gr50s)):
    for jj in range(len(estimated_gr50s[ii])):
        elem = estimated_gr50s[ii][jj]
        try:
            estimated_gr50_conc_indices[ii][jj] = next(idx for idx, value in enumerate(Conc) if value > elem)
        except:
            estimated_gr50_conc_indices[ii][jj] = -1
# Currently uses last interval if GR50 greater than max_conc.
bb = estimated_gr50_conc_indices
estimatedgr50regionlimits = [
[[Conc[bb[0][0]-1], Conc[bb[0][0]]], [Conc[bb[0][1]-1], Conc[bb[0][1]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]], [Conc[bb[1][1]-1], Conc[bb[1][1]]]],
]
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
gr50colors_inferred = [colordict["MM210819"][0:N_inferred_populations[0]], colordict["MM130120"][0:N_inferred_populations[1]]]
vlinerange = [0.3, N_cases - 0.3]
title = '   Mixture fractions                      GR50 values                                                         '#
concticks=[Conc[0], 1, 10**1, 10**2, 10**3, 10**4] # position values
savename = "./article/Fig 3 GR50 pie charts/MM-Ixazomib.png"
xlabel=r"$\rm Ixazomib \enspace concentration \enspace ( \mu M$)"
pctdistance=1.7
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [10,4]
#bbox_to_anchor = (1.01, N_cases*7/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.16, -0.25, 1, -0.25) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
logaxis=True
titlepadding = 10.0
wspace = 0.3
hspace = 0.0
ylim = [-0.25,N_cases+0.25]
legendcolors = mixturecolors[0:2]
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc=Conc, KNOW_TRUTH=False, concticks=concticks, pctdistance=pctdistance, textprops=textprops, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, figsize=figsize, xlabel=xlabel, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, logaxis=logaxis, titlepadding=titlepadding, wspace=wspace, labels_patients=labels_patients, hspace=hspace, ylim=ylim, legendcolors=legendcolors)

# Thalidomide
Conc = [0.01] + Conc_thalidomide[1:]
# Patients with this drug according to final_results.ods
# MM210819 nr 1
# MM130120 alias MM720 nr 1
labels_patients = [labels_all_patients[idx] for idx in [3,4]] # [1,2]]
true_mixtures = []
true_gr50_for_plotting = []
truegr50regionlimits = []
estimated_mixtures = [
estimated_mixtures_MM210819[1],
estimated_mixtures_MM130120[1],
]
estimated_gr50s = [
estimated_gr50s_MM210819[1],
estimated_gr50s_MM130120[1],
]
estimated_gr50_conc_indices = [[0,0],[0,0]]
for ii in range(len(estimated_gr50s)):
    for jj in range(len(estimated_gr50s[ii])):
        elem = estimated_gr50s[ii][jj]
        try:
            estimated_gr50_conc_indices[ii][jj] = next(idx for idx, value in enumerate(Conc) if value > elem)
        except:
            estimated_gr50_conc_indices[ii][jj] = -1
# Currently uses last interval if GR50 greater than max_conc.
bb = estimated_gr50_conc_indices
estimatedgr50regionlimits = [
[[Conc[bb[0][0]-1], Conc[bb[0][0]]], [Conc[bb[0][1]-1], Conc[bb[0][1]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]], [Conc[bb[1][1]-1], Conc[bb[1][1]]]],
]
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
gr50colors_inferred = [colordict["MM210819"][0:N_inferred_populations[0]], colordict["MM130120"][0:N_inferred_populations[1]]]
vlinerange = [0.3, N_cases - 0.3]
title = '   Mixture fractions                      GR50 values                                                         '#
concticks=[Conc[0], 0.1, 1, 10**1, 10**2, 10**3, 10**4] # position values
savename = "./article/Fig 3 GR50 pie charts/MM-Thalidomide.png"
xlabel=r"$\rm Thalidomide \enspace concentration \enspace ( \mu M$)"
pctdistance=1.7
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [10,4]
#bbox_to_anchor = (1.01, N_cases*7/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.16, -0.25, 1, -0.25) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
logaxis=True
titlepadding = 10.0
wspace = 0.3
ylim = [-0.25,N_cases+0.25]
legendcolors = mixturecolors[0:2]
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc=Conc, KNOW_TRUTH=False, concticks=concticks, pctdistance=pctdistance, textprops=textprops, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, figsize=figsize, xlabel=xlabel, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, logaxis=logaxis, titlepadding=titlepadding, wspace=wspace, labels_patients=labels_patients, ylim=ylim, legendcolors=legendcolors)

# Selinexor
Conc = [0.01] + Conc_selinexor[1:]
# Patients with this drug according to final_results.ods
# Everyone index nr 2
labels_patients = labels_all_patients
#labels_patients = [labels_all_patients[idx] for idx in [0,3,4,1,2]]
true_mixtures = []
true_gr50_for_plotting = []
truegr50regionlimits = []
estimated_mixtures = [
estimated_mixtures_MM1420[2],
estimated_mixtures_MM19520[2],
estimated_mixtures_MM3620[2],
estimated_mixtures_MM210819[2],
estimated_mixtures_MM130120[2],
]
estimated_gr50s = [
estimated_gr50s_MM1420[2],
estimated_gr50s_MM19520[2],
estimated_gr50s_MM3620[2],
estimated_gr50s_MM210819[2],
estimated_gr50s_MM130120[2],
]
estimated_gr50_conc_indices = [[0,0],[0,0],[0,0],[0,0],[0,0]]
for ii in range(len(estimated_gr50s)):
    for jj in range(len(estimated_gr50s[ii])):
        elem = estimated_gr50s[ii][jj]
        try:
            estimated_gr50_conc_indices[ii][jj] = next(idx for idx, value in enumerate(Conc) if value > elem)
        except:
            estimated_gr50_conc_indices[ii][jj] = -1
# Currently uses last interval if GR50 greater than max_conc.
bb = estimated_gr50_conc_indices
estimatedgr50regionlimits = [
[[Conc[bb[0][0]-1], Conc[bb[0][0]]], [Conc[bb[0][1]-1], Conc[bb[0][1]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]], [Conc[bb[1][1]-1], Conc[bb[1][1]]]],
[[Conc[bb[2][0]-1], Conc[bb[2][0]]], [Conc[bb[2][1]-1], Conc[bb[2][1]]]],
[[Conc[bb[3][0]-1], Conc[bb[3][0]]], [Conc[bb[3][1]-1], Conc[bb[3][1]]]],
[[Conc[bb[4][0]-1], Conc[bb[4][0]]], [Conc[bb[4][1]-1], Conc[bb[4][1]]]],
]
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
gr50colors_inferred = [colordict[patname][0:N_inferred_populations[idx]] for idx, patname in enumerate(patient_names)]
vlinerange = [0.3, N_cases - 0.3]
title = '   Mixture fractions                      GR50 values                                                         '#
concticks=[Conc[0], 0.1, 1, 10**1, 10**2, 10**3] # position values
savename = "./article/Fig 3 GR50 pie charts/MM-Selinexor.png"
xlabel=r"$\rm Selinexor \enspace concentration \enspace ( \mu M$)"
pctdistance=1.6
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [10,7]
#bbox_to_anchor = (1.01, N_cases*21/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.16, -0.45, 1, -0.45) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
logaxis=True
titlepadding = 40.0
wspace = 0.3
hspace = 0.3
legendcolors = mixturecolors[0:2]
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc=Conc, KNOW_TRUTH=False, concticks=concticks, pctdistance=pctdistance, textprops=textprops, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, figsize=figsize, xlabel=xlabel, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, logaxis=logaxis, titlepadding=titlepadding, wspace=wspace, hspace=hspace, labels_patients=labels_patients, legendcolors=legendcolors)

# Dexamethasone
Conc = [0.01] + Conc_dexamethasone[1:]
# Patients with this drug according to final_results.ods
# Everyone index nr 3
labels_patients = labels_all_patients
#labels_patients = [labels_all_patients[idx] for idx in [0,3,4,1,2]]
true_mixtures = []
true_gr50_for_plotting = []
truegr50regionlimits = []
estimated_mixtures = [
estimated_mixtures_MM1420[3],
estimated_mixtures_MM19520[3],
estimated_mixtures_MM3620[3],
estimated_mixtures_MM210819[3],
estimated_mixtures_MM130120[3],
]
estimated_gr50s = [
estimated_gr50s_MM1420[3],
estimated_gr50s_MM19520[3],
estimated_gr50s_MM3620[3],
estimated_gr50s_MM210819[3],
estimated_gr50s_MM130120[3],
]
estimated_gr50_conc_indices = [[0,0],[0,0],[0,0],[0,0],[0,0]]
for ii in range(len(estimated_gr50s)):
    for jj in range(len(estimated_gr50s[ii])):
        elem = estimated_gr50s[ii][jj]
        try:
            estimated_gr50_conc_indices[ii][jj] = next(idx for idx, value in enumerate(Conc) if value > elem)
        except:
            estimated_gr50_conc_indices[ii][jj] = -1
# Currently uses last interval if GR50 greater than max_conc.
bb = estimated_gr50_conc_indices
estimatedgr50regionlimits = [
[[Conc[bb[0][0]-1], Conc[bb[0][0]]], [Conc[bb[0][1]-1], Conc[bb[0][1]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]], [Conc[bb[1][1]-1], Conc[bb[1][1]]]],
[[Conc[bb[2][0]-1], Conc[bb[2][0]]], [Conc[bb[2][1]-1], Conc[bb[2][1]]]],
[[Conc[bb[3][0]-1], Conc[bb[3][0]]], [Conc[bb[3][1]-1], Conc[bb[3][1]]]],
[[Conc[bb[4][0]-1], Conc[bb[4][0]]], [Conc[bb[4][1]-1], Conc[bb[4][1]]]],
]
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
gr50colors_inferred = [colordict[patname][0:N_inferred_populations[idx]] for idx, patname in enumerate(patient_names)]
vlinerange = [0.3, N_cases - 0.3]
title = '   Mixture fractions                      GR50 values                                                         '#
concticks=[Conc[0], 0.1, 1, 10**1, 10**2, 10**3, 10**4] # position values
savename = "./article/Fig 3 GR50 pie charts/MM-Dexamethasone.png"
xlabel=r"$\rm Dexamethasone \enspace concentration \enspace ( \mu M$)"
pctdistance=1.5
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [10,7]
#bbox_to_anchor = (1.01, N_cases*21/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.16, -0.45, 1, -0.45) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
logaxis=True
titlepadding = 40.0
wspace = 0.3
hspace = 0.3
legendcolors = mixturecolors[0:2]
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc=Conc, KNOW_TRUTH=False, concticks=concticks, pctdistance=pctdistance, textprops=textprops, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, figsize=figsize, xlabel=xlabel, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, logaxis=logaxis, titlepadding=titlepadding, wspace=wspace, hspace=hspace, labels_patients=labels_patients, legendcolors=legendcolors)

