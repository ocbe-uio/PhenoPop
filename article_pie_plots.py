from make_pie_plot import *
from make_pie_plot_limitations_small_populations import *
from make_pie_plot_limitations_close_GR50s import *
from make_pie_plot_noise500 import *
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
mixturecolors = paul_tol

# Conventions: 
# Conc[0] is really 0 for all cases, but set Conc[0] to a custom value wher 0 fits nicely instead. 
# Then the tick text at that position is changed to 0

######################################################
# Fig 3A Dagim's data
######################################################

#Conc_dagim = [0, 62.5e-9, 125e-9, 250e-9, 375e-9, 500e-9, 1.25e-6, 2.5e-6, 3.75e-6, 5e-6, 7.5e-6] # Molar
Conc_dagim = np.array([30, 62.5, 125, 250, 375, 500, 1250, 2500, 3750, 5000, 7500])*1e-3 # micromolar
true_mixtures = [
[1],
[1],
[0.80, 0.20],
[2/3 ,1/3],
[0.50, 0.50]
]
true_gr50_for_plotting = [0.434916684805970, 3.98433644398336] # Sensitive 500 and Resistant 500
truegr50regionlimits = [
[Conc_dagim[4], Conc_dagim[5]],
[Conc_dagim[8], Conc_dagim[9]]
]
estimatedgr50regionlimits = [
[[Conc_dagim[4], Conc_dagim[5]]],
[[Conc_dagim[8], Conc_dagim[9]]],
[[Conc_dagim[4], Conc_dagim[5]], [Conc_dagim[8], Conc_dagim[9]]],
[[Conc_dagim[4], Conc_dagim[5]], [Conc_dagim[8], Conc_dagim[9]]],
[[Conc_dagim[5], Conc_dagim[6]], [Conc_dagim[8], Conc_dagim[9]]]
]
estimated_mixtures = [
[1],
[1],
[1-0.260201759430154, 0.260201759430154],
[1-0.356923337492292, 0.356923337492292],
[0.684759114684499, 1-0.684759114684499]
]
estimated_gr50s = [
[0.434916684805970],
[3.98433644398336],
[0.453631329448364, 4.06281089531045],
[0.408479555710222, 3.98879585161609], 
[1.13161521807390, 4.03271418945635]
]
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
true_gr50_indices = [
[0],
[1],
[0,1],
[0,1],
[0,1],
]
gr50colors_true = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
]
gr50colors_inferred = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
]
vlinerange = [0, N_cases] #[0.3, N_cases + 1 - 0.5]
title = 'True         Estimated                                          GR50 values                                 \nmixture        mixture                                                                                                 ' #'Mixture fractions                                                                   GR50 values                                                                    '#
#concticks=np.concatenate(np.array([]), np.array([Conc_dagim[1:len(Conc_dagim)]))
concticks=[Conc_dagim[0], 1e-1, 1, 1e1]  # position values
savename = "./plots/pie_plots/Dagims_pie.png"
xlabel=r"$\rm Imatinib \enspace concentration \enspace ( \mu M$)"
pctdistance=1.58
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [12,8]
#bbox_to_anchor = (1.01, N_cases*16.2/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.15, -0.3, 0.85, -0.3) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
titlepadding = 20.0
legendcolors = mixturecolors[0:2]
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc=Conc_dagim, KNOW_TRUTH=True, concticks=concticks, pctdistance=pctdistance, textprops=textprops, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, figsize=figsize, xlabel=xlabel, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, titlepadding=titlepadding, gr50colors_true=gr50colors_true, true_gr50_indices=true_gr50_indices, legendcolors=legendcolors, dagimflag=1)

######################################################
# Fig 3B Shannon's data
######################################################


Conc_shannon_erlotinib = np.array([0.0001, 0.1, 1, 10]) #micromolar
#Conc_shannon_paclitaxel = np.array([0, 1, 10, 100]) #NanoMolar
true_mixtures = [
[1],
[1],
[0.80, 0.20],
[0.50, 0.50]
]
true_gr50_for_plotting = [0.00013423, 1.1114] # cell lines 827 and 1975
truegr50regionlimits = [
[Conc_shannon_erlotinib[0], Conc_shannon_erlotinib[1]],
[Conc_shannon_erlotinib[2], Conc_shannon_erlotinib[3]]
]
estimatedgr50regionlimits = [
[[Conc_shannon_erlotinib[0], Conc_shannon_erlotinib[1]]],
[[Conc_shannon_erlotinib[2], Conc_shannon_erlotinib[3]]],
[[Conc_shannon_erlotinib[0], Conc_shannon_erlotinib[1]], [Conc_shannon_erlotinib[2], Conc_shannon_erlotinib[3]]],
[[Conc_shannon_erlotinib[0], Conc_shannon_erlotinib[1]], [Conc_shannon_erlotinib[2], Conc_shannon_erlotinib[3]]],
]
estimated_mixtures = [
[1],
[1],
[1-0.17241, 0.17241],
[1-0.47055, 0.47055]
]
estimated_gr50s = [
[0.00013423],
[1.1114],
[0.00019399, 1.0668],
[0.068014, 1.0845]
]
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
true_gr50_indices = [
[0],
[1],
[0,1],
[0,1],
]
gr50colors_true = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[0:2],
mixturecolors[0:2],
]
gr50colors_inferred = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[0:2],
mixturecolors[0:2]
]
vlinerange = [0, N_cases] #[0.3, N_cases - 0.3]
title = 'True              Estimated                                      GR50 values                               \nmixture             mixture                                                                                           ' #'Mixture fractions                                                                   GR50 values                                                                    '#
concticks=[Conc_shannon_erlotinib[0], 0.001, 0.01, 0.1, 1, 10] # position values
savename = "./plots/pie_plots/Shannons_pie.png"
xlabel=r"$\rm Erlotinib \enspace concentration \enspace ( \mu M$)"
pctdistance=1.55
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [12,7]
#bbox_to_anchor = (1.01, N_cases*15/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.2, -0.3, 0.85, -0.3) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
logaxis=True
titlepadding = 20.0
wspace = 0.3
ylim = [-1,5]
legendcolors = mixturecolors[0:2]
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc=Conc_shannon_erlotinib, KNOW_TRUTH=True, concticks=concticks, pctdistance=pctdistance, textprops=textprops, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, figsize=figsize, xlabel=xlabel, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, logaxis=logaxis, titlepadding=titlepadding, wspace=wspace, ylim=ylim, gr50colors_true=gr50colors_true, true_gr50_indices=true_gr50_indices, legendcolors=legendcolors, shannonflag=1, width_ratios=[1.1,1.1,6])


######################################################
# Fig 3D simulated data
######################################################

Noise = 10
Conc = [1e-6, 0.000005000000000, 0.000010772173450, 0.000023207944168, 0.000050000000000, 0.000107721734502, 0.000232079441681, 0.000500000000000, 0.001077217345016, 0.002320794416806, 0.005000000000000, 0.010772173450159, 0.023207944168064, 0.050000000000000, 0.107721734501594, 0.232079441680639, 0.500000000000000]
true_mixtures = [
[1],
[1],
[1],
[0.5, 0.5],
[0.7, 0.3],
[0.3, 0.7],
[1/3, 1/3, 1/3],
[0.4, 0.3, 0.3],
[0.6, 0.2, 0.2]
]
true_gr50_for_plotting = [1.22e-4, 1.16e-2, 1.12e-1]
true_gr50_conc_indices = np.zeros_like(true_gr50_for_plotting, dtype=int)
for ii in range(len(true_gr50_for_plotting)):
    elem = true_gr50_for_plotting[ii]
    try:
        true_gr50_conc_indices[ii] = next(idx for idx, value in enumerate(Conc) if value > elem)
    except:
        true_gr50_conc_indices[ii] = -1
# Currently uses last interval if GR50 greater than max_conc. Use this one if all is safe: 
#true_gr50_conc_indices = [next(idx for idx, value in enumerate(Conc) if value > elem) for elem in true_gr50_for_plotting]
aa = true_gr50_conc_indices
truegr50regionlimits = [
[Conc[aa[0]-1], Conc[aa[0]]],
[Conc[aa[1]-1], Conc[aa[1]]],
[Conc[aa[2]-1], Conc[aa[2]]]
]
estimated_mixtures = [
[1],
[1],
[1],
[0.5, 0.5],
[0.7, 0.3],
[0.3, 0.7],
[0.33, 0.33, 0.34],
[0.4, 0.3, 0.3],
[0.60, 0.19, 0.21]
]
estimated_gr50s = [
[1.30e-4],
[1.34e-2],
[1.11e-1],
[1.33e-4, 1.33e-2],
[1.23e-4, 1.35e-2],
[1.31e-4, 1.32e-2],
[1.37e-4, 9.95e-3, 9.36e-2],
[1.37e-4, 8.43e-3, 1.01e-1],
[1.42e-4, 1.21e-2, 4.64e-2]
]
estimated_gr50_conc_indices = [[0],[0],[0],[0,0],[0,0],[0,0],[0,0,0],[0,0,0],[0,0,0]]
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
[[Conc[bb[0][0]-1], Conc[bb[0][0]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]]],
[[Conc[bb[2][0]-1], Conc[bb[2][0]]]],
[[Conc[bb[3][0]-1], Conc[bb[3][0]]], [Conc[bb[3][1]-1], Conc[bb[3][1]]]],
[[Conc[bb[4][0]-1], Conc[bb[4][0]]], [Conc[bb[4][1]-1], Conc[bb[4][1]]]],
[[Conc[bb[5][0]-1], Conc[bb[5][0]]], [Conc[bb[5][1]-1], Conc[bb[5][1]]]],
[[Conc[bb[6][0]-1], Conc[bb[6][0]]], [Conc[bb[6][1]-1], Conc[bb[6][1]]], [Conc[bb[6][2]-1], Conc[bb[6][2]]]],
[[Conc[bb[7][0]-1], Conc[bb[7][0]]], [Conc[bb[7][1]-1], Conc[bb[7][1]]], [Conc[bb[7][2]-1], Conc[bb[7][2]]]],
[[Conc[bb[8][0]-1], Conc[bb[8][0]]], [Conc[bb[8][1]-1], Conc[bb[8][1]]], [Conc[bb[8][2]-1], Conc[bb[8][2]]]]
]
# Then rectify it if some of them are outside 
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
true_gr50_indices = [ # indices in the list of true populations
[0],
[1],
[2],
[0,1],
[0,1],
[0,1],
[0,1,2],
[0,1,2],
[0,1,2],
]
gr50colors_true = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
gr50colors_inferred = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
vlinerange = [0, N_cases]
title = 'True         Estimated                                       GR50 values                                    \nmixture        mixture                                                                                                 ' #'Mixture fractions                                                                   GR50 values                                                                    '#
savename = "./plots/pie_plots/simulated-noise-"+str(Noise)+".png"
pctdistance=1.59
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [12,14]
concticks=[Conc[0], 10**-5,10**-4,10**-3,10**-2,10**-1,1] # position values
#bbox_to_anchor = (1.01, N_cases*18/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.05, -0.3, 0.7, -0.3) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
legendcolors = mixturecolors[0:3]
titlepadding=15.0
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc, KNOW_TRUTH=True, concticks=concticks, pctdistance=pctdistance, textprops=textprops, figsize=figsize, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, gr50colors_true=gr50colors_true, true_gr50_indices=true_gr50_indices, legendcolors=legendcolors, titlepadding=titlepadding, skewflag=1)



Noise = 50
Conc = [1e-6, 0.000005000000000, 0.000010772173450, 0.000023207944168, 0.000050000000000, 0.000107721734502, 0.000232079441681, 0.000500000000000, 0.001077217345016, 0.002320794416806, 0.005000000000000, 0.010772173450159, 0.023207944168064, 0.050000000000000, 0.107721734501594, 0.232079441680639, 0.500000000000000]
true_mixtures = [
[1],
[1],
[1],
[0.5, 0.5],
[0.7, 0.3],
[0.3, 0.7],
[1/3, 1/3, 1/3],
[0.4, 0.3, 0.3],
[0.6, 0.2, 0.2]
]
true_gr50_for_plotting = [1.22e-4, 1.16e-2, 1.12e-1]
true_gr50_conc_indices = np.zeros_like(true_gr50_for_plotting, dtype=int)
for ii in range(len(true_gr50_for_plotting)):
    elem = true_gr50_for_plotting[ii]
    try:
        true_gr50_conc_indices[ii] = next(idx for idx, value in enumerate(Conc) if value > elem)
    except:
        true_gr50_conc_indices[ii] = -1
# Currently uses last interval if GR50 greater than max_conc. Use this one if all is safe: 
#true_gr50_conc_indices = [next(idx for idx, value in enumerate(Conc) if value > elem) for elem in true_gr50_for_plotting]
aa = true_gr50_conc_indices
truegr50regionlimits = [
[Conc[aa[0]-1], Conc[aa[0]]],
[Conc[aa[1]-1], Conc[aa[1]]],
[Conc[aa[2]-1], Conc[aa[2]]]
]
estimated_mixtures = [
[1],
[1],
[1],
[0.49, 0.51],
[0.69, 0.31],
[0.29, 0.71],
[0.33, 0.35, 0.32],
[0.39, 0.32, 0.29],
[0.59, 0.20, 0.21]
]
estimated_gr50s = [
[1.38e-4],
[1.33e-2],
[1.41e-1],
[1.44e-4, 1.26e-2],
[1.41e-4, 1.10e-2],
[1.54e-4, 1.33e-2],
[1.51e-4, 1.08e-2, 1.13e-1],
[1.52e-4, 1.32e-2, 8.92e-2],
[7.13e-5, 5.04e-3, 8.73e-2]
]
estimated_gr50_conc_indices = [[0],[0],[0],[0,0],[0,0],[0,0],[0,0,0],[0,0,0],[0,0,0]]
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
[[Conc[bb[0][0]-1], Conc[bb[0][0]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]]],
[[Conc[bb[2][0]-1], Conc[bb[2][0]]]],
[[Conc[bb[3][0]-1], Conc[bb[3][0]]], [Conc[bb[3][1]-1], Conc[bb[3][1]]]],
[[Conc[bb[4][0]-1], Conc[bb[4][0]]], [Conc[bb[4][1]-1], Conc[bb[4][1]]]],
[[Conc[bb[5][0]-1], Conc[bb[5][0]]], [Conc[bb[5][1]-1], Conc[bb[5][1]]]],
[[Conc[bb[6][0]-1], Conc[bb[6][0]]], [Conc[bb[6][1]-1], Conc[bb[6][1]]], [Conc[bb[6][2]-1], Conc[bb[6][2]]]],
[[Conc[bb[7][0]-1], Conc[bb[7][0]]], [Conc[bb[7][1]-1], Conc[bb[7][1]]], [Conc[bb[7][2]-1], Conc[bb[7][2]]]],
[[Conc[bb[8][0]-1], Conc[bb[8][0]]], [Conc[bb[8][1]-1], Conc[bb[8][1]]], [Conc[bb[8][2]-1], Conc[bb[8][2]]]]
]
# Then rectify it if some of them are outside 
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
true_gr50_indices = [
[0],
[1],
[2],
[0,1],
[0,1],
[0,1],
[0,1,2],
[0,1,2],
[0,1,2],
]
gr50colors_true = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
gr50colors_inferred = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
vlinerange = [0, N_cases]
title = 'True         Estimated                                       GR50 values                                    \nmixture        mixture                                                                                                 ' #'Mixture fractions                                                                   GR50 values                                                                    '#
savename = "./plots/pie_plots/simulated-noise-"+str(Noise)+".png"
pctdistance=1.59
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [12,14]
concticks=[Conc[0], 10**-5,10**-4,10**-3,10**-2,10**-1,1] # position values
#bbox_to_anchor = (1.01, N_cases*18/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.05, -0.3, 0.7, -0.3) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
legendcolors = mixturecolors[0:3]
titlepadding=15.0
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc, KNOW_TRUTH=True, concticks=concticks, pctdistance=pctdistance, textprops=textprops, figsize=figsize, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, gr50colors_true=gr50colors_true, true_gr50_indices=true_gr50_indices, legendcolors=legendcolors, titlepadding=titlepadding, skewflag=1)


Noise = 100
Conc = [1e-6, 0.000005000000000, 0.000010772173450, 0.000023207944168, 0.000050000000000, 0.000107721734502, 0.000232079441681, 0.000500000000000, 0.001077217345016, 0.002320794416806, 0.005000000000000, 0.010772173450159, 0.023207944168064, 0.050000000000000, 0.107721734501594, 0.232079441680639, 0.500000000000000]
true_mixtures = [
[1],
[1],
[1],
[0.5, 0.5],
[0.7, 0.3],
[0.3, 0.7],
[1/3, 1/3, 1/3],
[0.4, 0.3, 0.3],
[0.6, 0.2, 0.2]
]
true_gr50_for_plotting = [1.22e-4, 1.16e-2, 1.12e-1]
true_gr50_conc_indices = np.zeros_like(true_gr50_for_plotting, dtype=int)
for ii in range(len(true_gr50_for_plotting)):
    elem = true_gr50_for_plotting[ii]
    try:
        true_gr50_conc_indices[ii] = next(idx for idx, value in enumerate(Conc) if value > elem)
    except:
        true_gr50_conc_indices[ii] = -1
# Currently uses last interval if GR50 greater than max_conc. Use this one if all is safe: 
#true_gr50_conc_indices = [next(idx for idx, value in enumerate(Conc) if value > elem) for elem in true_gr50_for_plotting]
aa = true_gr50_conc_indices
truegr50regionlimits = [
[Conc[aa[0]-1], Conc[aa[0]]],
[Conc[aa[1]-1], Conc[aa[1]]],
[Conc[aa[2]-1], Conc[aa[2]]]
]
estimated_mixtures = [
[1],
[1],
[1],
[0.49, 0.51],
[0.69, 0.31],
[0.29, 0.71],
[0.33, 0.36, 0.31],
[0.38, 0.37, 0.25],
[0.59, 0.23, 0.18],
]
estimated_gr50s = [
[1.50e-4],
[1.31e-2],
[1.41e-1],
[1.62e-4, 7.11e-3],
[1.56e-4, 6.21e-3],
[1.72e-4, 1.30e-2],
[1.98e-4, 1.16e-2, 6.40e-2],
[1.32e-4, 1.58e-2, 1.52e-1],
[1.58e-4, 1.31e-2, 4.71e-2],
]
estimated_gr50_conc_indices = [[0],[0],[0],[0,0],[0,0],[0,0],[0,0,0],[0,0,0],[0,0,0]]
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
[[Conc[bb[0][0]-1], Conc[bb[0][0]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]]],
[[Conc[bb[2][0]-1], Conc[bb[2][0]]]],
[[Conc[bb[3][0]-1], Conc[bb[3][0]]], [Conc[bb[3][1]-1], Conc[bb[3][1]]]],
[[Conc[bb[4][0]-1], Conc[bb[4][0]]], [Conc[bb[4][1]-1], Conc[bb[4][1]]]],
[[Conc[bb[5][0]-1], Conc[bb[5][0]]], [Conc[bb[5][1]-1], Conc[bb[5][1]]]],
[[Conc[bb[6][0]-1], Conc[bb[6][0]]], [Conc[bb[6][1]-1], Conc[bb[6][1]]], [Conc[bb[6][2]-1], Conc[bb[6][2]]]],
[[Conc[bb[7][0]-1], Conc[bb[7][0]]], [Conc[bb[7][1]-1], Conc[bb[7][1]]], [Conc[bb[7][2]-1], Conc[bb[7][2]]]],
[[Conc[bb[8][0]-1], Conc[bb[8][0]]], [Conc[bb[8][1]-1], Conc[bb[8][1]]], [Conc[bb[8][2]-1], Conc[bb[8][2]]]]
]
# Then rectify it if some of them are outside 
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
true_gr50_indices = [
[0],
[1],
[2],
[0,1],
[0,1],
[0,1],
[0,1,2],
[0,1,2],
[0,1,2],
]
gr50colors_true = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
gr50colors_inferred = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
vlinerange = [0, N_cases]
title = 'True         Estimated                                       GR50 values                                    \nmixture        mixture                                                                                                 ' #'Mixture fractions                                                                   GR50 values                                                                    '#
savename = "./plots/pie_plots/simulated-noise-"+str(Noise)+".png"
pctdistance=1.59
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [12,14]
concticks=[Conc[0], 10**-5,10**-4,10**-3,10**-2,10**-1,1] # position values
#bbox_to_anchor = (1.01, N_cases*18/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.05, -0.3, 0.7, -0.3) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
legendcolors = mixturecolors[0:3]
titlepadding=15.0
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc, KNOW_TRUTH=True, concticks=concticks, pctdistance=pctdistance, textprops=textprops, figsize=figsize, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, gr50colors_true=gr50colors_true, true_gr50_indices=true_gr50_indices, legendcolors=legendcolors, titlepadding=titlepadding, skewflag=1)


Noise = 150
Conc = [1e-6, 0.000005000000000, 0.000010772173450, 0.000023207944168, 0.000050000000000, 0.000107721734502, 0.000232079441681, 0.000500000000000, 0.001077217345016, 0.002320794416806, 0.005000000000000, 0.010772173450159, 0.023207944168064, 0.050000000000000, 0.107721734501594, 0.232079441680639, 0.500000000000000]
true_mixtures = [
[1],
[1],
[1],
[0.5, 0.5],
[0.7, 0.3],
[0.3, 0.7],
[1/3, 1/3, 1/3],
[0.4, 0.3, 0.3],
[0.6, 0.2, 0.2]
]
true_gr50_for_plotting = [1.22e-4, 1.16e-2, 1.12e-1]
true_gr50_conc_indices = np.zeros_like(true_gr50_for_plotting, dtype=int)
for ii in range(len(true_gr50_for_plotting)):
    elem = true_gr50_for_plotting[ii]
    try:
        true_gr50_conc_indices[ii] = next(idx for idx, value in enumerate(Conc) if value > elem)
    except:
        true_gr50_conc_indices[ii] = -1
# Currently uses last interval if GR50 greater than max_conc. Use this one if all is safe: 
#true_gr50_conc_indices = [next(idx for idx, value in enumerate(Conc) if value > elem) for elem in true_gr50_for_plotting]
aa = true_gr50_conc_indices
truegr50regionlimits = [
[Conc[aa[0]-1], Conc[aa[0]]],
[Conc[aa[1]-1], Conc[aa[1]]],
[Conc[aa[2]-1], Conc[aa[2]]]
]
estimated_mixtures = [
[1],
[1],
[1],
[0.49, 0.51],
[0.69, 0.31],
[0.3,0.7],
[0.33, 0.30, 0.37], #[0.327, 0.297, 0.376],
[0.40, 0.33, 0.27],
[0.60, 0.22, 0.18], #[0.594, 0.223, 0.183],
]
estimated_gr50s = [
[8.08e-5],
[1.28e-2],
[1.41e-1],
[1.75e-4, 6.40e-3],
[1.69e-4, 5.39e-3],
[7.92e-5, 1.28e-2],
[1.22e-4, 1.10e-2, 1.38e-1],
[1.73e-4, 7.99e-3, 3.40e-2],
[1.76e-4, 6.23e-3, 4.15e-2],
]
estimated_gr50_conc_indices = [[0],[0],[0],[0,0],[0,0],[0,0],[0,0,0],[0,0,0],[0,0,0]]
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
[[Conc[bb[0][0]-1], Conc[bb[0][0]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]]],
[[Conc[bb[2][0]-1], Conc[bb[2][0]]]],
[[Conc[bb[3][0]-1], Conc[bb[3][0]]], [Conc[bb[3][1]-1], Conc[bb[3][1]]]],
[[Conc[bb[4][0]-1], Conc[bb[4][0]]], [Conc[bb[4][1]-1], Conc[bb[4][1]]]],
[[Conc[bb[5][0]-1], Conc[bb[5][0]]], [Conc[bb[5][1]-1], Conc[bb[5][1]]]],
[[Conc[bb[6][0]-1], Conc[bb[6][0]]], [Conc[bb[6][1]-1], Conc[bb[6][1]]], [Conc[bb[6][2]-1], Conc[bb[6][2]]]],
[[Conc[bb[7][0]-1], Conc[bb[7][0]]], [Conc[bb[7][1]-1], Conc[bb[7][1]]], [Conc[bb[7][2]-1], Conc[bb[7][2]]]],
[[Conc[bb[8][0]-1], Conc[bb[8][0]]], [Conc[bb[8][1]-1], Conc[bb[8][1]]], [Conc[bb[8][2]-1], Conc[bb[8][2]]]]
]
# Then rectify it if some of them are outside 
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
true_gr50_indices = [
[0],
[1],
[2],
[0,1],
[0,1],
[0,1],
[0,1,2],
[0,1,2],
[0,1,2],
]
gr50colors_true = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
gr50colors_inferred = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
vlinerange = [0, N_cases]
title = 'True         Estimated                                       GR50 values                                    \nmixture        mixture                                                                                                 ' #'Mixture fractions                                                                   GR50 values                                                                    '#
savename = "./plots/pie_plots/simulated-noise-"+str(Noise)+".png"
pctdistance=1.59
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [12,14]
concticks=[Conc[0], 10**-5,10**-4,10**-3,10**-2,10**-1,1] # position values
#bbox_to_anchor = (1.01, N_cases*18/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.05, -0.3, 0.7, -0.3) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
legendcolors = mixturecolors[0:3]
titlepadding=15.0
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc, KNOW_TRUTH=True, concticks=concticks, pctdistance=pctdistance, textprops=textprops, figsize=figsize, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, gr50colors_true=gr50colors_true, true_gr50_indices=true_gr50_indices, legendcolors=legendcolors, titlepadding=titlepadding, skewflag=1)



Noise = 200
Conc = [1e-6, 0.000005000000000, 0.000010772173450, 0.000023207944168, 0.000050000000000, 0.000107721734502, 0.000232079441681, 0.000500000000000, 0.001077217345016, 0.002320794416806, 0.005000000000000, 0.010772173450159, 0.023207944168064, 0.050000000000000, 0.107721734501594, 0.232079441680639, 0.500000000000000]
true_mixtures = [
[1],
[1],
[1],
[0.5, 0.5],
[0.7, 0.3],
[0.3, 0.7],
[1/3, 1/3, 1/3],
[0.4, 0.3, 0.3],
[0.6, 0.2, 0.2]
]
true_gr50_for_plotting = [1.22e-4, 1.16e-2, 1.12e-1]
true_gr50_conc_indices = np.zeros_like(true_gr50_for_plotting, dtype=int)
for ii in range(len(true_gr50_for_plotting)):
    elem = true_gr50_for_plotting[ii]
    try:
        true_gr50_conc_indices[ii] = next(idx for idx, value in enumerate(Conc) if value > elem)
    except:
        true_gr50_conc_indices[ii] = -1
# Currently uses last interval if GR50 greater than max_conc. Use this one if all is safe: 
#true_gr50_conc_indices = [next(idx for idx, value in enumerate(Conc) if value > elem) for elem in true_gr50_for_plotting]
aa = true_gr50_conc_indices
truegr50regionlimits = [
[Conc[aa[0]-1], Conc[aa[0]]],
[Conc[aa[1]-1], Conc[aa[1]]],
[Conc[aa[2]-1], Conc[aa[2]]]
]
estimated_mixtures = [
[1],
[1],
[1],
[0.5, 0.5],
[0.7, 0.3],
[0.3, 0.7],
[0.33, 0.35, 0.32],
[0.39, 0.37, 0.24],
[0.59, 0.28, 0.13],
]
estimated_gr50s = [
[7.75e-5],
[6.97e-3],
[1.42e-1],
[1.34e-4, 5.75e-3],
[1.72e-4, 4.80e-3],
[8.71e-5, 6.44e-3],
[1.04e-4, 4.36e-3, 1.06e-1],
[1.39e-4, 1.93e-2, 1.25e-1],
[1.75e-4, 1.44e-2, 2.61e-1],
]
estimated_gr50_conc_indices = [[0],[0],[0],[0,0],[0,0],[0,0],[0,0,0],[0,0,0],[0,0,0]]
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
[[Conc[bb[0][0]-1], Conc[bb[0][0]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]]],
[[Conc[bb[2][0]-1], Conc[bb[2][0]]]],
[[Conc[bb[3][0]-1], Conc[bb[3][0]]], [Conc[bb[3][1]-1], Conc[bb[3][1]]]],
[[Conc[bb[4][0]-1], Conc[bb[4][0]]], [Conc[bb[4][1]-1], Conc[bb[4][1]]]],
[[Conc[bb[5][0]-1], Conc[bb[5][0]]], [Conc[bb[5][1]-1], Conc[bb[5][1]]]],
[[Conc[bb[6][0]-1], Conc[bb[6][0]]], [Conc[bb[6][1]-1], Conc[bb[6][1]]], [Conc[bb[6][2]-1], Conc[bb[6][2]]]],
[[Conc[bb[7][0]-1], Conc[bb[7][0]]], [Conc[bb[7][1]-1], Conc[bb[7][1]]], [Conc[bb[7][2]-1], Conc[bb[7][2]]]],
[[Conc[bb[8][0]-1], Conc[bb[8][0]]], [Conc[bb[8][1]-1], Conc[bb[8][1]]], [Conc[bb[8][2]-1], Conc[bb[8][2]]]]
]
# Then rectify it if some of them are outside 
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
true_gr50_indices = [
[0],
[1],
[2],
[0,1],
[0,1],
[0,1],
[0,1,2],
[0,1,2],
[0,1,2],
]
gr50colors_true = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
gr50colors_inferred = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
vlinerange = [0, N_cases]
title = 'True         Estimated                                       GR50 values                                    \nmixture        mixture                                                                                                 ' #'Mixture fractions                                                                   GR50 values                                                                    '#
savename = "./plots/pie_plots/simulated-noise-"+str(Noise)+".png"
pctdistance=1.59
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [12,14]
concticks=[Conc[0], 10**-5,10**-4,10**-3,10**-2,10**-1,1] # position values
#bbox_to_anchor = (1.01, N_cases*18/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.05, -0.3, 0.7, -0.3) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
legendcolors = mixturecolors[0:3]
titlepadding=15.0
make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc, KNOW_TRUTH=True, concticks=concticks, pctdistance=pctdistance, textprops=textprops, figsize=figsize, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, gr50colors_true=gr50colors_true, true_gr50_indices=true_gr50_indices, legendcolors=legendcolors, titlepadding=titlepadding, skewflag=1)


Noise = 500
Conc = [1e-6, 0.000005000000000, 0.000010772173450, 0.000023207944168, 0.000050000000000, 0.000107721734502, 0.000232079441681, 0.000500000000000, 0.001077217345016, 0.002320794416806, 0.005000000000000, 0.010772173450159, 0.023207944168064, 0.050000000000000, 0.107721734501594, 0.232079441680639, 0.500000000000000]
true_mixtures = [
[1],
[1],
[1],
[0.5, 0.5],
[0.7, 0.3],
[0.3, 0.7],
[1/3, 1/3, 1/3],
[0.4, 0.3, 0.3],
[0.6, 0.2, 0.2]
]
true_gr50_for_plotting = [1.22e-4, 1.16e-2, 1.12e-1]
true_gr50_conc_indices = np.zeros_like(true_gr50_for_plotting, dtype=int)
for ii in range(len(true_gr50_for_plotting)):
    elem = true_gr50_for_plotting[ii]
    try:
        true_gr50_conc_indices[ii] = next(idx for idx, value in enumerate(Conc) if value > elem)
    except:
        true_gr50_conc_indices[ii] = -1
# Currently uses last interval if GR50 greater than max_conc. Use this one if all is safe: 
#true_gr50_conc_indices = [next(idx for idx, value in enumerate(Conc) if value > elem) for elem in true_gr50_for_plotting]
aa = true_gr50_conc_indices
truegr50regionlimits = [
[Conc[aa[0]-1], Conc[aa[0]]],
[Conc[aa[1]-1], Conc[aa[1]]],
[Conc[aa[2]-1], Conc[aa[2]]]
]
estimated_mixtures = [
[1],
[1],
[1],
[0.59, 0.41],
[0.75, 0.25],
[0.42, 0.58],
[0.5, 0.5],
[0.56, 0.44],
[0.7, 0.3],
#[0.41, 0.44, 0.15],
#[0.44, 0.44, 0.12],
#[0.62, 0.38, 0.00],
]
estimated_gr50s = [
[5.43e-5],
[5.38e-3],
[2.01e-1],
[8.36e-5, 4.49e-3],
[8.43e-5, 4.04e-3],
[6.79e-5, 4.92e-3],
[7.23e-5, 2.03e-2],
[7.94e-5, 1.86e-2],
[8.49e-5, 1.43e-2],
#[8.24e-5, 6.67e-3, 5.09e-2],
#[4.15e-5, 6.87e-3, 1.02e-1],
#[1.05e-4, 6.12e-3, 3.06e-1],
]
estimated_gr50_conc_indices = [[0],[0],[0],[0,0],[0,0],[0,0],[0,0,0],[0,0,0],[0,0,0]]
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
[[Conc[bb[0][0]-1], Conc[bb[0][0]]]],
[[Conc[bb[1][0]-1], Conc[bb[1][0]]]],
[[Conc[bb[2][0]-1], Conc[bb[2][0]]]],
[[Conc[bb[3][0]-1], Conc[bb[3][0]]], [Conc[bb[3][1]-1], Conc[bb[3][1]]]],
[[Conc[bb[4][0]-1], Conc[bb[4][0]]], [Conc[bb[4][1]-1], Conc[bb[4][1]]]],
[[Conc[bb[5][0]-1], Conc[bb[5][0]]], [Conc[bb[5][1]-1], Conc[bb[5][1]]]],
[[Conc[bb[6][0]-1], Conc[bb[6][0]]], [Conc[bb[6][1]-1], Conc[bb[6][1]]], [Conc[bb[6][2]-1], Conc[bb[6][2]]]],
[[Conc[bb[7][0]-1], Conc[bb[7][0]]], [Conc[bb[7][1]-1], Conc[bb[7][1]]], [Conc[bb[7][2]-1], Conc[bb[7][2]]]],
[[Conc[bb[8][0]-1], Conc[bb[8][0]]], [Conc[bb[8][1]-1], Conc[bb[8][1]]], [Conc[bb[8][2]-1], Conc[bb[8][2]]]]
]
# Then rectify it if some of them are outside 
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
true_gr50_indices = [
[0],
[1],
[2],
[0,1],
[0,1],
[0,1],
[0,1,2],
[0,1,2],
[0,1,2],
]
gr50colors_true = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
gr50colors_inferred = [
mixturecolors[0:1],
mixturecolors[1:2],
mixturecolors[2:3],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
vlinerange = [0, N_cases]
title = 'True         Estimated                                       GR50 values                                    \nmixture        mixture                                                                                                 ' #'Mixture fractions                                                                   GR50 values                                                                    '#
savename = "./plots/pie_plots/simulated-noise-"+str(Noise)+".png"
pctdistance=1.59
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [12,14]
concticks=[Conc[0], 10**-5,10**-4,10**-3,10**-2,10**-1,1] # position values
#bbox_to_anchor = (1.01, N_cases*18/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.05, -0.3, 0.7, -0.3) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
legendcolors = mixturecolors[0:3]

ylim = [-1,N_cases+1]
titlepadding=15.0
make_pie_plot_noise500(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc, KNOW_TRUTH=True, concticks=concticks, pctdistance=pctdistance, textprops=textprops, figsize=figsize, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, gr50colors_true=gr50colors_true, true_gr50_indices=true_gr50_indices, legendcolors=legendcolors, ylim=ylim, titlepadding=titlepadding, skewflag=1)

######################################################
# Limitations: small populations
######################################################

Noise = 50
Conc = [1e-6, 0.000005000000000, 0.000010772173450, 0.000023207944168, 0.000050000000000, 0.000107721734502, 0.000232079441681, 0.000500000000000, 0.001077217345016, 0.002320794416806, 0.005000000000000, 0.010772173450159, 0.023207944168064, 0.050000000000000, 0.107721734501594, 0.232079441680639, 0.500000000000000]
true_mixtures = [
#[0.70, 0.30], # is included in other plot
#[0.80, 0.20], # 20 is big
[0.90, 0.10],
[0.95, 0.05],
[0.97, 0.03],
[0.99, 0.01],
[0.50, 0.45, 0.05],
[0.50, 0.47, 0.03],
[0.50, 0.49, 0.01],
[0.80, 0.10, 0.10],
[0.90, 0.05, 0.05],
#[0.95, 0.02, 0.03] # 2 and 3 percent labels overlap each other when drawn 
]
true_gr50_for_plotting = [1.22e-4, 1.16e-2, 1.12e-1]
true_gr50_conc_indices = np.zeros_like(true_gr50_for_plotting, dtype=int)
for ii in range(len(true_gr50_for_plotting)):
    elem = true_gr50_for_plotting[ii]
    try:
        true_gr50_conc_indices[ii] = next(idx for idx, value in enumerate(Conc) if value > elem)
    except:
        true_gr50_conc_indices[ii] = -1
# Currently uses last interval if GR50 greater than max_conc. Use this one if all is safe: 
#true_gr50_conc_indices = [next(idx for idx, value in enumerate(Conc) if value > elem) for elem in true_gr50_for_plotting]
aa = true_gr50_conc_indices
truegr50regionlimits = [
[Conc[aa[0]-1], Conc[aa[0]]],
[Conc[aa[1]-1], Conc[aa[1]]],
[Conc[aa[2]-1], Conc[aa[2]]]
]
estimated_mixtures = [
#[0.69, 0.31],
#[0.80, 0.20],
[0.90, 0.10],
[0.95, 0.05], # elbow ok, difference ca 200
[0.97, 0.03], # elbow ok, difference ca 100
[0.99, 0.01], # elbow ok, difference ca 27 
[0.49, 0.46, 0.05], # elbow ok, difference ca 60 
[0.50, 0.48, 0.02], # elbow ok, 3038 vs 3009
[0.49, 0.51], # [0.49, 0.50, 0.01] # Elbow not ok, 3013 vs 3009
[0.80, 0.11, 0.09], # elbow ok, difference ca 70
[0.91, 0.09], #[0.91, 0.01, 0.08], # ELBOW difference BELOW 18, detect 2 instead of 3
#[0.95, 0.05] #[0.95, 0.02, 0.03] # ELBOW difference BELOW 18, detect 2 instead of 3
]
estimated_gr50s = [
#[1.41e-4, 1.10e-2],
#[1.42e-4, 6.38e-3],
[1.42e-4, 5.58e-3],
[1.45e-4, 4.32e-3],
[1.44e-4, 3.00e-3],
[1.41e-4, 5.06e-3],
[1.41e-4, 9.57e-3, 1.66e-1], # Two large and 1 small 
[1.50e-4, 1.33e-2, 1.25e-1], # Two large and 1 small 
[1.46e-4, 1.32e-2], # [1.46e-4, 1.36e-2, 5.52e-2] # Elbow not ok, 3013 vs 3009
[1.39e-4, 1.24e-2, 5.48e-2],
[1.44e-4, 1.59e-1],  #[1.46e-4, 9.94e-3, 7.71e-2],
#[1.37e-4, 9.88e-2] #[1.46e-4, 1.16e-2, 2.54e-1]
]
estimated_gr50_conc_indices = [[0,0],[0,0],[0,0],[0,0],[0,0,0],[0,0,0],[0,0],[0,0,0],[0,0],[0,0]] #[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0,0],[0,0,0],[0,0],[0,0,0],[0,0],[0,0]]
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
[[Conc[bb[4][0]-1], Conc[bb[4][0]]], [Conc[bb[4][1]-1], Conc[bb[4][1]]], [Conc[bb[4][2]-1], Conc[bb[4][2]]]],
[[Conc[bb[5][0]-1], Conc[bb[5][0]]], [Conc[bb[5][1]-1], Conc[bb[5][1]]], [Conc[bb[5][2]-1], Conc[bb[5][2]]]],
[[Conc[bb[6][0]-1], Conc[bb[6][0]]], [Conc[bb[6][1]-1], Conc[bb[6][1]]]], #, [Conc[bb[8][2]-1], Conc[bb[8][2]]]]
[[Conc[bb[7][0]-1], Conc[bb[7][0]]], [Conc[bb[7][1]-1], Conc[bb[7][1]]], [Conc[bb[7][2]-1], Conc[bb[7][2]]]],
[[Conc[bb[8][0]-1], Conc[bb[8][0]]], [Conc[bb[8][1]-1], Conc[bb[8][1]]]], #, [Conc[bb[10][2]-1], Conc[bb[10][2]]]]
#[[Conc[bb[11][0]-1], Conc[bb[11][0]]], [Conc[bb[11][1]-1], Conc[bb[11][1]]]], #, [Conc[bb[11][2]-1], Conc[bb[11][2]]]]
]
# Then rectify it if some of them are outside 
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
true_gr50_indices = [
[0,1],
[0,1],
[0,1],
[0,1],
[0,1,2],
[0,1,2],
[0,1,2],
[0,1,2],
[0,1,2],
]
gr50colors_true = [
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
]
gr50colors_inferred = [
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
vlinerange = [0, N_cases]
title = 'True         Estimated                                       GR50 values                                    \nmixture        mixture                                                                                                 ' #'Mixture fractions                                                                   GR50 values                                                                    '#
savename = "./plots/pie_plots/limitations_small_populations.png"
pctdistance=1.57
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [12,14]
concticks=[Conc[0], 10**-5,10**-4,10**-3,10**-2,10**-1,1] # position values
#bbox_to_anchor = (1.01, N_cases*18/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.05, -0.3, 0.7, -0.3) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
legendcolors = mixturecolors[0:3]
titlepadding=15.0
make_pie_plot_limitations_small_populations(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc, KNOW_TRUTH=True, concticks=concticks, pctdistance=pctdistance, textprops=textprops, figsize=figsize, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, gr50colors_true=gr50colors_true, true_gr50_indices=true_gr50_indices, legendcolors=legendcolors, titlepadding=titlepadding, skewflag=1)


######################################################
# Limitations: Close GR50 values
######################################################
#mixturecolors = [mixturecolors[0], mixturecolors[1], mixturecolors[2], mixturecolors[2]] #mixturecolors = limitationcolors
Conc = [1e-6, 0.000005000000000, 0.000010772173450, 0.000023207944168, 0.000050000000000, 0.000107721734502, 0.000232079441681, 0.000500000000000, 0.001077217345016, 0.002320794416806, 0.005000000000000, 0.010772173450159, 0.023207944168064, 0.050000000000000, 0.107721734501594, 0.232079441680639, 0.500000000000000]
true_mixtures = [
[0.5, 0.5],
[0.5, 0.5],
[0.5, 0.5],
[0.90, 0.10],
[0.90, 0.10],
[0.90, 0.10],
[0.95, 0.05],
[0.95, 0.05],
[0.95, 0.05],
#[0.99, 0.01],
#[0.99, 0.01],
#[0.99, 0.01],
]
true_gr50_for_plotting = [0.000122221, 0.00025099014, 0.00052100, 0.001222211745] # cell lines 1 2 3 4
# All true GR50s for all 7 populations: [0.000122221, 0.00025099014, 0.00052100, 0.001222211745, 0.0025099012, 0.00521000326, 0.0122219] 
true_gr50_conc_indices = np.zeros_like(true_gr50_for_plotting, dtype=int)
for ii in range(len(true_gr50_for_plotting)):
    elem = true_gr50_for_plotting[ii]
    try:
        true_gr50_conc_indices[ii] = next(idx for idx, value in enumerate(Conc) if value > elem)
    except:
        true_gr50_conc_indices[ii] = -1
# Currently uses last interval if GR50 greater than max_conc. Use this one if all is safe: 
#true_gr50_conc_indices = [next(idx for idx, value in enumerate(Conc) if value > elem) for elem in true_gr50_for_plotting]
aa = true_gr50_conc_indices
truegr50regionlimits = [
[Conc[aa[0]-1], Conc[aa[0]]],
[Conc[aa[1]-1], Conc[aa[1]]],
[Conc[aa[2]-1], Conc[aa[2]]],
[Conc[aa[3]-1], Conc[aa[3]]],
]
estimated_mixtures = [
# Elbow power: cases 1-9 are clear. Case 10 is also ok: The difference is 19.59
# Case 11 and 12 are also ok, the differences are over 20
[0.49, 0.51], #3
[0.50, 0.50], #2
[0.51, 0.49], #1
[0.90, 0.10], #6
[0.92, 0.08], #5
[0.97, 0.03], #4
[0.95, 0.05], #9
[0.97, 0.03], #8
[0.98, 0.02], #7
#[0.97, 0.03], #12
#[0.95, 0.05], #11
#[0.98, 0.02], #10
]
estimated_gr50s = [
[1.44e-4, 7.62e-4], #3
[1.49e-4, 3.77e-4], #2
[1.61e-4, 1.84e-4], #1
[1.43e-4, 6.66e-4], #6
[1.49e-4, 2.41e-4], #5
[1.56e-4, 1.00e-4], #4
[1.43e-4, 6.19e-4], #9
[1.49e-4, 1.95e-4], #8
[1.50e-4, 7.68e-5], #7
#[1.38e-4, 1.18e-5], #12
#[1.36e-4, 2.06e-5], #11
#[1.40e-4, 1.15e-5], #10
]
estimated_gr50_conc_indices = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
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
[[Conc[bb[5][0]-1], Conc[bb[5][0]]], [Conc[bb[5][1]-1], Conc[bb[5][1]]]],
[[Conc[bb[6][0]-1], Conc[bb[6][0]]], [Conc[bb[6][1]-1], Conc[bb[6][1]]]],
[[Conc[bb[7][0]-1], Conc[bb[7][0]]], [Conc[bb[7][1]-1], Conc[bb[7][1]]]],
[[Conc[bb[8][0]-1], Conc[bb[8][0]]], [Conc[bb[8][1]-1], Conc[bb[8][1]]]],
]
# Then rectify it if some of them are outside 
N_cases = len(estimated_gr50s)
N_inferred_populations = [len(elem) for elem in estimated_gr50s]
true_gr50_indices = [
[0,3], #3
[0,2], #2
[0,1], #1
[0,3], #3
[0,2], #2
[0,1], #1
[0,3], #3
[0,2], #2
[0,1], #1
]
gr50colors_true = [
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
]
gr50colors_inferred = [
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:2],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3],
mixturecolors[0:3]
]
vlinerange = [0, N_cases]
title = 'True         Estimated                                       GR50 values                                    \nmixture        mixture                                                                                                 ' #'Mixture fractions                                                                   GR50 values                                                                    '#
savename = "./plots/pie_plots/limitations_close_GR50s.png"
pctdistance=1.59
textprops = {'fontsize':16, 'weight' : 'bold'}
figsize = [12,14]
concticks=[Conc[0], 10**-5,10**-4,10**-3,10**-2,10**-1,1] # position values
#bbox_to_anchor = (1.01, N_cases*12/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
bbox_to_anchor = (-0.1, -0.3, 0.72, -0.3) # positions legend with space allotment in (x, y, width, height)
ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
ylim = [-1,N_cases+1]
legendcolors = mixturecolors[0:3]
legendcolors_inferred = mixturecolors[0:2]
titlepadding=10.0
startangles = [0,0,270]
make_pie_plot_limitations_close_GR50s(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc, KNOW_TRUTH=True, concticks=concticks, pctdistance=pctdistance, textprops=textprops, figsize=figsize, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, gr50colors_true=gr50colors_true, true_gr50_indices=true_gr50_indices, ylim=ylim, legendcolors=legendcolors, titlepadding=titlepadding, legendcolors_inferred=legendcolors_inferred, limitationflag=1)

# Old limitation plot with anecdotes at noise levels 200 and 1000: 

##mixturecolors = [mixturecolors[0], mixturecolors[1], mixturecolors[2], mixturecolors[2]] #mixturecolors = limitationcolors
#Conc = [1e-6, 0.000005000000000, 0.000010772173450, 0.000023207944168, 0.000050000000000, 0.000107721734502, 0.000232079441681, 0.000500000000000, 0.001077217345016, 0.002320794416806, 0.005000000000000, 0.010772173450159, 0.023207944168064, 0.050000000000000, 0.107721734501594, 0.232079441680639, 0.500000000000000]
#true_mixtures = [
#[1],
#[0.5, 0.5],
##[1/3, 1/3, 1/3],
##[0.95, 0.05],
#[0.05, 0.90, 0.05],
#]
#true_gr50_for_plotting = [0.000122221, 0.00052100, 0.001222211745, 0.0025099012, 0.0122219] # [0.000122221, 0.00025099014, 0.00052100, 0.001222211745, 0.0025099012, 0.00521000326, 0.0122219] 
#true_gr50_conc_indices = np.zeros_like(true_gr50_for_plotting, dtype=int)
#for ii in range(len(true_gr50_for_plotting)):
#    elem = true_gr50_for_plotting[ii]
#    try:
#        true_gr50_conc_indices[ii] = next(idx for idx, value in enumerate(Conc) if value > elem)
#    except:
#        true_gr50_conc_indices[ii] = -1
## Currently uses last interval if GR50 greater than max_conc. Use this one if all is safe: 
##true_gr50_conc_indices = [next(idx for idx, value in enumerate(Conc) if value > elem) for elem in true_gr50_for_plotting]
#aa = true_gr50_conc_indices
#truegr50regionlimits = [
#[Conc[aa[0]-1], Conc[aa[0]]],
#[Conc[aa[1]-1], Conc[aa[1]]],
#[Conc[aa[2]-1], Conc[aa[2]]],
#[Conc[aa[3]-1], Conc[aa[3]]],
#[Conc[aa[4]-1], Conc[aa[4]]],
#]
#estimated_mixtures = [
#[0.52, 0.48], #[1],
#[0.90, 0.10],
##[0.38, 0.50, 0.12], #[0.59, 0.41],
##[0.89, 0.11],
#[0.88, 0.12],
#]
#estimated_gr50s = [
#[2.96e-4, np.inf], #[2.15e-4],
#[4.48e-4, np.inf],#[1.19e-5, 4.48e-4],
##[2.92e-4, 8.90e-4, np.inf], #[2.92e-4, 8.90e-4, 1.92e-2], # [4.32e-4, 1.08e-3],
##[4.12e-4, np.inf], #[1.17e-5, 4.12e-4],
#[5.82e-4, np.inf], #[1.21e-5, 5.82e-4],
#]
#estimated_gr50_conc_indices = [[0,0],[0,0],[0,0],[0,0]] #[[0,0],[0,0],[0,0,0],[0,0],[0,0]]
#for ii in range(len(estimated_gr50s)):
#    for jj in range(len(estimated_gr50s[ii])):
#        elem = estimated_gr50s[ii][jj]
#        try:
#            estimated_gr50_conc_indices[ii][jj] = next(idx for idx, value in enumerate(Conc) if value > elem)
#        except:
#            estimated_gr50_conc_indices[ii][jj] = -1
## Currently uses last interval if GR50 greater than max_conc.
#bb = estimated_gr50_conc_indices
#estimatedgr50regionlimits = [
#[[Conc[bb[0][0]-1], Conc[bb[0][0]]]],
#[[Conc[bb[1][0]-1], Conc[bb[1][0]]], [Conc[bb[1][1]-1], Conc[bb[1][1]]]],
#[[Conc[bb[2][0]-1], Conc[bb[2][0]]], [Conc[bb[2][1]-1], Conc[bb[2][1]]]], #, [Conc[bb[2][2]-1], Conc[bb[2][2]]]],
##[[Conc[bb[3][0]-1], Conc[bb[3][0]]], [Conc[bb[3][1]-1], Conc[bb[3][1]]]],
##[[Conc[bb[4][0]-1], Conc[bb[4][0]]], [Conc[bb[4][1]-1], Conc[bb[4][1]]]],
#]
## Then rectify it if some of them are outside 
#N_cases = len(estimated_gr50s)
#N_inferred_populations = [len(elem) for elem in estimated_gr50s]
#true_gr50_indices = [
#[1],
#[1,2],
##[1,2,3],
##[1,4],
#[0,1,4],
#]
#gr50colors_true = [
#mixturecolors[0:2],
#mixturecolors[0:2],
#[mixturecolors[0], mixturecolors[1], mixturecolors[2], mixturecolors[2]],
#]
#gr50colors_inferred = [
#mixturecolors[0:2],
#mixturecolors[0:2],
#mixturecolors[0:2],
#]
##gr50colors_true = [
##mixturecolors[1:], #mixturecolors[0:], 
##mixturecolors[1:], #mixturecolors[0:], 
###mixturecolors[1:], #mixturecolors[0:], 
###[mixturecolors[1],mixturecolors[4]], #mixturecolors[0:], 
##[mixturecolors[0],mixturecolors[1],mixturecolors[4]], #mixturecolors[0:], 
##]
##gr50colors_inferred = [
##[mixturecolors[1], mixturecolors[4]],
##[mixturecolors[1], mixturecolors[4]],
###[mixturecolors[1], mixturecolors[2], mixturecolors[4]],
###[mixturecolors[1], mixturecolors[4]],
##[mixturecolors[1], mixturecolors[4]],
##]
#vlinerange = [0, N_cases]
#title = 'True         Estimated                                       GR50 values                                    \nmixture        mixture                                                                                                 ' #'Mixture fractions                                                                   GR50 values                                                                    '#
#savename = "./plots/pie_plots/limitations.png"
#pctdistance=1.59
#textprops = {'fontsize':16, 'weight' : 'bold'}
#figsize = [12,6]
#concticks=[Conc[0], 10**-5,10**-4,10**-3,10**-2,10**-1,1] # position values
##bbox_to_anchor = (1.01, N_cases*12/20, 0.7, 1) # positions legend with space allotment in (x, y, width, height)
#bbox_to_anchor = (-0.1, -0.3, 0.72, -0.3) # positions legend with space allotment in (x, y, width, height)
#ax_set_position=[] #[0.6,0.1,0.36,0.9] # Controls figure size for subplot in (x, y, width, height)
#ylim = [-1,N_cases+1]
#legendcolors = mixturecolors[0:3]
#legendcolors_inferred = mixturecolors[0:2]
#titlepadding=10.0
#startangles = [0,0,270]
#make_pie_plot(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc, KNOW_TRUTH=True, concticks=concticks, pctdistance=pctdistance, textprops=textprops, figsize=figsize, truegr50regionlimits=truegr50regionlimits, estimatedgr50regionlimits=estimatedgr50regionlimits, bbox_to_anchor=bbox_to_anchor, ax_set_position=ax_set_position, gr50colors_true=gr50colors_true, true_gr50_indices=true_gr50_indices, ylim=ylim, legendcolors=legendcolors, titlepadding=titlepadding, legendcolors_inferred=legendcolors_inferred, limitationflag=1)

