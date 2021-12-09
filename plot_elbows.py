import matplotlib.pyplot as plt
import numpy as np
import scipy.io


######################################## Simulated data ########################################
def get_savestr(N_true_populations, caseno, Noise, num_optim):
    if N_true_populations == 1:
        savestr = '-Noise-' + str(Noise) + '-num_optim-' + str(num_optim)
        MixParams = 1
    elif N_true_populations == 2:
        if caseno == 1:
            MixParams = [0.5]
        elif caseno == 2:
            MixParams = [0.7]
        elif caseno == 3:
            MixParams = [0.3]
        savestr = '-Noise-' + str(Noise) + '-mix-' + str(MixParams[0]) +  '-num_optim-' + str(num_optim)
    elif N_true_populations == 3:
        if caseno == 1:
            MixParams = [0.33333, 0.33333]
        elif caseno == 2:
            MixParams = [0.4, 0.3]
        elif caseno == 3:
            MixParams = [0.6, 0.2]
        savestr = '-Noise-' + str(Noise) + '-mix-' + str(MixParams[0]) + '-' + str(MixParams[1]) + '-num_optim-' + str(num_optim)
    return savestr, MixParams

Noise = 150
num_optim = 1000
max_N_populations = 4
saved_gr50s_at_runtime = True

max_neg_ll = 0
min_neg_ll = np.Inf
max_difference = 0

########## Find max values ##########
for N_true_populations in [1,2,3]:
    for caseno in [1,2,3]:
        savestr, MixParams = get_savestr(N_true_populations, caseno, Noise, num_optim)
        filename = './article/CLUSTERPLOTS/' + str(N_true_populations) + '/case' + str(caseno) + '/negative_loglikelihood_values' + savestr + '.mat'
        mat_contents = scipy.io.loadmat(filename)
        negative_loglikelihood_values = mat_contents['negative_loglikelihood_values'][0]
        if max(negative_loglikelihood_values) > max_neg_ll:
            max_neg_ll = max(negative_loglikelihood_values)
        if min(negative_loglikelihood_values) < min_neg_ll:
            min_neg_ll = min(negative_loglikelihood_values)
        if max(negative_loglikelihood_values) - min(negative_loglikelihood_values) > max_difference:
            max_difference = max(negative_loglikelihood_values) - min(negative_loglikelihood_values)

########## Plot ##########
fig, axes = plt.subplots(nrows = 3, ncols=3, figsize=(10,10))
plt.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0.2, hspace=0.2)
for N_true_populations in [1,2,3]:
    for caseno in [1,2,3]:
        savestr, MixParams = get_savestr(N_true_populations, caseno, Noise, num_optim)
        filename = './article/CLUSTERPLOTS/' + str(N_true_populations) + '/case' + str(caseno) + '/negative_loglikelihood_values' + savestr + '.mat'
        mat_contents = scipy.io.loadmat(filename)
        negative_loglikelihood_values = mat_contents['negative_loglikelihood_values'][0]
        if N_true_populations == 1:
            labell = "Truth: " + str(N_true_populations) + " clone"
        elif N_true_populations == 2:
            fullMixParams = MixParams + [1-sum(MixParams)]
            labell = "Truth: " + str(N_true_populations) + " clones\nMixtures {:.2f}".format(fullMixParams[0]) + ", {:.2f}".format(fullMixParams[1])
        elif N_true_populations == 3:
            fullMixParams = MixParams + [1-sum(MixParams)]
            labell = "Truth: " + str(N_true_populations) + " clones\nMixtures {:.2f}".format(fullMixParams[0]) + ", {:.2f}".format(fullMixParams[1]) + ", {:.2f}".format(fullMixParams[2])
        axes[N_true_populations-1, caseno-1].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
        if N_true_populations == 3:
            axes[N_true_populations-1, caseno-1].set_xlabel('Number of inferred populations')
        if caseno == 1:
            axes[N_true_populations-1, caseno-1].set_ylabel('Negative loglikelihood')
        axes[N_true_populations-1, caseno-1].set_xticks([1,2,3,4])
        #axes[N_true_populations-1, caseno-1].set_ylim([min_neg_ll, max_neg_ll])
        axes[N_true_populations-1, caseno-1].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
        axes[N_true_populations-1, caseno-1].legend()
plt.tight_layout()
plt.savefig("./article/elbowplots/elbow-noise-" + str(Noise) + ".png") #, bbox_inches="tight")
#plt.show()

if saved_gr50s_at_runtime:
    ########## GR50 and mixtures ##########
    #If new version, can print the GR50 and mixture values
    for N_true_populations in [1]: #,2,3]:
        print(N_true_populations)
        for caseno in [1,2,3]:
            print(caseno)
            savestr, MixParams = get_savestr(N_true_populations, caseno, Noise, num_optim)
            filename = './article/CLUSTERPLOTS/' + str(N_true_populations) + '/case' + str(caseno) + '/negative_loglikelihood_values' + savestr + '.mat'
            mat_contents = scipy.io.loadmat(filename)
            negative_loglikelihood_values = mat_contents['negative_loglikelihood_values'][0]
            print('negative_loglikelihood_values')
            print(mat_contents['negative_loglikelihood_values'][0])
            #print('x_finals_temp')
            #print(mat_contents['x_finals_temp'][0])
            #print('true_mixtures')
            #print(mat_contents['true_mixtures'][0])
            #print('inferred_mixtures')
            #print(mat_contents['inferred_mixtures'][0])
            #print('true_gr50values')
            #print(mat_contents['true_gr50values'][0])
            #print('inferred_GR50s')
            #print(["{:.15E}".format(elem) for elem in mat_contents['inferred_GR50s'][0]])
            ##print(f'Five plus ten is {a + b} and not {2 * (a + b)}.', mat_contents['inferred_GR50s'][0])

######################################## Limitations small populations ########################################
def get_savestr_limitations_small_populations(caseno, num_optim):
    Noise = 50
    if caseno > 4:
        N_true_populations = 3
    else:
        N_true_populations = 2
    if caseno == 1:
        MixParams = [0.9]
    elif caseno == 2:
        MixParams = [0.95]
    elif caseno == 3:
        MixParams = [0.97]
    elif caseno == 4:
        MixParams = [0.99]
    elif caseno == 5:
        MixParams = [0.5, 0.45]
    elif caseno == 6:
        MixParams = [0.5, 0.47]
    elif caseno == 7:
        MixParams = [0.5, 0.49]
    elif caseno == 8:
        MixParams = [0.8, 0.1]
    elif caseno == 9:
        MixParams = [0.9, 0.05]    
    if N_true_populations == 1:
        savestr = '-Noise-' + str(Noise) + '-num_optim-' + str(num_optim)
    elif N_true_populations == 2:
        savestr = '-Noise-' + str(Noise) + '-mix-' + str(MixParams[0]) +  '-num_optim-' + str(num_optim)
    elif N_true_populations == 3:
        savestr = '-Noise-' + str(Noise) + '-mix-' + str(MixParams[0]) + '-' + str(MixParams[1]) + '-num_optim-' + str(num_optim)
    return savestr, MixParams, N_true_populations

num_optim = 1000
max_N_populations = 4
saved_gr50s_at_runtime = True

max_neg_ll = 0
min_neg_ll = np.Inf
max_difference = 0
# caseno originally refrred to the slurm_array_id
# Based on this the negLL values were saved in folder 1,2,3 : case1/2/3
# But now this case also refers to place in the list from 1 to 9
# So we need a dict to find the folder name
caseno_lib = {
1 : 3,
2 : 1,
3 : 2,
4 : 3,
5 : 1,
6 : 2,
7 : 3,
8 : 1,
9 : 2,
}

########## Find max values ##########
for caseno in [1,2,3,4,5,6,7,8,9]:
    savestr, MixParams, N_true_populations = get_savestr_limitations_small_populations(caseno, num_optim)
    casefolder_name = caseno_lib[caseno]
    filename = './article/CLUSTERPLOTS/' + str(N_true_populations) + '/case' + str(casefolder_name) + '/negative_loglikelihood_values' + savestr + '.mat'
    mat_contents = scipy.io.loadmat(filename)
    negative_loglikelihood_values = mat_contents['negative_loglikelihood_values'][0]
    if max(negative_loglikelihood_values) > max_neg_ll:
        max_neg_ll = max(negative_loglikelihood_values)
    if min(negative_loglikelihood_values) < min_neg_ll:
        min_neg_ll = min(negative_loglikelihood_values)
    if max(negative_loglikelihood_values) - min(negative_loglikelihood_values) > max_difference:
        max_difference = max(negative_loglikelihood_values) - min(negative_loglikelihood_values)

########## Plot ##########
fig, axes = plt.subplots(nrows = 3, ncols=3, figsize=(10,10))
plt.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0.2, hspace=0.2)
plotno = 0
for caseno in [1,2,3,4,5,6,7,8,9]:
    plotno = plotno + 1
    savestr, MixParams, N_true_populations = get_savestr_limitations_small_populations(caseno, num_optim)
    casefolder_name = caseno_lib[caseno]
    filename = './article/CLUSTERPLOTS/' + str(N_true_populations) + '/case' + str(casefolder_name) + '/negative_loglikelihood_values' + savestr + '.mat'
    #filename = './article/limitations/CLUSTERPLOTS/case-' + str(caseno) + '-negative_loglikelihood_values' + savestr + '.mat'
    mat_contents = scipy.io.loadmat(filename)
    negative_loglikelihood_values = mat_contents['negative_loglikelihood_values'][0]
    if N_true_populations == 1:
        labell = "Truth: " + str(N_true_populations) + " clone"
    elif N_true_populations == 2:
        fullMixParams = MixParams + [1-sum(MixParams)]
        labell = "Case " + str(caseno) + ": " + str(N_true_populations) + " clones\nMixtures {:.2f}".format(fullMixParams[0]) + ", {:.2f}".format(fullMixParams[1])
    elif N_true_populations == 3:
        fullMixParams = MixParams + [1-sum(MixParams)]
        labell = "Case " + str(caseno) + ": " + str(N_true_populations) + " clones\nMixtures {:.2f}".format(fullMixParams[0]) + ", {:.2f}".format(fullMixParams[1]) + ", {:.2f}".format(fullMixParams[2])
    # Find x and y coordinates
    ycoord = int((plotno-1)//3)
    xcoord = np.mod(plotno-1,3)
    axes[ycoord, xcoord].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
    if ycoord == 2:
        axes[ycoord, xcoord].set_xlabel('Number of inferred populations')
    if xcoord == 0:
        axes[ycoord, xcoord].set_ylabel('Negative loglikelihood')
    axes[ycoord, xcoord].set_xticks([1,2,3,4])
    #axes[ycoord, xcoord].set_ylim([min_neg_ll, max_neg_ll])
    axes[ycoord, xcoord].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
    axes[ycoord, xcoord].legend()
    #axes[plotno-1].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
    #if N_true_populations == 3:
    #    axes[plotno-1].set_xlabel('Number of inferred populations')
    #if plotno == 1:
    #    axes[plotno-1].set_ylabel('Negative loglikelihood')
    #axes[plotno-1].set_xticks([1,2,3,4])
    ##axes[plotno-1].set_ylim([min_neg_ll, max_neg_ll])
    #axes[plotno-1].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
    #axes[plotno-1].legend()
plt.tight_layout()
plt.savefig("./article/elbowplots/elbow-limitations-small-populations.pdf") #, bbox_inches="tight")
#plt.show()

######################################## Limitations close GR50 values ########################################
def get_savestr_limitations_close_GR50s(caseno, num_optim):
    Noise = 50
    N_true_populations = 2
    if caseno == 1:
        MixParams = [0.5]
    elif caseno == 2:
        MixParams = [0.5]
    elif caseno == 3:
        MixParams = [0.5]
    elif caseno == 4:
        MixParams = [0.9]
    elif caseno == 5:
        MixParams = [0.9]
    elif caseno == 6:
        MixParams = [0.9]
    elif caseno == 7:
        MixParams = [0.95]
    elif caseno == 8:
        MixParams = [0.95]
    elif caseno == 9:
        MixParams = [0.95]
    savestr = '-Noise-' + str(Noise) + '-mix-' + str(MixParams[0]) +  '-num_optim-' + str(num_optim)
    return savestr, MixParams, N_true_populations

num_optim = 1000
max_N_populations = 4
saved_gr50s_at_runtime = True

max_neg_ll = 0
min_neg_ll = np.Inf
max_difference = 0

########## Find max values ##########
for caseno in [3,2,1,6,5,4,9,8,7]: #[1,2,3,4,5,6,7,8,9]:
    savestr, MixParams, N_true_populations = get_savestr_limitations_close_GR50s(caseno, num_optim)
    filename = './article/limitations/CLUSTERPLOTS/case-' + str(caseno) + '-negative_loglikelihood_values' + savestr + '.mat'
    mat_contents = scipy.io.loadmat(filename)
    negative_loglikelihood_values = mat_contents['negative_loglikelihood_values'][0]
    if max(negative_loglikelihood_values) > max_neg_ll:
        max_neg_ll = max(negative_loglikelihood_values)
    if min(negative_loglikelihood_values) < min_neg_ll:
        min_neg_ll = min(negative_loglikelihood_values)
    if max(negative_loglikelihood_values) - min(negative_loglikelihood_values) > max_difference:
        max_difference = max(negative_loglikelihood_values) - min(negative_loglikelihood_values)

########## Plot ##########
fig, axes = plt.subplots(nrows = 3, ncols=3, figsize=(10,10))
plt.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0.2, hspace=0.2)
plotno = 0
for caseno in [3,2,1,6,5,4,9,8,7]: #[1,2,3,4,5,6,7,8,9]:
    plotno = plotno + 1
    savestr, MixParams, N_true_populations = get_savestr_limitations_close_GR50s(caseno, num_optim)
    filename = './article/limitations/CLUSTERPLOTS/case-' + str(caseno) + '-negative_loglikelihood_values' + savestr + '.mat'
    mat_contents = scipy.io.loadmat(filename)
    negative_loglikelihood_values = mat_contents['negative_loglikelihood_values'][0]
    if N_true_populations == 1:
        labell = "Truth: " + str(N_true_populations) + " clone"
    elif N_true_populations == 2:
        fullMixParams = MixParams + [1-sum(MixParams)]
        labell = "Case " + str(plotno) + ": " + str(N_true_populations) + " clones\nMixtures {:.2f}".format(fullMixParams[0]) + ", {:.2f}".format(fullMixParams[1])
    elif N_true_populations == 3:
        fullMixParams = MixParams + [1-sum(MixParams)]
        labell = "Case " + str(plotno) + ": " + str(N_true_populations) + " clones\nMixtures {:.2f}".format(fullMixParams[0]) + ", {:.2f}".format(fullMixParams[1]) + ", {:.2f}".format(fullMixParams[2])
    # Find x and y coordinates
    ycoord = int((plotno-1)//3)
    xcoord = np.mod(plotno-1,3)
    axes[ycoord, xcoord].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
    if ycoord == 2:
        axes[ycoord, xcoord].set_xlabel('Number of inferred populations')
    if xcoord == 0:
        axes[ycoord, xcoord].set_ylabel('Negative loglikelihood')
    axes[ycoord, xcoord].set_xticks([1,2,3,4])
    #axes[ycoord, xcoord].set_ylim([min_neg_ll, max_neg_ll])
    axes[ycoord, xcoord].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
    axes[ycoord, xcoord].legend()
    #axes[plotno-1].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
    #if N_true_populations == 3:
    #    axes[plotno-1].set_xlabel('Number of inferred populations')
    #if plotno == 1:
    #    axes[plotno-1].set_ylabel('Negative loglikelihood')
    #axes[plotno-1].set_xticks([1,2,3,4])
    ##axes[plotno-1].set_ylim([min_neg_ll, max_neg_ll])
    #axes[plotno-1].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
    #axes[plotno-1].legend()
plt.tight_layout()
plt.savefig("./article/elbowplots/elbow-limitations-close-GR50s.pdf") #, bbox_inches="tight")
#plt.show()

# Old version limitations #############################################
#num_optim = 1000
#max_N_populations = 4
#saved_gr50s_at_runtime = True
#
#max_neg_ll = 0
#min_neg_ll = np.Inf
#max_difference = 0
#
########### Find max values ##########
#for caseno in [1,2,5]: #[1,2,3,4,5]:
#    savestr, MixParams, N_true_populations = get_savestr_limitations(caseno, num_optim)
#    filename = './article/limitations/CLUSTERPLOTS/case-' + str(caseno) + '-negative_loglikelihood_values' + savestr + '.mat'
#    mat_contents = scipy.io.loadmat(filename)
#    negative_loglikelihood_values = mat_contents['negative_loglikelihood_values'][0]
#    if max(negative_loglikelihood_values) > max_neg_ll:
#        max_neg_ll = max(negative_loglikelihood_values)
#    if min(negative_loglikelihood_values) < min_neg_ll:
#        min_neg_ll = min(negative_loglikelihood_values)
#    if max(negative_loglikelihood_values) - min(negative_loglikelihood_values) > max_difference:
#        max_difference = max(negative_loglikelihood_values) - min(negative_loglikelihood_values)
#
########### Plot ##########
#fig, axes = plt.subplots(nrows = 3, ncols=1, figsize=(4,9))
#plt.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0.2, hspace=0.2)
#plotno = 0
#for caseno in [1,2,5]: #[1,2,3,4,5]:
#    plotno = plotno + 1
#    savestr, MixParams, N_true_populations = get_savestr_limitations(caseno, num_optim)
#    filename = './article/limitations/CLUSTERPLOTS/case-' + str(caseno) + '-negative_loglikelihood_values' + savestr + '.mat'
#    mat_contents = scipy.io.loadmat(filename)
#    negative_loglikelihood_values = mat_contents['negative_loglikelihood_values'][0]
#    if N_true_populations == 1:
#        labell = "Truth: " + str(N_true_populations) + " clone"
#    elif N_true_populations == 2:
#        fullMixParams = MixParams + [1-sum(MixParams)]
#        labell = "Truth: " + str(N_true_populations) + " clones\nMixtures {:.2f}".format(fullMixParams[0]) + ", {:.2f}".format(fullMixParams[1])
#    elif N_true_populations == 3:
#        fullMixParams = MixParams + [1-sum(MixParams)]
#        labell = "Truth: " + str(N_true_populations) + " clones\nMixtures {:.2f}".format(fullMixParams[0]) + ", {:.2f}".format(fullMixParams[1]) + ", {:.2f}".format(fullMixParams[2])
#    axes[plotno-1].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
#    if N_true_populations == 3:
#        axes[plotno-1].set_xlabel('Number of inferred populations')
#    if plotno == 1:
#        axes[plotno-1].set_ylabel('Negative loglikelihood')
#    axes[plotno-1].set_xticks([1,2,3,4])
#    #axes[plotno-1].set_ylim([min_neg_ll, max_neg_ll])
#    axes[plotno-1].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
#    axes[plotno-1].legend()
#plt.tight_layout()
#plt.savefig("./article/elbowplots/elbow-limitations.png") #, bbox_inches="tight")
##plt.show()


######################################## Shannon ########################################

cell_line_825 = [303.563, 301.7188, 301.7188, 301.7188]
cell_line_1975 = [350.4950, 350.4950, 350.4950, 350.4950]
#eightytwenty_shannon = 100*np.array([3.453952579665839, 3.452474898215619, 3.452477200002798, 3.452721855864264]) #, 3.453194893597295, 3.452523656154327] # 14
eightytwenty_shannon = 100*np.array([3.225027912916597, 3.094875476185684, 3.084815650911908, 3.084820237386532]) #, 3.086597651098093, 3.096281523059681] # 41
fiftyfifty_shannon = 100*np.array([3.287210629087783, 3.109920381548262, 3.101590402418301, 3.101591160234995]) #, 3.100805124104205, 3.106348364096599]

max_difference = 20.0
max_N_populations = 4
fig, axes = plt.subplots(nrows = 4, ncols=1, figsize=(3,10))
########## Plot ##########
plt.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0.2, hspace=0.2)
labell = "Truth: 1 clone"

negative_loglikelihood_values = cell_line_1975
axes[0].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
axes[0].set_xlabel('Number of inferred populations')
axes[0].set_ylabel('Negative loglikelihood')
axes[0].set_xticks([1,2,3,4])
axes[0].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
axes[0].legend()

negative_loglikelihood_values = cell_line_825
axes[1].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
axes[1].set_xlabel('Number of inferred populations')
axes[1].set_ylabel('Negative loglikelihood')
axes[1].set_xticks([1,2,3,4])
axes[1].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
axes[1].legend()

negative_loglikelihood_values = eightytwenty_shannon
axes[2].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
axes[2].set_xlabel('Number of inferred populations')
axes[2].set_ylabel('Negative loglikelihood')
axes[2].set_xticks([1,2,3,4])
axes[2].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
axes[2].legend()

negative_loglikelihood_values = fiftyfifty_shannon
axes[3].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
axes[3].set_xlabel('Number of inferred populations')
axes[3].set_ylabel('Negative loglikelihood')
axes[3].set_xticks([1,2,3,4])
axes[3].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
axes[3].legend()

plt.tight_layout()
plt.savefig("./article/elbowplots/elbow-shannon.png") #, bbox_inches="tight")
plt.show()


######################################## Dagim ########################################

sens_500 = [1.3305e4, 1.3238e4, 1.3234e4, 1.3234e4]
res_500 = [1.2524e4, 1.2413e4, 1.2368e4, 1.2365e4]
eighty_twenty = [13325.5074194264, 12900.2574627944, 12853.4825068622, 12841.3786236734]
sixtyseven_thirtythree = [12687.5198352985, 12063.6612358055, 12051.5944422424, 12049.0141573515]
fifty_fifty = [12481.1898524187, 12081.484858823, 12013.488249187, 12002.7119377631]

#[1.1757e4, 1.1338e4, 1.126e4, 1.1212e4] # Resistant 250
#33 67 : [13308.2462670243, 13024.6077094253, 12938.5640479326, 12922.2316240557]

########## Find max values ##########
max_neg_ll = 0
min_neg_ll = np.Inf
max_difference = 0
negllarray = [sens_500, res_500, eighty_twenty, sixtyseven_thirtythree, fifty_fifty]
for ii in range(len(negllarray)):
    negative_loglikelihood_values = negllarray[ii]
    if max(negative_loglikelihood_values) > max_neg_ll:
        max_neg_ll = max(negative_loglikelihood_values)
    if min(negative_loglikelihood_values) < min_neg_ll:
        min_neg_ll = min(negative_loglikelihood_values)
    if max(negative_loglikelihood_values) - min(negative_loglikelihood_values) > max_difference:
        max_difference = max(negative_loglikelihood_values) - min(negative_loglikelihood_values)

fig, axes = plt.subplots(nrows = 5, ncols=1, figsize=(3,12))
max_N_populations = 4
########## Plot ##########
plt.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0.2, hspace=0.2)

labell = "Only sensitive"
negative_loglikelihood_values = sens_500
axes[0].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
axes[0].set_xlabel('Number of inferred populations')
axes[0].set_ylabel('Negative loglikelihood')
axes[0].set_xticks([1,2,3,4])
axes[0].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
axes[0].legend()

labell = "Only resistant"
negative_loglikelihood_values = res_500
axes[1].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
axes[1].set_xlabel('Number of inferred populations')
axes[1].set_ylabel('Negative loglikelihood')
axes[1].set_xticks([1,2,3,4])
axes[1].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
axes[1].legend()

labell = "80 / 20"
negative_loglikelihood_values = eighty_twenty
axes[2].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
axes[2].set_xlabel('Number of inferred populations')
axes[2].set_ylabel('Negative loglikelihood')
axes[2].set_xticks([1,2,3,4])
axes[2].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
axes[2].legend()

labell = "67 / 33"
negative_loglikelihood_values = sixtyseven_thirtythree
axes[3].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
axes[3].set_xlabel('Number of inferred populations')
axes[3].set_ylabel('Negative loglikelihood')
axes[3].set_xticks([1,2,3,4])
axes[3].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
axes[3].legend()

labell = "50 / 50"
negative_loglikelihood_values = fifty_fifty
axes[4].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
axes[4].set_xlabel('Number of inferred populations')
axes[4].set_ylabel('Negative loglikelihood')
axes[4].set_xticks([1,2,3,4])
axes[4].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
axes[4].legend()

plt.tight_layout()
plt.savefig("./article/elbowplots/elbow-dagim.png") #, bbox_inches="tight")
plt.show()

######################################## Real MM patients ########################################

#patno_array = ["MM210819", "MM130120", "MM19520", "MM3620", "MM1420"] # from MM inference scripts
patno_array = ["MM1420", "MM210819", "MM130120", "MM19520", "MM3620"] # Patients 1,2,3,4,5 with ids 0,1,2,3,4

########## Find max values ##########
max_neg_ll = 0
min_neg_ll = np.Inf
max_difference = 0
for patient_index in [1,2,3,4,5]:
    patient_number = patno_array[patient_index-1]
    for drugnumber in [1,2,3,4]:
        filename = './data/real-MM-cellines/negative_loglikelihood_values-' + str(patient_number) + '-drug-' + str(drugnumber) + '.mat'
        mat_contents = scipy.io.loadmat(filename)
        negative_loglikelihood_values = mat_contents['negative_loglikelihood_values'][0]
        if max(negative_loglikelihood_values) > max_neg_ll:
            max_neg_ll = max(negative_loglikelihood_values)
        if min(negative_loglikelihood_values) < min_neg_ll:
            min_neg_ll = min(negative_loglikelihood_values)
        if max(negative_loglikelihood_values) - min(negative_loglikelihood_values) > max_difference:
            max_difference = max(negative_loglikelihood_values) - min(negative_loglikelihood_values)

########## Plot ##########
max_N_populations = 5
fig, axes = plt.subplots(nrows = 4, ncols=5, figsize=(15,12))
plt.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0.2, hspace=0.2)
for patient_index in [1,2,3,4,5]:
    patient_number = patno_array[patient_index-1]
    if patient_number == "MM210819":
        drug_list = ["Ixazomib", "Thalidomide", "Selinexor", "Dexamethasone"]
    elif patient_number == "MM130120":
        drug_list = ["Ixazomib", "Thalidomide", "Selinexor", "Dexamethasone"]
    elif patient_number == "MM19520":
        drug_list = ["Venetoclax", "Melflufen", "Selinexor", "Dexamethasone"]
    elif patient_number == "MM3620":
        drug_list = ["Venetoclax", "Melflufen", "Selinexor", "Dexamethasone"]
    elif patient_number == "MM1420":
        drug_list = ["Venetoclax", "Melflufen", "Selinexor", "Dexamethasone"]
    for drugnumber in [1,2,3,4]:
        filename = './data/real-MM-cellines/negative_loglikelihood_values-' + str(patient_number) + '-drug-' + str(drugnumber) + '.mat'
        mat_contents = scipy.io.loadmat(filename)
        negative_loglikelihood_values = mat_contents['negative_loglikelihood_values'][0]
        labell = drug_list[drugnumber-1]
        axes[drugnumber-1, patient_index-1].plot(range(1,max_N_populations+1), negative_loglikelihood_values, '.-k', label=labell)
        if drugnumber == 1: 
            axes[drugnumber-1, patient_index-1].set_title("Patient " + patient_number)
        if drugnumber == 4:
            axes[drugnumber-1, patient_index-1].set_xlabel('Number of inferred populations')
        if patient_index == 1:
            axes[drugnumber-1, patient_index-1].set_ylabel('Negative loglikelihood')
        axes[drugnumber-1, patient_index-1].set_xticks([1,2,3,4])
        axes[drugnumber-1, patient_index-1].set_ylim([min(negative_loglikelihood_values) - 0.1*max_difference, min(negative_loglikelihood_values) + max_difference*1.2])
        axes[drugnumber-1, patient_index-1].legend()
plt.tight_layout()
plt.savefig("./article/elbowplots/elbow-MM-patients.png") #, bbox_inches="tight")
#plt.show()

