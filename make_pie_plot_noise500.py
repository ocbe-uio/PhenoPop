import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import colorsys
import matplotlib.colors as mc
import matplotlib.patches as patches
import matplotlib
from moodboard import *
import matplotlib.patheffects as pe
from matplotlib.patches import Patch

def make_pie_plot_noise500(N_cases, N_inferred_populations, true_mixtures, estimated_mixtures, estimated_gr50s, true_gr50_for_plotting, mixturecolors, gr50colors_inferred, vlinerange, title, savename, Conc, KNOW_TRUTH=True, concticks=[], pctdistance=1.5, textprops=None, figsize=[], xlabel="None", truegr50regionlimits=[], estimatedgr50regionlimits=[], bbox_to_anchor=(1.03, 0, 0.57, 1), ax_set_position=[0.6,0.1,0.35,0.9], logaxis=True, titlepadding=0.0, wspace=0.2, hspace=0.0, labels_patients = [], mixturecolors_light=[], arrow_multiplier=3, ylim=[], gr50colors_true=[], true_gr50_indices=[], legendcolors=[], legendcolors_inferred=[], limitationflag=0, skewflag=0, dagimflag=0, shannonflag=0, width_ratios=[1,1,7]):
    if len(gr50colors_true) < 1:
        gr50colors_true = gr50colors_inferred
    if len(legendcolors_inferred) < 1:
        legendcolors_inferred = legendcolors
    if KNOW_TRUTH:
        if len(figsize) > 0:
            fig, axes = plt.subplots(nrows=N_cases, ncols = 3, figsize=(figsize[0], figsize[1]), gridspec_kw={"width_ratios":width_ratios})
        else:
            fig, axes = plt.subplots(nrows=N_cases, ncols = 3, gridspec_kw={"width_ratios":width_ratios})
        plt.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=wspace, hspace=hspace)
        axes[0,2].set_title(title, loc="right", pad=titlepadding, fontweight='bold', fontsize=19)
        labels_est = ['Estimated GR50 clone 1', 'Estimated GR50 clone 2', 'Estimated GR50 clone 3']
        labels_conc = "Observed concentration"
        labels_true_gr50 = ['True GR50 clone 1', 'True GR50 clone 2', 'True GR50 clone 3']
        region_labels = ['True GR50 region clone 1', 'True GR50 region clone 2', 'True GR50 region clone 3', 'True GR50 region clone 4', 'True GR50 region clone 5']
        max_x_value = 0

        for ii in range(N_cases):
            # Plot true pies 
            if ii == 3:
                gr50colors_pie_true = [desaturate(color, 0.5) for color in mixturecolors]
                axes[ii,0].pie(true_mixtures[ii], colors=gr50colors_pie_true, autopct='%1.0f%%', pctdistance=pctdistance, startangle=60, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            elif ii == 4:
                gr50colors_pie_true = [desaturate(color, 0.5) for color in mixturecolors]
                axes[ii,0].pie(true_mixtures[ii], colors=gr50colors_pie_true, autopct='%1.0f%%', pctdistance=pctdistance, startangle=22, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            elif ii == 5:
                gr50colors_pie_true = [desaturate(color, 0.5) for color in mixturecolors]
                axes[ii,0].pie(true_mixtures[ii], colors=gr50colors_pie_true, autopct='%1.0f%%', pctdistance=pctdistance, startangle=90, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            elif ii == 6:
                gr50colors_pie_true = [desaturate(color, 0.5) for color in mixturecolors]
                axes[ii,0].pie(true_mixtures[ii], colors=gr50colors_pie_true, autopct='%1.0f%%', pctdistance=pctdistance, startangle=90, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            elif ii == 7:
                gr50colors_pie_true = [desaturate(color, 0.5) for color in mixturecolors]
                axes[ii,0].pie(true_mixtures[ii], colors=gr50colors_pie_true, autopct='%1.0f%%', pctdistance=pctdistance, startangle=109-40, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            elif ii == 8:
                gr50colors_pie_true = [desaturate(color, 0.5) for color in mixturecolors]
                axes[ii,0].pie(true_mixtures[ii], colors=gr50colors_pie_true, autopct='%1.0f%%', pctdistance=pctdistance, startangle=74-35, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            else:
                gr50colors_pie_true = [desaturate(color, 0.5) for color in gr50colors_true[ii]]
                axes[ii,0].pie(true_mixtures[ii], colors=gr50colors_pie_true, autopct='%1.0f%%', pctdistance=pctdistance, startangle=90, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            # Plot inferred pies 
            if ii == 3:
                axes[ii,1].pie(estimated_mixtures[ii], colors=gr50colors_inferred[ii], autopct='%1.0f%%', pctdistance=pctdistance, startangle=45, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            elif ii == 4:
                axes[ii,1].pie(estimated_mixtures[ii], colors=gr50colors_inferred[ii], autopct='%1.0f%%', pctdistance=pctdistance, startangle=13, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            elif ii == 5:
                axes[ii,1].pie(estimated_mixtures[ii], colors=gr50colors_inferred[ii], autopct='%1.0f%%', pctdistance=pctdistance, startangle=70, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            elif ii == 6:
                axes[ii,1].pie(estimated_mixtures[ii], colors=gr50colors_inferred[ii], autopct='%1.0f%%', pctdistance=pctdistance, startangle=50, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            elif ii == 7:
                axes[ii,1].pie(estimated_mixtures[ii], colors=gr50colors_inferred[ii], autopct='%1.0f%%', pctdistance=pctdistance, startangle=35, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            elif ii == 8:
                axes[ii,1].pie(estimated_mixtures[ii], colors=gr50colors_inferred[ii], autopct='%1.0f%%', pctdistance=pctdistance, startangle=10, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)
            else:
                axes[ii,1].pie(estimated_mixtures[ii], colors=gr50colors_inferred[ii], autopct='%1.0f%%', pctdistance=pctdistance, startangle=90, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)

            # Plot GR50 
            if xlabel != "None":
                axes[ii,2].set_xlabel(xlabel, size=19)
            else:
                axes[ii,2].set_xlabel(r"$\rm Simulated \enspace drug \enspace concentration \enspace (\mu M$)", size=19)

            # Mark true GR50s only for those populations present in this case
            # With regions based on concentrations:
            Npoints = 2
            dummy_index = 0
            for jj in true_gr50_indices[ii]: #range(len(truegr50regionlimits)):
                if limitationflag == 1:
                    color = mixturecolors[dummy_index]
                    dummy_index = dummy_index + 1
                else:
                    color = mixturecolors[jj]
                c = desaturate(color, 0.4)
                fill_x = np.linspace(truegr50regionlimits[jj][0], truegr50regionlimits[jj][1], num=Npoints)
                y1 = np.repeat(vlinerange[0], Npoints)
                y2 = np.repeat(vlinerange[1], Npoints)
                legends1 = axes[ii,2].fill_between(fill_x, y1, y2, color=c, label=region_labels[jj])

            ## With lines:
            #for jj in range(len(truegr50regionlimits)):
            #    legends2 = axes[ii,2].vlines(true_gr50_for_plotting[jj], ymin=vlinerange[0], ymax=vlinerange[1], colors=mixturecolors[jj], linestyles='-', label=labels_true_gr50[jj], linewidth=3, zorder=1000) #path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()],
            # Scatter estimated values
            for kk in range(N_inferred_populations[ii]):
                est = estimated_gr50s[ii][kk]
                if est == np.Inf:
                    # GR50 outside observed concentration range
                    bar = [[0], [arrow_multiplier*max(Conc)]]
                    # Line:
                    axes[ii,2].errorbar(max(Conc), N_cases/2, xerr=bar, color=gr50colors_inferred[ii][kk], marker='', linestyle='', elinewidth=30, linewidth=1, zorder=4)
                    # Arrow:
                    tweak = 1.95
                    axes[ii,2].plot(tweak*arrow_multiplier*max(Conc), N_cases/2+0.012, color=gr50colors_inferred[ii][kk], marker='>', linestyle='', markersize=28, markeredgecolor=gr50colors_inferred[ii][kk], zorder=20000)
                    if tweak*arrow_multiplier*max(Conc) > max_x_value:
                        max_x_value = tweak*arrow_multiplier*max(Conc)
                    #legends3 = axes[ii,2].plot(Conc[-1]*arrow_multiplier, (N_cases)/2, color=gr50colors_inferred[ii][kk], marker='>', linestyle='', markersize=20, markeredgecolor=gr50colors_inferred[ii][kk], zorder=20000)
                else:
                    est = estimated_gr50s[ii][kk]
                    bar_lo = abs(est - estimatedgr50regionlimits[ii][kk][0])
                    bar_hi = abs(est - estimatedgr50regionlimits[ii][kk][1])
                    bar = [[bar_lo], [bar_hi]]
                    if ii==N_cases-1:
                        legends3 = axes[ii,2].errorbar(est, N_cases/2, xerr=bar, color=gr50colors_inferred[ii][kk], marker='o', markeredgewidth=2, markeredgecolor='k', markerfacecolor='white', linestyle='', markersize=10, elinewidth=30, linewidth=1, zorder=4, label=labels_est[kk])
                    else:
                        axes[ii,2].errorbar(est, N_cases/2, xerr=bar, color=gr50colors_inferred[ii][kk], marker='o', markeredgewidth=2, markeredgecolor='k', markerfacecolor='white', linestyle='', markersize=10, elinewidth=30, linewidth=1, zorder=4)
            axes[ii,2].spines['top'].set_visible(False)
            axes[ii,2].spines['right'].set_visible(False)
            axes[ii,2].spines['bottom'].set_visible(False)
            axes[ii,2].spines['left'].set_visible(False)
            if logaxis==True:
                axes[ii,2].set_xscale('log')
            if len(Conc) > 0:
                legends4 = axes[ii,2].vlines(Conc, ymin=vlinerange[0], ymax=vlinerange[1], colors='grey', alpha=0.7, linestyles='-', linewidth=1.7) #, label=labels_conc
            if len(concticks) > 0:
                axes[ii,2].set_xticks(concticks)
                ticklabels = [0] + ['$10^{' + format(np.log10(elem), ".0f") + '}$' for elem in concticks[1:len(concticks)]]
                axes[ii,2].set_xticklabels(ticklabels)
            plt.minorticks_off()
            axes[ii,2].get_yaxis().set_visible(False)
            if ii != N_cases-1:
                axes[ii,2].get_xaxis().set_visible(False)
            # Set x limit xlim set_xlim
            axes[ii,1].plot(max_x_value, N_cases/2, marker='.', color="white", zorder=0)
            # Set y limit
            if len(ylim) > 0:
                axes[ii,2].set_ylim(bottom=ylim[0], top=ylim[1])

        ## Legend
        #if len(legendcolors)>1:
        #    legend_elements_1 = [Patch(label='True clone '+ format(ii+1, ".0f"), color=desaturate(color, 0.4)) for ii, color in enumerate(legendcolors)]
        #    legend_elements_2 = [Patch(label='Estimated clone '+ format(ii+1, ".0f"), color=elem) for ii,elem in enumerate(legendcolors_inferred)]
        #    legend_elements = legend_elements_1 + legend_elements_2
        #    axes[N_cases-1,2].legend(handles=legend_elements, bbox_to_anchor=bbox_to_anchor, loc='upper left', ncol=2, mode="expand", borderaxespad=0.0, frameon=False) # bbox:(x y width height) # upper left means upper left position is fixed. 
        #else:
        #    axes[N_cases-1,2].legend(bbox_to_anchor=bbox_to_anchor, loc='upper left', borderaxespad=0.0, frameon=False) # bbox:(x y width height) # upper left means upper left position is fixed. 
    else:
        legend_elements = []
        if len(figsize) > 0:
            fig, axes = plt.subplots(nrows=N_cases, ncols = 2, figsize=(figsize[0], figsize[1]), gridspec_kw={"width_ratios":[1,7]})
        else:
            fig, axes = plt.subplots(nrows=N_cases, ncols = 2, gridspec_kw={"width_ratios":[1,7]})
        plt.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=wspace, hspace=hspace)
        axes[0,1].set_title('Estimated mixture                            GR50 values                                       ', loc="right", pad=titlepadding, fontweight='bold', fontsize=19)
        labels_conc = "Observed concentration"
        max_x_value = 0
        for ii in range(N_cases):
            this_max_x = Conc[-1]*arrow_multiplier
            if this_max_x > max_x_value:
                max_x_value = this_max_x
        for ii in range(N_cases):
            # Plot inferred pies 
            axes[ii,0].pie(estimated_mixtures[ii], colors=gr50colors_inferred[ii], autopct='%1.0f%%', pctdistance=pctdistance, startangle=90, wedgeprops = {'linewidth': 1, 'edgecolor' : 'k'}, textprops=textprops)

            # Plot GR50 
            if xlabel != "None":
                axes[ii,1].set_xlabel(xlabel)
            else:
                axes[ii,1].set_xlabel("Drug concentration")
            # Scatter estimated GR50 values
            for kk in range(N_inferred_populations[ii]):
                legend_elements = legend_elements + [Patch(label=labels_patients[ii][kk], color=gr50colors_inferred[ii][kk])]
                est = estimated_gr50s[ii][kk]
                if est == np.Inf:
                    # GR50 outside observed concentration range
                    bar = [[0], [arrow_multiplier*max(Conc)]]
                    # Line:
                    axes[ii,1].errorbar(max(Conc), N_cases/2, xerr=bar, color=gr50colors_inferred[ii][kk], marker='', linestyle='', elinewidth=30, linewidth=1, zorder=4, label=labels_patients[ii][kk])
                    # Arrow:
                    tweak = 2.08
                    axes[ii,1].plot(tweak*arrow_multiplier*max(Conc), N_cases/2+0.012, color=gr50colors_inferred[ii][kk], marker='>', linestyle='', markersize=28, markeredgecolor=gr50colors_inferred[ii][kk], zorder=20000, label=labels_patients[ii][kk])
                    if tweak*arrow_multiplier*max(Conc) > max_x_value:
                        max_x_value = tweak*arrow_multiplier*max(Conc)
                else:
                    bar_lo = abs(est - estimatedgr50regionlimits[ii][kk][0])
                    bar_hi = abs(est - estimatedgr50regionlimits[ii][kk][1])
                    bar = [[bar_lo], [bar_hi]]
                    #if ii==N_cases-1:
                    legends3 = axes[ii,1].errorbar(est, N_cases/2, xerr=bar, color=gr50colors_inferred[ii][kk], marker='o', markeredgewidth=2, markeredgecolor='k', markerfacecolor='white', linestyle='', markersize=10, elinewidth=30, linewidth=1, zorder=4, label=labels_patients[ii][kk])
                    #else:
                    #    axes[ii,1].errorbar(est, N_cases/2, xerr=bar, color=gr50colors_inferred[ii][kk], marker='o', linestyle='', markersize=8, markeredgecolor='black', elinewidth=3, linewidth=1, zorder=4)
            axes[ii,1].spines['top'].set_visible(False)
            axes[ii,1].spines['right'].set_visible(False)
            axes[ii,1].spines['bottom'].set_visible(False)
            axes[ii,1].spines['left'].set_visible(False)
            if logaxis==True:
                axes[ii,1].set_xscale('log')
            if len(Conc) > 0:
                legends4 = axes[ii,1].vlines(Conc, ymin=vlinerange[0], ymax=vlinerange[1], colors='grey', alpha=0.7, linestyles='-', linewidth=1.7) #label=labels_conc
            if len(concticks) > 0:
                axes[ii,1].set_xticks(concticks)
                ticklabels = [0] + ['$10^{' + format(np.log10(elem*10**-3), ".0f") + '}$' for elem in concticks[1:len(concticks)]]
                axes[ii,1].set_xticklabels(ticklabels)
            plt.minorticks_off()
            axes[ii,1].get_yaxis().set_visible(False)
            # Set x limit xlim set_xlim
            axes[ii,1].plot(max_x_value, N_cases/2, marker='.', color="white", zorder=0)
            # Set y limit
            if len(ylim) > 0:
                axes[ii,1].set_ylim(bottom=ylim[0], top=ylim[1])
            if ii != N_cases-1:
                axes[ii,1].get_xaxis().set_visible(False)
        ## Legend
        #axes[N_cases-1,1].legend(handles=legend_elements, bbox_to_anchor=bbox_to_anchor, loc='upper left', ncol=2, mode="expand", borderaxespad=0.0, frameon=False) # bbox:(x y width height) # upper left means upper left position is fixed. 
    #plt.tight_layout()
    plt.savefig(savename, bbox_inches="tight")
    #plt.show()
