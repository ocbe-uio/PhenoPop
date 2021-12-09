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
from matplotlib import lines

labels_all_patients = [
"True clone 1",
"Estimated clone 1",
"True clone 2",
"Estimated clone 2",
"True clone 3",
"Estimated clone 3",
]

legendcolors = [
desaturate(paul_tol[0], 0.4),
paul_tol[0],
desaturate(paul_tol[1], 0.4),
paul_tol[1],
desaturate(paul_tol[2], 0.4),
paul_tol[2],
]

legend_elements = [Patch(label=labels_all_patients[ii], color=legendcolors[ii]) for ii in range(len(labels_all_patients))]
bbox_to_anchor = (0, 0.5, 1, 0.5) # bbox:(x y width height) 

# Create the figure for true and estimated clones
fig, ax = plt.subplots()
fig.set_figheight(0.6)
fig.set_figwidth(6)
ax.legend(handles=legend_elements, bbox_to_anchor=bbox_to_anchor, loc='center', ncol=3, mode="expand", borderaxespad=0.0, frameon=False) # upper left means upper left position is fixed. 

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)
plt.tight_layout()
plt.show()


# Create the figure for observed concentrations
vertical_line = lines.Line2D([], [], color='grey', marker='|', linestyle='None',
                          markersize=10, markeredgewidth=1.5, label='Observed concentration')
bbox_to_anchor = (0, 0.5, 0.5, 0.5) # bbox:(x y width height) 

fig, ax = plt.subplots()
fig.set_figheight(0.6)
fig.set_figwidth(6)
ax.legend(handles=[vertical_line], bbox_to_anchor=bbox_to_anchor, loc='center', ncol=3, mode="expand", borderaxespad=0.0, frameon=False) # upper left means upper left position is fixed. 

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)
plt.tight_layout()
plt.show()
