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

labels_all_patients = [
"Patient MM1420 clone 1", 
"Patient MM1420 clone 2", 
"Patient MM210819 clone 1", 
"Patient MM130120 clone 1", 
"Patient MM19520 clone 1", 
"Patient MM19520 clone 2", 
"Patient MM3620 clone 1", 
"Patient MM3620 clone 2",
]

legendcolors = [
paul_tol_darker[3],
paul_tol_lighter[3],
paul_tol_darker[4],
paul_tol_darker[5],
paul_tol_darker[6],
paul_tol_lighter[6],
paul_tol_darker[7],
paul_tol_lighter[7],
]

legend_elements = [Patch(label=labels_all_patients[ii], color=legendcolors[ii]) for ii in range(len(labels_all_patients))]
bbox_to_anchor = (0, 0.5, 1, 0.5) # bbox:(x y width height) 


# Create the figure
fig, ax = plt.subplots()
fig.set_figheight(4)
fig.set_figwidth(10)
ax.legend(handles=legend_elements, bbox_to_anchor=bbox_to_anchor, loc='center', ncol=4, mode="expand", borderaxespad=0.0, frameon=False) # upper left means upper left position is fixed. 

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)
plt.tight_layout()
plt.show()
