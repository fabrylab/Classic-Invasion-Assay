# coding: utf-8

from __future__ import division, print_function
import numpy as np
from glob import glob
import csv
import matplotlib.pyplot as plt
import os
from collections import defaultdict

# Insert key facts here
celltype = 'U87'
cond = 'Alginate'
# Insert path of cell-type folder here
folder=r'H:\Experiment_data\B01_AnWi_Invasion_2020-01-31\U87'




data = []
dates = []
dist_to_gel=20  # in um
z_slice_thickness=5*1.33
# Searches for the indices txt file
files = glob(folder+'/*/*/' + celltype + '_' + cond + '*.txt')

old_date = ''
old_pos = ''
#files2=["a.txt","b.txt"]
#files=[os.path.join(folder,x) for x in files2]
files = [r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\CNTRL\pos00\xyz_positions.txt",
         r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\CNTRL\pos01\xyz_positions.txt",
         r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\CNTRL\pos02\xyz_positions.txt",
         r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\CNTRL\pos03\xyz_positions.txt",
         r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\CNTRL\pos04\xyz_positions.txt",
         r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\CNTRL\pos05\xyz_positions.txt",
         r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\CNTRL\pos06\xyz_positions.txt",
         r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\CNTRL\pos07\xyz_positions.txt",
         r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\CNTRL\pos08\xyz_positions.txt"]


lables=["pos00", "pos01", "pos02", "pos03", "pos04", "pos05", "pos05", "pos05", "pos05", "pos05", "pos05", "pos05", "pos05", "pos05"]

for f in files:
    raw=[]
    with open(f, 'r') as txt_file:
        next(txt_file)
        for l in txt_file.readlines():
            raw.append(float(l.strip().split(" ")[-1]))
    # projecting cells above the gel surface to the gel surface
    raw = np.array([z if z <= 100 else 100 for z in raw])
    # adding to dat dictionary
    data.append(np.sort(raw))





colors = [p['color'] for p in plt.rcParams['axes.prop_cycle']]

for i,(data, label) in enumerate(zip(data,lables)):
    data=-(data - 100) * z_slice_thickness
    p = np.linspace(1, 0, len(data))
    ps = p * 100 # percentiles
    pooled_ps = np.array([np.percentile(data, i, interpolation="midpoint") for i in ps])
    plt.plot(p[::-1], pooled_ps, label=label, color=colors[i])

plt.axhline(0, c='k', lw=1)
plt.grid()

#plt.xscale('log')
plt.xlim([10**0, 10**-4])
plt.ylim([201*2*1.33, -50])
plt.title(celltype)
#plt.gca().spines["right"].set_visible(False)
#plt.gca().spines["top"].set_visible(False)
plt.ylabel('Invasion depth D [µm]')
plt.xlabel('Probability of (Invasion depth$\mathregular{\geq}$D)')
plt.legend(loc='lower left')
plt.tight_layout()
#plt.show()
plt.savefig(os.path.join("H:\Experiment_data\B01_AnWi_Invasion_2020-03-04", 'test.png'), dpi=200)

### merging all