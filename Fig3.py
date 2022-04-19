#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 13:42:19 2022
Script to plot figure 3 in article of JGR: Atmospheres, titled:

"Experimental radial profiles of early time ($<$ 4 $\mu$s) neutral and ion spectroscopic signatures in lightning-like discharges"

Author: 
    María Passas Varo
Email:
    passasv_at_iaa.es
Affiliation: 
    Instituto de Astrofísica de Andalucía IAA-CSIC
    Glorieta de la Astronomía sn
    18016 Granada
    Spain
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


dirname = os.getcwd() + '/Data'
alpha = 0.58                # Spatial dispersion 0.58 mm/px 


VI = pd.read_csv(dirname + '/Fig3.csv') 

t = VI['Time(s)']
V = VI['V(V)']
I = VI['I(A)']

fig, ax1 = plt.subplots(figsize=(8,5))    
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'
color2='darkslategrey'

ax2.set_ylabel('Intensity (A)', color=color, fontsize=14)  # we already handled the x-label with ax1
ax2.plot(t*1e6,I,'.-', alpha=0.5, markersize= 2, color=color)


ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
ax1.set_xlim(-1,4)
ax2.set_xlim(-1,4)

color = 'tab:blue'
ax1.set_xlabel('Time ($\mu$s)', fontsize=14)
ax1.set_ylabel('Voltage (V)', color=color, fontsize=14)
ax1.plot(t*1e6,V, '.-', alpha=0.5,markersize= 2, color=color)
ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
ax1.tick_params(axis='x', labelsize=12)

idx_ini=np.argwhere(t<0.72e-6)[-1][0]

a=np.arange(0.72, 1.51, 0.001)
ax1.fill_between(a , 0, 30000, facecolor= color2, hatch='..', alpha=0.1,zorder=1)
a=np.arange(1.83, 2.62, 0.001)
ax1.fill_between(a , 0, 30000, facecolor= color2, hatch='..', alpha=0.1,zorder=1)
a=np.arange(2.94, 3.73, 0.001)
ax1.fill_between(a , 0, 30000, facecolor= color2, hatch='..', alpha=0.1,zorder=1)

plt.tight_layout()
plt.show()


