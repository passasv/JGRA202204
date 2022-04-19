#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 13:42:19 2022
Script to plot figure 10 in article of JGR: Atmospheres, titled:

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
from scipy.interpolate import interp1d


dirname = os.getcwd() + '/Data'

colors = plt.cm.viridis(np.linspace(0,1,5))
colors2 = ['black', 'red', 'orange']

#=======================================================
# Panel 10a
#=======================================================

data = pd.read_csv(dirname + '/Fig10a.csv') 


r = np.asarray(data['Radius(mm)'])
T1 = np.asarray(data['T_1.11us(K)'])
T2 = np.asarray(data['T_2.22us(K)'])
T3  = np.asarray(data['T_3.33us(K)'])
B1 = np.asarray(data['B_1.11us(a.u.)'])
B2  = np.asarray(data['B_2.22us(a.u.)'])
B3 = np.asarray(data['B_3.33us(a.u.)'])

#=======================================================
# Panel 10b
#=======================================================

data2 = pd.read_csv(dirname + '/Fig10b.csv') 

N1 = np.asarray(data2['Ne_1.11us(K)'])
N2 = np.asarray(data2['Ne_2.22us(K)'])
N3 = np.asarray(data2['Ne_3.33us(K)'])


#=======================================================
# Panel 10c
#=======================================================
data3 = pd.read_csv(dirname + '/Fig10c.csv') 

C1 = np.asarray(data3['C_1.11us(K)'])
C2 = np.asarray(data3['C_2.22us(K)'])
C3  = np.asarray(data3['C_3.33us(K)'])


#=======================================================
# Panel 10d
#=======================================================

data4 = pd.read_csv(dirname + '/Fig10d.csv') 

D1 = np.asarray(data4['d_1.11us(K)'])
D2 = np.asarray(data4['d_2.22us(K)'])
D3  = np.asarray(data4['d_3.33us(K)'])

#=======================================================
# Plots
#=======================================================

fig, axs = plt.subplots(2,2, figsize=(9,8))

ax0 = axs[0,0].twinx() 
ax0.fill_between(r, 0, B1, facecolor=colors[0],alpha=0.1, label='0.00 $\mu$s $<$ t $\leq$ 1.11 $\mu$s',zorder=4) 
ax0.fill_between(r, 0, B2, facecolor=colors[1],alpha=0.2, label='1.11 $\mu$s $<$ t $\leq$ 2.22 $\mu$s',zorder=4) 
ax0.fill_between(r, 0, B3, facecolor=colors[2],alpha=0.3, label='2.22 $\mu$s $<$ t $\leq$ 3.33 $\mu$s',zorder=4) 

ax0.set_ylabel('Brightness (a.u.)', fontsize=14)

ax1 = axs[0,1].twinx() 
ax1.fill_between(r, 0, B1, facecolor=colors[0], alpha=0.1, label='0.00 $\mu$s $<$ t $\leq$ 1.11 $\mu$s',zorder=4) 
ax1.fill_between(r, 0, B2, facecolor=colors[1], alpha=0.2, label='1.11 $\mu$s $<$ t $\leq$ 2.22 $\mu$s',zorder=4) 
ax1.fill_between(r, 0, B3, facecolor=colors[2], alpha=0.3, label='2.22 $\mu$s $<$ t $\leq$ 3.33 $\mu$s',zorder=4) 

ax1.set_ylabel('Brightness (a.u.)', fontsize=14)

ax2 = axs[1,0].twinx() 
ax2.fill_between(r, 0, B1, facecolor=colors[0], alpha=0.1, label='0.00 $\mu$s $<$ t $\leq$ 1.11 $\mu$s',zorder=4) 
ax2.fill_between(r, 0, B2, facecolor=colors[1], alpha=0.2, label='1.11 $\mu$s $<$ t $\leq$ 2.22 $\mu$s',zorder=4) 
ax2.fill_between(r, 0, B3, facecolor=colors[2], alpha=0.3, label='2.22 $\mu$s $<$ t $\leq$ 3.33 $\mu$s',zorder=4) 

ax2.set_ylabel('Brightness (a.u.)', fontsize=14)

ax3 = axs[1,1].twinx() 
ax3.fill_between(r, 0, B1, facecolor=colors[0],alpha=0.1, label='0.00 $\mu$s $<$ t $\leq$ 1.11 $\mu$s',zorder=4) 
ax3.fill_between(r, 0, B2, facecolor=colors[1],alpha=0.2, label='1.11 $\mu$s $<$ t $\leq$ 2.22 $\mu$s',zorder=4) 
ax3.fill_between(r, 0, B3, facecolor=colors[2],alpha=0.3, label='2.22 $\mu$s $<$ t $\leq$ 3.33 $\mu$s',zorder=4) 

ax3.set_ylabel('Brightness (a.u.)', fontsize=14)



axs[0,0].plot(r[:7], T1[:7],'.-', color=colors2[0], label='0.00 $\mu$s $<$ t $\leq$ 1.11 $\mu$s',zorder=1) 
axs[0,0].plot(r[:8], T2[:8],'.-', color=colors2[1], label='1.11 $\mu$s $<$ t $\leq$ 2.22 $\mu$s',zorder=1) 
axs[0,0].plot(r[:6], T3[:6],'.-', color=colors2[2], label='2.22 $\mu$s $<$ t $\leq$ 3.33 $\mu$s',zorder=1) 

r1=np.arange(r[6],5,0.01)
f1 = interp1d([r[6], 5], [T1[6], 0])

r2=np.arange(r[7],5,0.01)
f2 = interp1d([r[7], 5], [T2[7], 0])

r3=np.arange(r[5],5,0.01)
f3 = interp1d([r[5], 5], [T3[5], 0])

line1 = f1(r1)
line2 = f2(r2)
line3 = f3(r3)


axs[0,0].plot(r1, line1,'--', color=colors2[0],alpha=0.5, zorder=1) 
axs[0,0].plot(r2, line2,'--', color=colors2[1],alpha=0.5, zorder=1) 
axs[0,0].plot(r3, line3,'--', color=colors2[2],alpha=0.5, zorder=1) 
axs[0,0].set_ylabel('Temperature (K)', fontsize=14)
axs[0,0].set_xlim(0,6)
axs[0,0].set_ylim(0, 35000)

axs[0,1].plot(r[:7], N1[:7],'.-', color=colors2[0], label='0.00 $\mu$s $<$ t $\leq$ 1.11 $\mu$s',zorder=1) 
axs[0,1].plot(r[:7], N2[:7],'.-', color=colors2[1], label='1.11 $\mu$s $<$ t $\leq$ 2.22 $\mu$s',zorder=1) 
axs[0,1].plot(r[:3], N3[:3],'.-', color=colors2[2], label='2.22 $\mu$s $<$ t $\leq$ 3.33 $\mu$s',zorder=1)       
axs[0,1].set_ylabel('Ne (cm$^{-3}}$)', fontsize=14)
axs[0,1].set_xlim(0,6)
axs[0,1].set_ylim(1e16,3e18)
axs[0,1].set_yscale('log')

axs[1,0].plot(r[:7], C1[:7],'.-', color=colors2[0], label='0.00 $\mu$s $<$ t $\leq$ 1.11 $\mu$s',zorder=1) 
axs[1,0].plot(r[:7], C2[:7],'.-', color=colors2[1], label='1.11 $\mu$s $<$ t $\leq$ 2.22 $\mu$s',zorder=1) 
axs[1,0].plot(r[:3], C3[:3],'.-', color=colors2[2], label='2.22 $\mu$s $<$ t $\leq$ 3.33 $\mu$s',zorder=1)   
axs[1,0].set_xlabel('Radial position (mm)', fontsize=14)
axs[1,0].set_ylabel('Conductivity (S/m$^{-1}$)', fontsize=14)
axs[1,0].set_xlim(0,6)
axs[1,0].set_ylim(0, 20000)

axs[1,1].plot(r[:7], D1[:7],'.-', color=colors2[0], label='0.00 $\mu$s $<$ t $\leq$ 1.11 $\mu$s',zorder=1) 
axs[1,1].plot(r[:7], D2[:7],'.-', color=colors2[1], label='1.11 $\mu$s $<$ t $\leq$ 2.22 $\mu$s',zorder=1) 
axs[1,1].plot(r[:3], D3[:3],'.-', color=colors2[2], label='2.22 $\mu$s $<$ t $\leq$ 3.33 $\mu$s',zorder=1)   
axs[1,1].set_xlabel('Radial position (mm)', fontsize=14)
axs[1,1].set_ylabel('Overpressure', fontsize=14)
axs[1,1].set_xlim(0,6)
axs[1,1].set_ylim(0, 10)

ax0.set_ylim(0,500)
ax1.set_ylim(0,500)
ax2.set_ylim(0,500)
ax3.set_ylim(0,500)  



axs[1,0].legend(loc='upper right')
ax3.legend(loc='upper right')

plt.tight_layout() 
plt.show()


