#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 13:42:19 2022
Script to plot figure 8 in article of JGR: Atmospheres, titled:

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
# Panel 8a
#=======================================================

data = pd.read_csv(dirname + '/Fig8a.csv') 

r = np.asarray(data['Radius(mm)'])
T1 = np.asarray(data['T_0.72us(K)'])
T2 = np.asarray(data['T_1.83us(K)'])
T3  = np.asarray(data['T_2.943us(K)'])
B1 = np.asarray(data['B_0.72us(a.u.)'])
B2  = np.asarray(data['B_1.83us(a.u.)'])
B3 = np.asarray(data['B_2.94us(a.u.)'])

#=======================================================
# Panel 8b
#=======================================================

data2 = pd.read_csv(dirname + '/Fig8b.csv') 

N1 = np.asarray(data2['Ne_0.72us(K)'])
N2 = np.asarray(data2['Ne_1.83us(K)'])
N3 = np.asarray(data2['Ne_2.94us(K)'])

#=======================================================
# Panel 8c
#=======================================================

data3 = pd.read_csv(dirname + '/Fig8c.csv') 

C1 = np.asarray(data3['C_0.72us(K)'])
C2 = np.asarray(data3['C_1.83us(K)'])
C3  = np.asarray(data3['C_2.94us(K)'])

#=======================================================
# Panel 8d
#=======================================================

data4 = pd.read_csv(dirname + '/Fig8d.csv') 

D1 = np.asarray(data4['d_0.72us(K)'])
D2 = np.asarray(data4['d_1.83us(K)'])
D3  = np.asarray(data4['d_2.94us(K)'])

#=======================================================
# Plots
#=======================================================

fig, axs = plt.subplots(2,2, figsize=(9,8))

ax0 = axs[0,0].twinx() 
ax0.fill_between(r, 0, B1, facecolor=colors[0],alpha=0.1, label='0.72 $\mu$s',zorder=4) 
ax0.fill_between(r, 0, B2, facecolor=colors[1],alpha=0.2, label='1.83 $\mu$s',zorder=4) 
ax0.fill_between(r, 0, B3, facecolor=colors[2],alpha=0.3, label='2.94 $\mu$s',zorder=4) 
ax0.set_ylim(0,np.amax([np.amax(B1),np.amax(B2),np.amax(B3)]))

ax0.text(0.3, 0.9*np.amax([np.amax(B1),np.amax(B2),np.amax(B3)]), r'(a)',  fontsize=12)
ax0.set_ylabel('Brightness (a.u.)', fontsize=14)

ax1 = axs[0,1].twinx() 
ax1.fill_between(r, 0, B1, facecolor=colors[0], alpha=0.1, label='0.72 $\mu$s',zorder=4) 
ax1.fill_between(r, 0, B2, facecolor=colors[1], alpha=0.2, label='1.83 $\mu$s',zorder=4) 
ax1.fill_between(r, 0, B3, facecolor=colors[2], alpha=0.3, label='2.94 $\mu$s',zorder=4) 
ax1.set_ylim(0,np.amax([np.amax(B1),np.amax(B2),np.amax(B3)]))

ax1.text(0.3, 0.9*np.amax([np.amax(B1),np.amax(B2),np.amax(B3)]), r'(b)',  fontsize=12)
ax1.set_ylabel('Brightness (a.u.)', fontsize=14)

ax2 = axs[1,0].twinx() 
ax2.fill_between(r, 0, B1, facecolor=colors[0], alpha=0.1, label='0.72 $\mu$s',zorder=4) 
ax2.fill_between(r, 0, B2, facecolor=colors[1], alpha=0.2, label='1.83 $\mu$s',zorder=4) 
ax2.fill_between(r, 0, B3, facecolor=colors[2], alpha=0.3, label='2.94 $\mu$s',zorder=4) 
ax2.set_ylim(0,np.amax([np.amax(B1),np.amax(B2),np.amax(B3)]))

ax2.text(0.3, 0.9*np.amax([np.amax(B1),np.amax(B2),np.amax(B3)]), r'(c)',  fontsize=12)
ax2.set_ylabel('Brightness (a.u.)', fontsize=14)

ax3 = axs[1,1].twinx() 
ax3.fill_between(r, 0, B1, facecolor=colors[0],alpha=0.1, label='0.72 $\mu$s',zorder=4) 
ax3.fill_between(r, 0, B2, facecolor=colors[1],alpha=0.2, label='1.83 $\mu$s',zorder=4) 
ax3.fill_between(r, 0, B3, facecolor=colors[2],alpha=0.3, label='2.94 $\mu$s',zorder=4) 
ax3.set_ylim(0,np.amax([np.amax(B1),np.amax(B2),np.amax(B3)]))

ax3.text(0.3, 0.9*np.amax([np.amax(B1),np.amax(B2),np.amax(B3)]), r'(d)',  fontsize=12)
ax3.set_ylabel('Brightness (a.u.)', fontsize=14)


                           
axs[0,0].plot(r[:7], T1[:7],'.-', color=colors2[0], label='0.72 $\mu$s',zorder=1) 
axs[0,0].plot(r[:8], T2[:8],'.-', color=colors2[1], label='1.83 $\mu$s',zorder=1) 
axs[0,0].plot(r[:4], T3[:4],'.-', color=colors2[2], label='2.94 $\mu$s',zorder=1) 

r1=np.arange(r[6],5,0.01)
f1 = interp1d([r[6], 5], [T1[6], 0])

r2=np.arange(r[7],5,0.01)
f2 = interp1d([r[7], 5], [T2[7], 0])

r3=np.arange(r[3],5,0.01)
f3 = interp1d([r[3], 5], [T3[3], 0])

line1 = f1(r1)
line2 = f2(r2)
line3 = f3(r3)


axs[0,0].plot(r1, line1,'--', color=colors2[0],alpha=0.5, zorder=1) 
axs[0,0].plot(r2, line2,'--', color=colors2[1],alpha=0.5, zorder=1) 
axs[0,0].plot(r3, line3,'--', color=colors2[2],alpha=0.5, zorder=1) 

axs[0,0].set_ylabel('Temperature (K)', fontsize=14)
axs[0,0].set_xlim(0,6)
axs[0,0].set_ylim(0, 35000)

axs[0,1].plot(r[:7], N1[:7],'.-', color=colors2[0], label='0.72 $\mu$s',zorder=1) 
axs[0,1].plot(r[:5], N2[:5],'.-', color=colors2[1], label='1.83 $\mu$s',zorder=1) 
axs[0,1].plot(r[:3], N3[:3],'.-', color=colors2[2], label='2.94 $\mu$s',zorder=1)       
axs[0,1].set_ylabel('Ne (cm$^{-3}}$)', fontsize=14)
axs[0,1].set_xlim(0,6)
axs[0,1].set_ylim(1e16,3e18)
axs[0,1].set_yscale('log')

axs[1,0].plot(r[:7], C1[:7],'.-', color=colors2[0], label='0.72 $\mu$s',zorder=1) 
axs[1,0].plot(r[:5], C2[:5],'.-', color=colors2[1], label='1.83 $\mu$s',zorder=1) 
axs[1,0].plot(r[:3], C3[:3],'.-', color=colors2[2], label='2.94 $\mu$s',zorder=1)   
axs[1,0].set_xlabel('Radial position (mm)', fontsize=14)
axs[1,0].set_ylabel('Conductivity (S/m$^{-1}$)', fontsize=14)
axs[1,0].set_xlim(0,6)
axs[1,0].set_ylim(0, 20000)

axs[1,1].plot(r[:7], D1[:7],'.-', color=colors2[0], label='0.72 $\mu$s',zorder=1) 
axs[1,1].plot(r[:5], D2[:5],'.-', color=colors2[1], label='1.83 $\mu$s',zorder=1) 
axs[1,1].plot(r[:3], D3[:3],'.-', color=colors2[2], label='2.94 $\mu$s',zorder=1) 
axs[1,1].set_xlabel('Radial position (mm)', fontsize=14)
axs[1,1].set_ylabel('Overpressure', fontsize=14)
axs[1,1].set_xlim(0,6)
axs[1,1].set_ylim(0, 10)
        
plt.sca(axs[0,0])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.sca(axs[0,1])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.sca(axs[1,0])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.sca(axs[1,1])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.sca(ax0)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.sca(ax1)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.sca(ax2)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.sca(ax3)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
  
axs[1,0].legend(loc='upper right')
ax3.legend(loc='upper right')

plt.tight_layout()
plt.show()
