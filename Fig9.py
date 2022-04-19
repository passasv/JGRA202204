#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 13:42:19 2022
Script to plot figure 9 in article of JGR: Atmospheres, titled:

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


data2 = pd.read_csv(dirname + '/Fig8a.csv') 
B1 = np.asarray(data2['B_0.72us(a.u.)']) # Brightness
B2  = np.asarray(data2['B_1.83us(a.u.)'])
B3 = np.asarray(data2['B_2.94us(a.u.)'])
r2 = np.asarray(data2['Radius(mm)'])

#=======================================================
# Panel 9a
#=======================================================

data = pd.read_csv(dirname + '/Fig9a.csv') 


r = np.asarray(data['Radius(mm)'])
Ne = np.asarray(data['Ne(cm^-3)'])
Neq = np.asarray(data['Neq(cm^-3)'])
N2 = np.asarray(data['N2(cm^-3)'])
NO = np.asarray(data['NO(cm^-3)'])
O2 = np.asarray(data['O2(cm^-3)'])
OH = np.asarray(data['OH(cm^-3)'])
H2 = np.asarray(data['H2(cm^-3)'])
N2O = np.asarray(data['N2O(cm^-3)'])
NO2 = np.asarray(data['NO2(cm^-3)'])
HO2 = np.asarray(data['HO2(cm^-3)'])
O3 = np.asarray(data['O3(cm^-3)'])
H2O = np.asarray(data['H2O(cm^-3)'])
T = np.asarray(data['T(K)'])


colors = plt.cm.viridis(np.linspace(0,1,5))
color = 'black'
color3= 'grey'
alphaval = 0.1
fss=9

plt.figure(figsize=(7, 10), dpi=100)
ax0 = plt.subplot(311)
ax1= ax0.twinx()

ax2 = plt.subplot(312)
ax3 = ax2.twinx()

ax4 = plt.subplot(313)
ax5 = ax4.twinx()

ax0.plot(r,Ne,'-.', zorder=2, label='N$_{e}$') 
ax0.plot(r,Neq,'-.', zorder=2, label='N$_{e}^{EQ}$')  
ax0.plot(r,N2,zorder=2, label='N$_2$')     
ax0.plot(r,NO,zorder=2, label='NO')
ax0.plot(r,O2,zorder=2, label='O$_2$')
ax0.plot(r,OH,zorder=2, label='OH')
ax0.plot(r,H2,zorder=2, label='H$_2$') 
ax0.plot(r,N2O,zorder=2, label='N$_2$O') 
ax0.plot(r,NO2,zorder=2, label='NO$_2$') 
ax0.plot(r,HO2,zorder=2, label='HO$_2$') 
ax0.plot(r,O3,zorder=2, label='O$_3$') 
ax0.plot(r,H2O,zorder=2, label='H$_2$O') 


rpx=np.arange(6,16)
ax0.fill_between(rpx*alpha, 1e-10, 1e20, facecolor=color3, hatch='//', alpha=alphaval, zorder=2)
ax0.set_yscale('log') 
ax0.set_xlim(0,6)
ax0.set_ylim((1e-10,1e20))
ax0.legend(loc='right', framealpha=1, fontsize=fss) 
ax0.set_ylabel('Concentration (cm$^{-3}$)', fontsize=14)
ax0.text(0.1, 1e15, r'(a)',  fontsize=12)
ax0.text(0.12, 5e-9, '0.72 $\mu$s' , bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'), alpha=1)

ax1.plot(r,T, '--o', color=color, markersize=4, zorder=1) 
ax1.fill_between(r2,0,B1/np.amax(B1)*30000, facecolor=colors[1], alpha=0.1)
ax1.set_xlim(0,6)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylabel('Temperature (K)', color=color, fontsize=14)
ax1.set_ylim(0,35000)

#=======================================================
# Panel 9b
#=======================================================


data = pd.read_csv(dirname + '/Fig9b.csv') 

r = np.asarray(data['Radius(mm)'])
Ne = np.asarray(data['Ne(cm^-3)'])
Neq = np.asarray(data['Neq(cm^-3)'])
N2 = np.asarray(data['N2(cm^-3)'])
NO = np.asarray(data['NO(cm^-3)'])
O2 = np.asarray(data['O2(cm^-3)'])
OH = np.asarray(data['OH(cm^-3)'])
H2 = np.asarray(data['H2(cm^-3)'])
N2O = np.asarray(data['N2O(cm^-3)'])
NO2 = np.asarray(data['NO2(cm^-3)'])
HO2 = np.asarray(data['HO2(cm^-3)'])
O3 = np.asarray(data['O3(cm^-3)'])
H2O = np.asarray(data['H2O(cm^-3)'])
T = np.asarray(data['T(K)'])

ax2.plot(r,Ne, '-.', zorder=2, label='N$_{e}$') 
ax2.plot(r,Neq,'-.', zorder=2, label='N$_{e}^{EQ}$')  
ax2.plot(r,N2,zorder=2, label='N$_2$')     
ax2.plot(r,NO,zorder=2, label='NO')
ax2.plot(r,O2,zorder=2, label='O$_2$')
ax2.plot(r,OH,zorder=2, label='OH')
ax2.plot(r,H2,zorder=2, label='H$_2$') 
ax2.plot(r,N2O,zorder=2, label='N$_2$O') 
ax2.plot(r,NO2,zorder=2, label='NO$_2$') 
ax2.plot(r,HO2,zorder=2, label='HO$_2$') 
ax2.plot(r,O3,zorder=2, label='O$_3$') 
ax2.plot(r,H2O,zorder=2, label='H$_2$O') 
ax2.set_yscale('log') 

ax2.legend(loc='right', framealpha=1, fontsize=fss) 
ax2.set_xlim(0,6)
ax2.set_ylabel('Concentration (cm$^{-3}$)', fontsize=14)
ax2.set_ylim((1e-10,1e20))
rpx=np.arange(4,16)
ax2.fill_between(rpx*alpha, 1e-10, 1e20, facecolor=color3, hatch='//', alpha=alphaval, zorder=2)
ax2.text(0.12, 5e-9, '1.83 $\mu$s' , bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'), alpha=1)
ax2.text(0.1, 1e15, r'(b)',  fontsize=12)
ax3.plot(r,T, '--o', markersize=4,color=color, zorder=1) 
ax3.fill_between(r2,0,B2/np.amax(B1)*30000, facecolor=colors[1], alpha=0.1)
ax3.tick_params(axis='y', labelcolor=color)
ax3.set_ylabel('Temperature (K)', color=color, fontsize=14)
ax3.set_ylim(0,35000)

#=======================================================
# Panel 9c
#=======================================================

data = pd.read_csv(dirname + '/Fig9c.csv') 

r = np.asarray(data['Radius(mm)'])
Ne = np.asarray(data['Ne(cm^-3)'])
Neq = np.asarray(data['Neq(cm^-3)'])
N2 = np.asarray(data['N2(cm^-3)'])
NO = np.asarray(data['NO(cm^-3)'])
O2 = np.asarray(data['O2(cm^-3)'])
OH = np.asarray(data['OH(cm^-3)'])
H2 = np.asarray(data['H2(cm^-3)'])
N2O = np.asarray(data['N2O(cm^-3)'])
NO2 = np.asarray(data['NO2(cm^-3)'])
HO2 = np.asarray(data['HO2(cm^-3)'])
O3 = np.asarray(data['O3(cm^-3)'])
H2O = np.asarray(data['H2O(cm^-3)'])
T = np.asarray(data['T(K)'])

ax4.plot(r,Ne, '-.', zorder=2, label='N$_{e}$') 
ax4.plot(r,Neq,'-.', zorder=2, label='N$_{e}^{EQ}$')  
ax4.plot(r,N2,zorder=2, label='N$_2$')     
ax4.plot(r,NO,zorder=2, label='NO')
ax4.plot(r,O2,zorder=2, label='O$_2$')
ax4.plot(r,OH,zorder=2, label='OH')
ax4.plot(r,H2,zorder=2, label='H$_2$') 
ax4.plot(r,N2O,zorder=2, label='N$_2$O') 
ax4.plot(r,NO2,zorder=2, label='NO$_2$') 
ax4.plot(r,HO2,zorder=2, label='HO$_2$') 
ax4.plot(r,O3,zorder=2, label='O$_3$') 
ax4.plot(r,H2O,zorder=2, label='H$_2$O') 
ax4.set_yscale('log') 

ax4.legend(loc='right', framealpha=1, fontsize=fss) 
ax4.set_xlim(0,6)
ax4.set_ylabel('Concentration (cm$^{-3}$)', fontsize=14)
ax4.set_ylim((1e-10,1e20))
rpx=np.arange(2,16)
ax4.fill_between(rpx*alpha, 1e-10, 1e20, facecolor=color3, hatch='//', alpha=alphaval, zorder=2)
ax4.text(0.12, 5e-9, '2.94 $\mu$s' , bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'), alpha=1)
ax4.text(0.1, 1e15, r'(c)',  fontsize=12)
ax5.plot(r,T, '--o', markersize=4,color=color, zorder=1) 
ax5.fill_between(r2,0,B3/np.amax(B1)*30000, facecolor=colors[1], alpha=0.1)
ax5.tick_params(axis='y', labelcolor=color)
ax5.set_ylabel('Temperature (K)', color=color, fontsize=14)
ax5.set_ylim(0,35000)


plt.rcParams.update({'hatch.color': 'lightgrey'})
  

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

plt.sca(ax4)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.sca(ax5)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.tight_layout() 
plt.show()


