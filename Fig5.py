#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 13:42:19 2022
Script to plot figure 5 in article of JGR: Atmospheres, titled:

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

arrowprops_k=dict(arrowstyle='->', color='k')

fig, ax = plt.subplots(3,1, gridspec_kw={'height_ratios': [3, 1, 1]}, figsize=(15,12))

data = pd.read_csv(dirname + '/Fig5a.csv') 

L = data['Wavelength(nm)']
a1= data['0.00_mm']
b1= data['1.16_mm']
c1= data['2.32_mm']
d1= data['3.48_mm']

ax[0].plot(L,a1, label='0.00 mm')
ax[0].plot(L,b1, label='1.16 mm')
ax[0].plot(L,c1, label='2.32 mm')
ax[0].plot(L,d1, label='3.48 mm')

data = pd.read_csv(dirname + '/Fig5b.csv') 

L = data['Wavelength(nm)']
a2= data['0.00_mm']
b2= data['1.16_mm']
c2= data['2.32_mm']
d2= data['3.48_mm']

ax[1].plot(L,a2, label='0.00 mm')
ax[1].plot(L,b2, label='1.16 mm')
ax[1].plot(L,c2, label='2.32 mm')
ax[1].plot(L,d2, label='3.48 mm')

data = pd.read_csv(dirname + '/Fig5c.csv') 

L = data['Wavelength(nm)']
a3= data['0.00_mm']
b3= data['1.16_mm']
c3= data['2.32_mm']
d3= data['3.48_mm']

ax[2].plot(L,a3, label='0.00 mm')
ax[2].plot(L,b3, label='1.16 mm')
ax[2].plot(L,c3, label='2.32 mm')
ax[2].plot(L,d3, label='3.48 mm')


j=0
    
ax[0].annotate('N III 647.88 nm', xy=(647.876, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('N II 648.20 nm', xy=(648.20, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k, fontweight="bold")

ax[0].annotate('O$_2⁺$ 648.9 nm + Cu II 648.88 nm', xy=(648.9, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('O II 648.65 nm', xy=(648.65, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k, fontweight="bold")

ax[0].annotate('N$_2$ 650.0 nm', xy=(650.0, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('Cu II 649.40 nm', xy=(649.40, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('Cu II 652.38 nm', xy=(652.38, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)
 
ax[0].annotate('N$_2$ FPS (7, 4) 653.00 nm', xy=(653.00, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('N$_2$ 655.5 nm', xy=(655.5, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('O II 655.90 nm', xy=(655.90, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
               rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('Cu II 654.16 nm', xy=(654.16, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
               rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('H' + r'$_\alpha$' + ' 656.27 nm', xy=(656.27, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k, fontweight="bold")

ax[0].annotate('O II 656.539 nm', xy=(656.539, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k, fontweight="bold")

ax[0].annotate('Cu II 656.79 nm + C I 656.87 nm' , xy=(656.84, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
               rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('O II 657.11 nm', xy=(657.11, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('N$_2$ + O$_2$$^+$ 657.7 nm', xy=(657.7, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('C I 658.76 nm', xy=(658.76, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
                rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('Cu II 659.29 nm', xy=(659.29, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
                rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('C I 659.52 nm + N II 659.567 nm', xy=(659.52, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
                rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('Cu I 659.96 nm', xy=(659.96, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
                rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('N I 660.30 nm + Cu II 660.35 nm' , xy=(660.35, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k)

ax[0].annotate('N$_2$ FPS (6, 3) 660.80 nm', xy=(660.80, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k, fontweight="bold")

ax[0].annotate('N II 661.05 nm', xy=(661.05, np.amax(a1)), xytext=(0, 25), textcoords='offset points', color = 'k',
              rotation=90, va='bottom', ha='center', annotation_clip=False, arrowprops=arrowprops_k, fontweight="bold")
 

for j in range(3):
    if j==0:
        a = a1
    elif j == 1: 
        a = a2
    else:
        a = a3

    ax[j].vlines(647.88, 0, np.amax(a), 'grey', 'dashed', alpha=0.5 )
    ax[j].vlines(648.20, 0, np.amax(a), 'grey', 'dashed' )
    ax[j].vlines(648.65, 0, np.amax(a), 'grey', 'dashed' )
    ax[j].vlines(648.9, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5) # O2
    ax[j].vlines(649.40, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5)
    ax[j].vlines(650.0, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5)
    ax[j].vlines(652.38, 0, np.amax(a), 'grey', 'dashed', alpha=0.5 )
    ax[j].vlines(653.00, 0, np.amax(a), 'grey', 'dashed', alpha=0.5 )
    ax[j].vlines(654.16, 0, np.amax(a), 'grey', 'dashed', alpha=0.5 ) # Cu II
    ax[j].vlines(655.5, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5) # N2   
    ax[j].vlines(655.90, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5)
    ax[j].vlines(656.27, 0, np.amax(a), 'grey', 'dashed' )
    ax[j].vlines(656.84, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5) # doblete C I y Cu II
    ax[j].vlines(657.11, 0, np.amax(a), 'grey', 'dashed', alpha=0.5)
    ax[j].vlines(656.539, 0, np.amax(a), 'grey', 'dashed' )
    ax[j].vlines(657.7, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5)
    ax[j].vlines(658.76, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5) # C I
    ax[j].vlines(659.29, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5) # Cu II
    ax[j].vlines(659.52, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5) # C I
    ax[j].vlines(659.96, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5) # Cu II
    ax[j].vlines(660.35, 0, np.amax(a), 'grey', 'dashed' , alpha=0.5) # Cu II
    ax[j].vlines(660.80, 0, np.amax(a), 'grey', 'dashed' )
    ax[j].vlines(661.05, 0, np.amax(a), 'grey', 'dashed' )
         
    
    ax[j].set_ylabel('Brightness (a.u.)', fontsize=18)


ax[0].legend()
    
ax[0].set_ylim(-10,np.amax(a1)*2.5)
ax[1].set_ylim(-10,np.amax(a2)*1.1)
ax[2].set_ylim(-10,np.amax(a3)*1.1)

    
ax[0].text(645.1, np.amax(a1)*0.91*2.61, '0.72 $\mu$s',  fontsize=12, bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'))
ax[1].text(645.1, np.amax(a2)*0.92, '1.83 $\mu$s',  fontsize=12, bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'))
ax[2].text(645.1, np.amax(a3)*0.92, '2.94 $\mu$s',  fontsize=12, bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'))

ax[2].set_xlabel('Wavelength (nm)', fontsize=18)

plt.sca(ax[2])
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.sca(ax[0])
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.sca(ax[1])
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.tight_layout()
plt.show()



