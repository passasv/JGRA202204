#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 13:42:19 2022
Script to plot figure 7 in article of JGR: Atmospheres, titled:

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

#=======================================================
# Panel 7a
#=======================================================

data = pd.read_csv(dirname + '/Fig7a.csv') 

r = np.asarray(data['Radius(mm)'])
NII_648 = np.asarray(data['NII(648.20nm)'])
OII_648 = np.asarray(data['OII(648.65nm)'])
HI_656  = np.asarray(data['HI(656.27nm)'])
OII_656 = np.asarray(data['OII(656.54nm)'])
N2_660  = np.asarray(data['N2(660.80nm)'])
NII_661 = np.asarray(data['NII(661.05nm)'])
T = np.asarray(data['T(K)'])
B = np.asarray(data['B(a.u.)'])

fig= plt.figure(figsize=(7, 10), dpi=100)
ax0 = plt.subplot(311)
ax1= ax0.twinx()

ax2 = plt.subplot(312)
ax3 = ax2.twinx()

ax4 = plt.subplot(313)
ax5 = ax4.twinx()

ax0.plot(r,NII_648, linestyle='dotted', label ='N II 648.20 nm' )
ax0.plot(r,OII_648,linestyle='dotted', label ='O II 648.65 nm' )
ax0.plot(r,HI_656,linestyle='dashed', label ='H' + r'$\alpha$' + ' 656.27 nm' )
ax0.plot(r,OII_656,linestyle='dotted', label ='O II 656.54 nm' )
ax0.plot(r,N2_660, label ='N$_2$ 660.80 nm' )
ax0.plot(r,NII_661, linestyle='dotted', label ='N II 661.05 nm' ) 
 
ax0.fill_between(r, 0, B/np.amax(B*1.1)*np.amax(HI_656)*1.2, alpha=0.2)

idxt = np.argwhere((B/np.amax(B)*np.amax(N2_660)>10) & (T>0))
ax1.plot(r[idxt],T[idxt], '-k')
ax1.plot([r[idxt[-1]],4.6],[T[idxt[-1]],300], '--k')

idx222=np.argwhere((r>=0) & (r<=16))

nii648 = NII_648[idx222]
oii648 = OII_648[idx222]
h = HI_656[idx222]
oii656 = OII_656[idx222]
n2 = N2_660[idx222]
nii661 = NII_661[idx222]
temp2 = T[idx222]
temp2[-6:]=0

b=B[idx222]

#===============================================================
# Plot box inside main axes

xini = 0
xend = 16
sxini = 4
sxend = 16

facecolor = 'w'
rect = [0.32,0.52,0.625,0.6]
box = ax0.get_position()
width = box.width
height = box.height
inax_position  = ax0.transAxes.transform(rect[0:2])
transFigure = fig.transFigure.inverted()
infig_position = transFigure.transform(inax_position)    
x = infig_position[0]
y = infig_position[1]
width *= rect[2]
height *= rect[3]  
subax = fig.add_axes([x,y,width,height],facecolor=facecolor)  # matplotlib 2.0+

x_labelsize = subax.get_xticklabels()[0].get_size()
y_labelsize = subax.get_yticklabels()[0].get_size()
x_labelsize *= rect[2]**0.5
y_labelsize *= rect[3]**0.5

subax.xaxis.set_tick_params(labelsize=x_labelsize)
subax.yaxis.set_tick_params(labelsize=y_labelsize)

idxz = np.argwhere((r>sxini) & (r<sxend))

subax.plot(r[idxz],NII_648[idxz], linestyle='dotted')
subax.plot(r[idxz],OII_648[idxz], linestyle='dotted')
subax.plot(r[idxz],HI_656[idxz],linestyle='dashed')
subax.plot(r[idxz],OII_656[idxz], linestyle='dotted')
subax.plot(r[idxz],N2_660[idxz])
subax.plot(r[idxz],NII_661[idxz], linestyle='dotted')

subax.patch.set_alpha(0)
subax.spines['right'].set_visible(False)
subax.spines['top'].set_visible(False)
subax.set_xlim(4.6,16)
subax.set_ylim(0,60)


#=======================================================
# Panel 7b
#=======================================================
data = pd.read_csv(dirname + '/Fig7b.csv') 

r = np.asarray(data['Radius(mm)'])
NII_648 = np.asarray(data['NII(648.20nm)'])
OII_648 = np.asarray(data['OII(648.65nm)'])
HI_656  = np.asarray(data['HI(656.27nm)'])
OII_656 = np.asarray(data['OII(656.54nm)'])
N2_660  = np.asarray(data['N2(660.80nm)'])
NII_661 = np.asarray(data['NII(661.05nm)'])
T = np.asarray(data['T(K)'])
B = np.asarray(data['B(a.u.)'])


ax2.plot(r,NII_648, linestyle='dotted', label ='N II 648.20 nm' )
ax2.plot(r,OII_648,linestyle='dotted', label ='O II 648.65 nm' )
ax2.plot(r,HI_656,linestyle='dashed', label ='H' + r'$\alpha$' + ' 656.27 nm' )
ax2.plot(r,OII_656,linestyle='dotted', label ='O II 656.54 nm' )
ax2.plot(r,N2_660, label ='N$_2$ 660.80 nm' )
ax2.plot(r,NII_661, linestyle='dotted', label ='N II 661.05 nm' ) 
 
ax2.fill_between(r, 0, B/np.amax(B*1.1)*np.amax(HI_656)*1.2, alpha=0.2)

idxt = np.argwhere((B/np.amax(B)*np.amax(N2_660)>10) & (T>0))
ax3.plot(r[idxt],T[idxt], '-k')
ax3.plot([r[idxt[-1]],4.6],[T[idxt[-1]],300], '--k')

idx222=np.argwhere((r>=0) & (r<=16))

nii648 = NII_648[idx222]
oii648 = OII_648[idx222]
h = HI_656[idx222]
oii656 = OII_656[idx222]
n2 = N2_660[idx222]
nii661 = NII_661[idx222]
temp2 = T[idx222]
temp2[-6:]=0

b=B[idx222]

#===============================================================
# Plot box inside main axes

facecolor = 'w'
rect = [0.32,0.32,0.625,0.6]
box = ax2.get_position()
width = box.width
height = box.height
inax_position  = ax2.transAxes.transform(rect[0:2])
transFigure = fig.transFigure.inverted()
infig_position = transFigure.transform(inax_position)    
x = infig_position[0]
y = infig_position[1]
width *= rect[2]
height *= rect[3]  
subax = fig.add_axes([x,y,width,height],facecolor=facecolor)  # matplotlib 2.0+

x_labelsize = subax.get_xticklabels()[0].get_size()
y_labelsize = subax.get_yticklabels()[0].get_size()
x_labelsize *= rect[2]**0.5
y_labelsize *= rect[3]**0.5

subax.xaxis.set_tick_params(labelsize=x_labelsize)
subax.yaxis.set_tick_params(labelsize=y_labelsize)

idxz = np.argwhere((r>sxini) & (r<sxend))

subax.plot(r[idxz],NII_648[idxz], linestyle='dotted')
subax.plot(r[idxz],OII_648[idxz], linestyle='dotted')
subax.plot(r[idxz],HI_656[idxz],linestyle='dashed')
subax.plot(r[idxz],OII_656[idxz], linestyle='dotted')
subax.plot(r[idxz],N2_660[idxz])
subax.plot(r[idxz],NII_661[idxz], linestyle='dotted')

subax.patch.set_alpha(0)
subax.spines['right'].set_visible(False)
subax.spines['top'].set_visible(False)
subax.set_xlim(4.6,16)
subax.set_ylim(0,30)


#=======================================================
# Panel 7c
#=======================================================

data = pd.read_csv(dirname + '/Fig7c.csv') 

r = np.asarray(data['Radius(mm)'])
NII_648 = np.asarray(data['NII(648.20nm)'])
OII_648 = np.asarray(data['OII(648.65nm)'])
HI_656  = np.asarray(data['HI(656.27nm)'])
OII_656 = np.asarray(data['OII(656.54nm)'])
N2_660  = np.asarray(data['N2(660.80nm)'])
NII_661 = np.asarray(data['NII(661.05nm)'])
T = np.asarray(data['T(K)'])
B = np.asarray(data['B(a.u.)'])


ax4.plot(r,NII_648, linestyle='dotted', label ='N II 648.20 nm' )
ax4.plot(r,OII_648,linestyle='dotted', label ='O II 648.65 nm' )
ax4.plot(r,HI_656,linestyle='dashed', label ='H' + r'$\alpha$' + ' 656.27 nm' )
ax4.plot(r,OII_656,linestyle='dotted', label ='O II 656.54 nm' )
ax4.plot(r,N2_660, label ='N$_2$ 660.80 nm' )
ax4.plot(r,NII_661, linestyle='dotted', label ='N II 661.05 nm' ) 
 
ax4.fill_between(r, 0, B/np.amax(B*1.1)*np.amax(HI_656)*1.2, alpha=0.2)

idxt = np.argwhere((B/np.amax(B)*np.amax(N2_660)>10) & (T>0))
ax5.plot(r[idxt],T[idxt], '-k')
ax5.plot([r[idxt[-1]],5.2],[T[idxt[-1]],300], '--k')

idx222=np.argwhere((r>=0) & (r<=16))

nii648 = NII_648[idx222]
oii648 = OII_648[idx222]
h = HI_656[idx222]
oii656 = OII_656[idx222]
n2 = N2_660[idx222]
nii661 = NII_661[idx222]
temp2 = T[idx222]
temp2[-6:]=0

b=B[idx222]

#===============================================================
# Plot box inside main axes

facecolor = 'w'
rect = [0.32,0.12,0.625,0.6]
box = ax4.get_position()
width = box.width
height = box.height
inax_position  = ax4.transAxes.transform(rect[0:2])
transFigure = fig.transFigure.inverted()
infig_position = transFigure.transform(inax_position)    
x = infig_position[0]
y = infig_position[1]
width *= rect[2]
height *= rect[3]  
subax = fig.add_axes([x,y,width,height],facecolor=facecolor)  # matplotlib 2.0+

x_labelsize = subax.get_xticklabels()[0].get_size()
y_labelsize = subax.get_yticklabels()[0].get_size()
x_labelsize *= rect[2]**0.5
y_labelsize *= rect[3]**0.5

subax.xaxis.set_tick_params(labelsize=x_labelsize)
subax.yaxis.set_tick_params(labelsize=y_labelsize)

idxz = np.argwhere((r>sxini) & (r<sxend))

subax.plot(r[idxz],NII_648[idxz], linestyle='dotted')
subax.plot(r[idxz],OII_648[idxz], linestyle='dotted')
subax.plot(r[idxz],HI_656[idxz],linestyle='dashed')
subax.plot(r[idxz],OII_656[idxz], linestyle='dotted')
subax.plot(r[idxz],N2_660[idxz])
subax.plot(r[idxz],NII_661[idxz], linestyle='dotted')

subax.patch.set_alpha(0)
subax.spines['right'].set_visible(False)
subax.spines['top'].set_visible(False)
subax.set_xlim(4.6,16)
subax.set_ylim(0,1)
#=======================================================

ylim = ax0.get_ylim()[1]
ax0.text(10.7, ylim*0.05, r'0 $\mu$s $<$ t $\leq$ 1.11 $\mu$s',  fontsize=12,  bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'))    
ax0.legend(loc='upper right')
ax0.axvline(0.5, color='k', linestyle= 'dashed',alpha=0.3)
ax0.axvline(2, color='k', linestyle= 'dashed',alpha=0.3)
ax0.axvline(4.6, color='k', linestyle= 'dashed',alpha=0.3)
ax0.set_ylabel('Brightness (a.u.)', fontsize=14)
ax0.set_xlim(xini,xend)
ax0.set_ylim(0,ylim*1.1)
ax0.text(0.15, ylim*0.93*1.1, r'I',  fontsize=12, color='navy')
ax0.text(1.1, ylim*0.93*1.1, r'II',  fontsize=12, color='navy')
ax0.text(2.8, ylim*0.93*1.1, r'III',  fontsize=12, color='navy')
ax0.text(4.8, ylim*0.93*1.1, r'IV',  fontsize=12, color='navy')


ax1.set_ylabel('Temperature (K)', fontsize=14)
ax1.set_ylim(0,ax1.get_ylim()[1]*1.1)

ylim = ax2.get_ylim()[1] 
ax2.text(9.9, ylim*0.05, r'1.11 $\mu$s $<$ t $\leq$ 2.22 $\mu$s',  fontsize=12,  bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'))
ax2.text(0.15, ylim*0.93*1.1, r'I',  fontsize=12, color='navy')
ax2.text(1.1, ylim*0.93*1.1, r'II',  fontsize=12, color='navy')
ax2.text(2.8, ylim*0.93*1.1, r'III',  fontsize=12, color='navy')
ax2.text(4.8, ylim*0.93*1.1, r'IV',  fontsize=12, color='navy')
ax2.set_ylim(0,ylim*1.1)

ax2.set_ylabel('Brightness (a.u.)', fontsize=14)
ax2.set_xlim(xini,xend)
ax2.axvline(0.5, color='k', linestyle= 'dashed',alpha=0.3)
ax2.axvline(2, color='k', linestyle= 'dashed',alpha=0.3)
ax2.axvline(4.6, color='k', linestyle= 'dashed',alpha=0.3)
ax2.legend(loc='upper right')
ax3.set_ylabel('Temperature (K)', fontsize=14)
ax3.set_ylim(0,ax3.get_ylim()[1]*1.1)

ylim = ax4.get_ylim()[1]
ax4.text(9.9, ylim*0.05, r'2.22 $\mu$s $<$ t $\leq$ 3.33 $\mu$s',  fontsize=12,  bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'))
ax4.axvline(0.5, color='k', linestyle= 'dashed',alpha=0.3)
ax4.axvline(2, color='k', linestyle= 'dashed',alpha=0.3)
ax4.axvline(4.6, color='k', linestyle= 'dashed',alpha=0.3)
ax4.set_ylabel('Brightness (a.u.)', fontsize=14)
ax4.set_xlabel('Radial position (mm)',fontsize=14)
ax4.set_xlim(xini,xend)
ax4.legend(loc='upper right')
ax4.text(0.15, ylim*0.93*1.1, r'I',  fontsize=12, color='navy')
ax4.text(1.1, ylim*0.93*1.1, r'II',  fontsize=12, color='navy')
ax4.text(2.8, ylim*0.93*1.1, r'III',  fontsize=12, color='navy')
ax4.text(4.8, ylim*0.93*1.1, r'IV',  fontsize=12, color='navy')
ax4.set_ylim(0,ylim*1.1)
ax5.set_ylabel('Temperature (K)',fontsize=14)
ax5.set_ylim(0,ax5.get_ylim()[1]*1.1)
plt.tight_layout()

plt.show()



