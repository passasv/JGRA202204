#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 13:42:19 2022
Script to plot figures 3 to 11 in article of JGR: Atmospheres, titled:

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
from PIL import Image

dirname = os.getcwd() + '/Data'
alpha = 0.58                # Spatial dispersion 0.58 mm/px 

#%%===========================================================================
#  FIGURE 3
#=============================================================================

VI = pd.read_csv(dirname + '/Fig3.csv') 

t = VI['Time(s)']
V = VI['V(V)']
I = VI['I(A)']

fig, ax1 = plt.subplots()    
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

plt.show()

#%%===========================================================================
#  FIGURE 4
#=============================================================================

spectrum0 = dirname + '/Fig4a.png'
spectrum1 = dirname + '/Fig4b.png'
spectrum2 = dirname + '/Fig4c.png'

# spectrum0 = r'/home/passasv/Documentos/JGR_Atm_2022/JGR_Atm_RADIAL_042022/Data/Fig4a.png'

img0 = np.asarray(Image.open(spectrum0))
img1 = np.asarray(Image.open(spectrum1))
img2 = np.asarray(Image.open(spectrum2))

rpx= np.arange(-39,52,1)    # Radial positions in px   
           

data = pd.read_csv(dirname + '/Fig4_wavelengths.csv') 

L = data['Wavelength(nm)']

fig, axs = plt.subplots(3,1, figsize=(8,8))
im0 = axs[0].imshow(img0,extent=[L[0],L[127],rpx[0]*alpha,rpx[-1]*alpha], aspect='auto')
im1 = axs[1].imshow(img1,extent=[L[0],L[127],rpx[0]*alpha,rpx[-1]*alpha], aspect='auto')
im2 = axs[2].imshow(img2,extent=[L[0],L[127],rpx[0]*alpha,rpx[-1]*alpha], aspect='auto')

axs[0].tick_params(
    axis='x',          
    which='both',      
    bottom=False,      
    top=False,         
    labelbottom=False) 

axs[1].tick_params(
    axis='x',          
    which='both',      
    bottom=False,      
    top=False,         
    labelbottom=False)

fs=12
axs[2].set_xlabel('Wavelength (nm)', fontsize=fs)
axs[0].set_ylabel('Radial position (mm)', fontsize=fs)
axs[1].set_ylabel('Radial position (mm)', fontsize=fs)
axs[2].set_ylabel('Radial position (mm)', fontsize=fs)


x0 = 645.5
y0 = -19
axs[0].text(x0, y0, '0.72 $\mu$s', bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey', alpha=0.8), alpha=1)
axs[1].text(x0, y0, '1.83 $\mu$s', bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey', alpha=0.8), alpha=1)
axs[2].text(x0, y0, '2.94 $\mu$s', bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey', alpha=0.8), alpha=1)


fig.colorbar(im0, ax=axs[0])
fig.colorbar(im1, ax=axs[1])
fig.colorbar(im2, ax=axs[2])
plt.show()

#%%===========================================================================
#  FIGURE 5
#=============================================================================

arrowprops_k=dict(arrowstyle='->', color='k')

fig, ax = plt.subplots(3,1, gridspec_kw={'height_ratios': [3, 1, 1]})

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
plt.show()

#%%===========================================================================
#  FIGURE 6
#=============================================================================

data = pd.read_csv(dirname + '/Fig6a.csv') 

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
ax1.plot([r[idxt[-1]],r[idxt[-1]+2]],[T[idxt[-1]],300], '--k')

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

axisbg='w'
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

#=======================================================
data = pd.read_csv(dirname + '/Fig6b.csv') 

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
ax3.plot([r[idxt[-1]],4],[T[idxt[-1]],300], '--k')

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
rect = [0.32,0.52,0.625,0.6]
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

#=======================================================

data = pd.read_csv(dirname + '/Fig6c.csv') 

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
ax5.plot([r[idxt[-1]],4],[T[idxt[-1]],300], '--k')

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
rect = [0.32,0.22,0.625,0.6]
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
subax.set_ylim(0,0.1)
#=======================================================

ylim = ax0.get_ylim()[1]
ax0.text(13.7, ylim*0.05, r'0.72 $\mu$s',  fontsize=12,  bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'))    
ax0.legend(loc='upper right')
ax0.axvline(0.5, color='k', linestyle= 'dashed',alpha=0.3)
ax0.axvline(2, color='k', linestyle= 'dashed',alpha=0.3)
ax0.axvline(4, color='k', linestyle= 'dashed',alpha=0.3)
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
ax2.text(13.7, ylim*0.05, r'1.83 $\mu$s',  fontsize=12,  bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'))
ax2.text(0.15, ylim*0.93*1.1, r'I',  fontsize=12, color='navy')
ax2.text(1.1, ylim*0.93*1.1, r'II',  fontsize=12, color='navy')
ax2.text(2.8, ylim*0.93*1.1, r'III',  fontsize=12, color='navy')
ax2.text(4.8, ylim*0.93*1.1, r'IV',  fontsize=12, color='navy')
ax2.set_ylim(0,ylim*1.1)

ax2.set_ylabel('Brightness (a.u.)', fontsize=14)
ax2.set_xlim(xini,xend)
ax2.axvline(0.5, color='k', linestyle= 'dashed',alpha=0.3)
ax2.axvline(2, color='k', linestyle= 'dashed',alpha=0.3)
ax2.axvline(4, color='k', linestyle= 'dashed',alpha=0.3)
ax2.legend(loc='upper right')
ax3.set_ylabel('Temperature (K)', fontsize=14)
ax3.set_ylim(0,ax3.get_ylim()[1]*1.1)

ylim = ax4.get_ylim()[1]
ax4.text(13.7, ylim*0.05, r'2.94 $\mu$s',  fontsize=12,  bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'))
ax4.axvline(0.5, color='k', linestyle= 'dashed',alpha=0.3)
ax4.axvline(2, color='k', linestyle= 'dashed',alpha=0.3)
ax4.axvline(4, color='k', linestyle= 'dashed',alpha=0.3)
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

#%%===========================================================================
#  FIGURE 7
#=============================================================================

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
rect = [0.32,0.52,0.625,0.6]
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
rect = [0.32,0.22,0.625,0.6]
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


#%%===========================================================================
#  FIGURE 8
#=============================================================================

colors = plt.cm.viridis(np.linspace(0,1,5))
colors2 = ['black', 'red', 'orange', 'pink' ,'yellow']

data = pd.read_csv(dirname + '/Fig8a.csv') 

r = np.asarray(data['Radius(mm)'])
T1 = np.asarray(data['T_0.72us(K)'])
T2 = np.asarray(data['T_1.83us(K)'])
T3  = np.asarray(data['T_2.943us(K)'])
B1 = np.asarray(data['B_0.72us(a.u.)'])
B2  = np.asarray(data['B_1.83us(a.u.)'])
B3 = np.asarray(data['B_2.94us(a.u.)'])

data2 = pd.read_csv(dirname + '/Fig8b.csv') 

N1 = np.asarray(data2['Ne_0.72us(K)'])
N2 = np.asarray(data2['Ne_1.83us(K)'])
N3 = np.asarray(data2['Ne_2.94us(K)'])

data3 = pd.read_csv(dirname + '/Fig8c.csv') 

C1 = np.asarray(data3['C_0.72us(K)'])
C2 = np.asarray(data3['C_1.83us(K)'])
C3  = np.asarray(data3['C_2.94us(K)'])

data4 = pd.read_csv(dirname + '/Fig8d.csv') 

D1 = np.asarray(data4['d_0.72us(K)'])
D2 = np.asarray(data4['d_1.83us(K)'])
D3  = np.asarray(data4['d_2.94us(K)'])

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
  

plt.tight_layout() 

axs[1,0].legend(loc='upper right')
ax3.legend(loc='upper right')

plt.show()

#%%===========================================================================
#  FIGURE 9
#=============================================================================

data2 = pd.read_csv(dirname + '/Fig8a.csv') 
B1 = np.asarray(data2['B_0.72us(a.u.)'])
B2  = np.asarray(data2['B_1.83us(a.u.)'])
B3 = np.asarray(data2['B_2.94us(a.u.)'])
r2 = np.asarray(data2['Radius(mm)'])

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
plt.tight_layout()   

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
plt.show()


#%%===========================================================================
#  FIGURE 3
#=============================================================================

colors = plt.cm.viridis(np.linspace(0,1,5))
colors2 = ['black', 'red', 'orange', 'pink' ,'yellow']

data = pd.read_csv(dirname + '/Fig10a.csv') 
alpha=0.58

r = np.asarray(data['Radius(mm)'])
T1 = np.asarray(data['T_1.11us(K)'])
T2 = np.asarray(data['T_2.22us(K)'])
T3  = np.asarray(data['T_3.33us(K)'])
B1 = np.asarray(data['B_1.11us(a.u.)'])
B2  = np.asarray(data['B_2.22us(a.u.)'])
B3 = np.asarray(data['B_3.33us(a.u.)'])

data2 = pd.read_csv(dirname + '/Fig10b.csv') 

N1 = np.asarray(data2['Ne_1.11us(K)'])
N2 = np.asarray(data2['Ne_2.22us(K)'])
N3 = np.asarray(data2['Ne_3.33us(K)'])

data3 = pd.read_csv(dirname + '/Fig10c.csv') 

C1 = np.asarray(data3['C_1.11us(K)'])
C2 = np.asarray(data3['C_2.22us(K)'])
C3  = np.asarray(data3['C_3.33us(K)'])

data4 = pd.read_csv(dirname + '/Fig10d.csv') 

D1 = np.asarray(data4['d_1.11us(K)'])
D2 = np.asarray(data4['d_2.22us(K)'])
D3  = np.asarray(data4['d_3.33us(K)'])

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

plt.tight_layout() 

axs[1,0].legend(loc='upper right')
ax3.legend(loc='upper right')

plt.show()

#%%===========================================================================
#  FIGURE 11
#=============================================================================

data2 = pd.read_csv(dirname + '/Fig10a.csv') 
r2 = np.asarray(data2['Radius(mm)'])
B1 = np.asarray(data2['B_1.11us(a.u.)'])
B2  = np.asarray(data2['B_2.22us(a.u.)'])
B3 = np.asarray(data2['B_3.33us(a.u.)'])

data = pd.read_csv(dirname + '/Fig11a.csv') 

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
ax0.text(0.12, 5e-9, '0 $\mu$s $<$ t $\leq$ 1.11 $\mu$s' , bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'), alpha=1)

ax1.plot(r,T, '--o', color=color, markersize=4, zorder=1) 
ax1.fill_between(r2,0,B1/np.amax(B1)*30000, facecolor=colors[1], alpha=0.1)

ax1.set_xlim(0,6)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylabel('Temperature (K)', color=color, fontsize=14)
ax1.set_ylim(0,35000)

data = pd.read_csv(dirname + '/Fig11b.csv') 

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

rpx=np.arange(6,16)
ax2.fill_between(rpx*alpha, 1e-10, 1e20, facecolor=color3, hatch='//', alpha=alphaval, zorder=2)

ax2.text(0.12, 5e-9, '1.11 $\mu$s $<$ t $\leq$ 2.22 $\mu$s' , bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'), alpha=1)
ax2.text(0.1, 1e15, r'(b)',  fontsize=12)

ax3.plot(r,T, '--o', markersize=4,color=color, zorder=1) 
ax3.fill_between(r2,0,B2/np.amax(B1)*30000, facecolor=colors[1], alpha=0.1)
ax3.tick_params(axis='y', labelcolor=color)
ax3.set_ylabel('Temperature (K)', color=color, fontsize=14)
ax3.set_ylim(0,35000)

data = pd.read_csv(dirname + '/Fig11c.csv') 


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
ax4.text(0.12, 5e-9, '2.22 $\mu$s $<$ t $\leq$ 3.33 $\mu$s' , bbox=dict(boxstyle="square,pad=0.3", facecolor='white', ec='grey'), alpha=1)
ax4.text(0.1, 1e15, r'(c)',  fontsize=12)

rpx=np.arange(2,16)
ax4.fill_between(rpx*alpha, 1e-10, 1e20, facecolor=color3, hatch='//', alpha=alphaval, zorder=2)


ax5.plot(r,T, '--o', markersize=4,color=color, zorder=1) 
ax5.fill_between(r2,0,B3/np.amax(B1)*30000, facecolor=colors[1], alpha=0.1)
ax5.tick_params(axis='y', labelcolor=color)
ax5.set_ylabel('Temperature (K)', color=color, fontsize=14)
ax5.set_ylim(0,35000)

plt.rcParams.update({'hatch.color': 'lightgrey'})
plt.tight_layout()   

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
plt.show()


