#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 13:42:19 2022
Script to plot figure 4 in article of JGR: Atmospheres, titled:

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
from PIL import Image

dirname = os.getcwd() + '/Data'


spectrum0 = dirname + '/Fig4a.png'
spectrum1 = dirname + '/Fig4b.png'
spectrum2 = dirname + '/Fig4c.png'

img0 = np.asarray(Image.open(spectrum0))
img1 = np.asarray(Image.open(spectrum1))
img2 = np.asarray(Image.open(spectrum2))

rpx= np.arange(-39,52,1)    # Radial positions in px   
alpha = 0.58                # Spatial dispersion 0.58 mm/px           

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

plt.tight_layout()
plt.show()
