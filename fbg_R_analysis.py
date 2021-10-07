# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 22:15:55 2021

@author: Anya Houk
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import chirp, find_peaks, peak_widths

df = pd.read_csv("C:/Users/Anya/Desktop/Anritzu OSA/Data/2021-10-06/T_cav1nm_bare.csv")
#dfr = pd.read_csv("C:/Users/Anya/Desktop/Anritzu OSA/Data/2021-09-22/R_33966_bare.csv")

lam = df["Wavelength(nm)"]
power = df["Power(dBm)"]

#lamr = dfr["Wavelength(nm)"]
#powerr = dfr["Power(dBm)"]

#power = (np.array(power)/np.array(powerr))**2
#lam_r = np.linspace(1038.5, 1041.5, 1001)
#power_r = dfr["Power(dBm)"]

mx = np.where(power == max(power))[0][0]
mn = np.where(power == min(power))[0][0]
height = power[mx] - power[mn]
diff = np.abs(np.array(power) - (power[mx] - 3))
print(power[mx], power[mx]/0.707)
l = np.where(diff[:mx] == min(diff[:mx]))[-1][0]
r = np.where(diff[mx:] == min(diff[mx:]))[0][0] + mx
c = int((r - l)/2) + l
#print(round(lam[r]-lam[l], 2), lam[c], r-l, lam[1]-lam[0])

#results_half = peak_widths(power, [c], rel_height=0.2)
#l = int(results_half[2])
#r = int(results_half[3])
#c = int(0.5*(r-l)) + l 

plt.plot(lam, power, color='blue', marker='o', markersize=1, linewidth=1)
#plt.plot(lam_r, power_r, color='yellow', marker='o', markersize=1, linewidth=1)
plt.vlines(lam[c], ymin=min(power)-1, ymax=max(power)+1, color='red', linewidth=1, label=r'$\lambda_{Bragg}$='+str(lam[c]))
plt.vlines(lam[l], ymin=min(power)-1, ymax=max(power)+1, color='green', linewidth=1)
plt.vlines(lam[r], ymin=min(power)-1, ymax=max(power)+1, color='green', linewidth=1)
plt.hlines(power[r], xmin=lam[l], xmax=lam[r], color='green', linewidth=1, label='FWHM='+str(round(lam[r]-lam[l], 2))+'nm')
plt.legend()
plt.title('Transmission')
plt.ylabel('Power (dB)')
plt.xlabel('$\lambda$(nm)')
plt.xlim([lam[0], lam[len(lam)-1]])
plt.ylim([min(power)+10, max(power)])
plt.text(lam[20], max(power)-1, 'Points: '+ str(len(lam))+r', $\Delta\lambda$: '+str(round(lam[1]-lam[0], 4))+'nm')
plt.show()