# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 21:18:28 2021

@author: Anya

"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import chirp, find_peaks, peak_widths

def toWatt(power):
    watt = []
    for p in power:
        watt.append(10**((p-30)/10))
    return np.array(watt)

df1 = pd.read_csv("C:/Users/Anya/Desktop/Anritzu OSA/Data/2021-09-22/Loss/R_source_1.csv")
df2 = pd.read_csv("C:/Users/Anya/Desktop/Anritzu OSA/Data/2021-09-22/Loss/R_source_2.csv")
df3 = pd.read_csv("C:/Users/Anya/Desktop/Anritzu OSA/Data/2021-09-22/Loss/R_source_3.csv")

dft1 = pd.read_csv("C:/Users/Anya/Desktop/Anritzu OSA/Data/2021-09-22/Loss/T_source_1.csv")
dft2 = pd.read_csv("C:/Users/Anya/Desktop/Anritzu OSA/Data/2021-09-22/Loss/T_source_2.csv")

lam = df1["Wavelength(nm)"]
power1 = df1["Power(dBm)"]
#lam2 = df2["Wavelength(nm)"]
power2 = df2["Power(dBm)"]
#lam3 = df3["Wavelength(nm)"]
power3 = df3["Power(dBm)"]

lamt = df1["Wavelength(nm)"]
powert1 = dft1["Power(dBm)"]
powert2 = dft2["Power(dBm)"]
losst12 = powert2 - powert1

loss12 = power1 - power2
loss23 = power2 - power3
print(np.average(loss12), np.average(loss23), np.average(losst12))
#plt.plot(lamt, powert1, color='blue', marker='o', markersize=1, linewidth=1)
#plt.plot(lamt, powert2, color='blue', marker='o', markersize=1, linewidth=1)

plt.plot(lam, loss12, color='blue', marker='o', markersize=1, linewidth=1)
plt.plot(lam, loss23, color='green', marker='o', markersize=1, linewidth=1)
plt.plot(lamt, losst12, color='red', marker='o', markersize=1, linewidth=1)

plt.show()