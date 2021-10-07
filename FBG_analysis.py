# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 22:15:55 2021

@author: Anya Houk
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from scipy.signal import savgol_filter

def toWatt(power):
    watt = []
    for p in power:
        watt.append(10**((p-30)/10))
    return np.array(watt)

def downsize(arr, new_len):
    res = int(np.ceil(len(arr)/new_len))
    print(res)
    new_arr = []
    for i in range(0, len(arr), res):
        new_arr.append(np.mean(arr[i:i+res]))
    
    print(len(new_arr))
    return new_arr
        
def find_fwhm(power):
    mn = np.where(power == min(power))[0][0]
    mx_l = np.where(power[:mn] == max(power[:mn]))[0][0]
    mx_r = np.where(power[mn:] == max(power[mn:]))[0][0] + mn
    mx = mx_l if power[mx_l] < power[mx_r] else mx_r
    if(power[mx] <= 1):
        fwhm = power[mn] + (power[mx] - power[mn])/2
    else:
        fwhm = power[mn] + (1-power[mn])/2
    #diff = np.abs(np.array(power) - ((1-power[mn])/2))
    #l = np.where(diff[:mn] == min(diff[:mn]))[-1][-1]
    #r = np.where(diff[mn:] == min(diff[mn:]))[-1][-1] + mn
    l = np.where(power[:mn] >= fwhm)[0][-1]
    r = np.where(power[mn:] >= fwhm)[0][0] + mn
    c = int((r - l)/2) + l    
    
    return l, r, c

def find_fwhm_power(power):
    power = savgol_filter(power, 5, 2)
    mn = np.where(power == min(power))[0][0]
    mx_l = np.where(power[:mn] == max(power[:mn]))[0][0]
    mx_r = np.where(power[mn:] == max(power[mn:]))[0][0] + mn
    mx = mx_l if power[mx_l] < power[mx_r] else mx_r
    fwhm = (power[mx] - 3 )
    
    l = np.where(power[:mn] >= fwhm)[0][-1]
    r = np.where(power[mn:] >= fwhm)[0][0] + mn
    c = int((r - l)/2) + l    
    
    return l, r, c


def plot_spectrum(lam, power, percent, path, sn):
    #l, r, c = find_fwhm(percent)
    ll, rr, cc = find_fwhm_power(power)
    fig, ax = plt.subplots(2,1)
    
    ax[0].plot(lam, input_power, color='orange', marker='o', markersize=1, linewidth=1, label='Input Power')
    ax[0].plot(lam, power, color='blue', marker='o', markersize=1, linewidth=1, label='Transmitted Power')
    ax[0].vlines(lam[cc], ymin=min(power), ymax=max(power), color='red', linewidth=1, label=r'$\lambda_{Bragg}$='+str(lam[cc]))
    ax[0].vlines(lam[ll], ymin=min(power), ymax=max(power), color='green', linewidth=1)
    ax[0].vlines(lam[rr], ymin=min(power), ymax=max(power), color='green', linewidth=1)
    ax[0].hlines(power[rr], xmin=lam[ll], xmax=lam[rr], color='green', linewidth=1, label='FWHM='+str(round(lam[rr]-lam[ll], 2))+'nm')
    ax[0].text(lam[20], min(power), 'Points: '+ str(len(lam))+r', $\Delta\lambda$: '+str(round(lam[1]-lam[0], 4))+'nm')
    
    percent *=100
    ax[1].plot(lam, percent, color='blue', marker='o', markersize=1, linewidth=1, label='R='+str(round(100-min(percent), 4))+'%')
    #ax[1].vlines(lam[c], ymin=min(percent), ymax=max(percent), color='red', linewidth=1, label=r'$\lambda_{Bragg}$='+str(lam[c]))
    #ax[1].vlines(lam[l], ymin=min(percent), ymax=max(percent), color='green', linewidth=1)
    #ax[1].vlines(lam[r], ymin=min(percent), ymax=max(percent), color='green', linewidth=1)
    #ax[1].hlines(percent[r], xmin=lam[l], xmax=lam[r], color='green', linewidth=1, label='FWHM='+str(round(lam[r]-lam[l], 2))+'nm')
    ax[1].hlines(min(percent), xmin=min(lam), xmax=max(lam), color='green', linewidth=1, label='Tmin='+str(round(min(percent), 4))+'%')
    ax[1].hlines(max(percent), xmin=min(lam), xmax=max(lam), color='green', linewidth=1, label='T_max='+str(round(max(percent), 2))+'%')

    
    ax[0].set_title('Transmission SN#' + sn)
    ax[0].set_ylabel('Power (dBm)')
    ax[0].set_xlabel('$\lambda$(nm)')
    ax[1].set_ylabel('% Transmission')
    ax[1].set_xlabel('$\lambda$(nm)')
    
    ax[0].set_xlim([lam[0], lam[len(lam)-1]])
    #ax[0].set_ylim([min(power)-1, max(input_power)+2])
    ax[1].set_xlim([lam[0], lam[len(lam)-1]])
    ax[1].set_ylim([-5, 105])
    ax[0].legend(loc='lower left', bbox_to_anchor=(0.8, 0),
              fancybox=True, ncol=1)
    ax[1].legend(loc='lower left', bbox_to_anchor=(0.8, 0),
              fancybox=True, ncol=1)
    
    fig.savefig(path+"plots/T_"+sn+'.png', dpi=fig.dpi, bbox_inches='tight', pad_inches=0.5)


path = "C:/Users/Anya/Desktop/Anritzu OSA/Data/2021-10-06/"
files = glob.glob(path + '*.csv')
fbg_num = list(set([f.split('_')[1].split('.')[0] for f in files]))
#print(fbg_num)
#fbg_num = ['33971']
for sn in fbg_num:
    df = pd.read_csv(path+"T_"+sn+".csv")
    dfb = pd.read_csv(path+"T_"+sn+"_bare.csv")
    
    lam = np.array(df["Wavelength(nm)"])
    power = np.array(df["Power(dBm)"])
    input_power = np.array(dfb["Power(dBm)"])

    try:
        percent = (toWatt(power)/toWatt(input_power))**2
        #percent = (input_power/power)**2
    except:
        if len(power) > len(input_power):
            power = downsize(power, len(input_power))
            lam = downsize(lam, len(input_power))
        else:
            input_power = downsize(input_power, len(power))
            
    percent = (toWatt(power)/toWatt(input_power))**2
        
    plot_spectrum(lam, power, percent, path, sn)
    try:
        dfr = pd.read_csv(path+"T_"+sn+"_rev.csv")
        power = dfr["Power(dBm)"]
    
        percent = (toWatt(power)/toWatt(input_power))**2
        #percent = (input_power/power)**2
        plot_spectrum(lam, power, percent, path, sn+"_rev")
    except:
        print('no rev')
    
