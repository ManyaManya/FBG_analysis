import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, pi, hbar, m_e
import pandas as pd
from scipy.interpolate import interp1d

from FBG_analysis import fbg_analysis as a

if __name__ == "__main__":
#---Moglabs----
    path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-28/"
    dat = np.loadtxt(path+"T_33768_0.dat")
    dat = np.transpose(dat)
    intensity = dat[1]

    fwhm = 2.37 #nm
    bragg = 1039.9 #nm
    percent = intensity/max(intensity)
    l, r, c = a.find_fwhm(intensity/max(intensity))

    diff = r-l
    res = 4.8e-3#fwhm/diff
    start = bragg - c*res
    end = bragg + (len(intensity) - c)*res
    lam = np.linspace(start, end, len(intensity))
    start = np.where(lam <= 1038)[-1][-1]
    end = np.where(lam <= 1040)[-1][-1]

    path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-27/"
    fbg_num = ["T_33966_33971_Yt_"]#,"R_33966_33971_Yt_0_", "R_33966_33971_Yt_1_", "R_33966_33971_Yt_2_"]#["thresh_"]
    pers = []
    lam = lam[start:end]
    #label=["Above Threshold", "Below Threshold"]
    for sn in range(len(fbg_num)):
        i=0
        idx = True
        while(idx):
            try:
                dat = np.loadtxt(path+fbg_num[sn]+str(i)+".dat", skiprows=2)
                dat = np.transpose(dat)
                #lam = dat[0]
                intensity = dat[1]
                intensity = intensity[start:end]
            except:
                print(fbg_num[0]+str(i)+" not existing")
                idx = False
                break
            ip_avg = np.ones(len(intensity))*1600#max(intensity)
            percent = (intensity/ip_avg)**2
            pers.append(percent)

            #plt.plot(lam, intensity/max(intensity), linewidth=1, label=label[i])
            i += 1

    per = np.min(pers, axis=0)
    l, r, c = a.find_fwhm(np.log10(per))
    plt.plot(lam, per, linewidth=2, label="R="+str(round(1-min(per), 2)))
    plt.vlines(lam[l], ymin=-0.1, ymax=1.1, color='green', linewidth=1)
    plt.vlines(lam[r], ymin=-0.1, ymax=1.1, color='green', linewidth=1)
    plt.hlines(per[r], xmin=lam[l], xmax=lam[r], color='green', linewidth=1, label='FWHM='+str(round(lam[r]-lam[l], 2))+'nm')
    #a.plot_spectrum(lam, intensity, ip_avg, np.max(pers, axis=0), path, "Reflectance")

    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Intensity (Normalized)")
    #plt.xlim([900, 1100])
    plt.xlim([1038, 1040])
    plt.ylim([-0.1, 1.1])
    plt.legend()
    plt.show()



    path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-28/stellarnet/"
    fbg_num = [["910nm_above_768.SSM", "910nm_above_T_768.SSM"], \
               ["910nm_above_T.SSM", "910nm_below_T.SSM"],\
               ["910nm_above_768.SSM", "940nm_above_768.SSM","950nm_above_768.SSM", "975nm_above_768.SSM"]]
    pers = []
    #lam = lam[start:end]
    label=[["Above Threshold (Reflectance)", "Above Threshold (Transmission)"],\
           ["Above Threshold", "Below Threshold"],\
           ["910nm", "940nm", "950nm", "975nm"]]
    for sn in range(len(fbg_num)):
        for i in range(len(fbg_num[sn])):
            try:
                dat = np.loadtxt(path+fbg_num[sn][i], skiprows=2)
                dat = np.transpose(dat)
                lam = dat[0]
                intensity = dat[1]
            except:
                print(fbg_num[sn][i]+" not existing")
                break
            ip_avg = np.ones(len(intensity))*max(intensity)
            percent = ((np.array(intensity))/(ip_avg))**2
            pers.append(percent)

            plt.plot(lam, intensity/max(intensity), linewidth=1, label=label[sn][i])
            i += 1

        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Intensity (Normalized)")
        plt.xlim([850, 1100])
        plt.ylim([-0.1, 1.1])
        plt.legend()
        plt.show()
