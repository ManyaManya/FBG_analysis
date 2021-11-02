
from fbg_analysis import *

path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-28/"

fbg_num = ["33768_", "33768_33163_"]
for sn in fbg_num:
    i=0
    idx = True
    in_pows = []
    pows = []
    pers = []
    '''
    while(idx):
        try:
            dfb = pd.read_csv(path+"T_"+sn+"bare_"+str(i))
        except:
            print("T_"+sn+"bare_"+str(i)+" not existing")
            idx = False
        lam = np.array(dfb["Wavelength(nm)"])
        input_power = np.array(dfb["Power(dBm)"])
        i += 1
        in_pows.append(np.array(input_power))
    ip_avg = np.mean(np.array(in_pows), axis=0)
    '''
    i=0
    idx = True
    while(idx):
        try:
            df = pd.read_csv(path+"T_"+sn+str(i))
        except:
            print("T_"+sn+str(i)+" not existing")
            idx = False
        lam = np.array(df["Wavelength(nm)"])
        power = np.array(df["Power(dBm)"])
        i += 1
        pows.append(np.array(power))
        ip_avg = np.ones(len(power))*max(power)
        diff = np.abs(min(ip_avg) - max(power))
        try:
            percent = (toWatt(power)/toWatt(ip_avg-diff))**2
        except:
            if len(power) > len(ip_avg):
                power = downsize(power, len(ip_avg))
                lam = downsize(lam, len(ip_avg))
            else:
                ip_avg = downsize(ip_avg, len(power))
            percent = (toWatt(power)/toWatt(ip_avg-diff))**2
        pers.append(percent)
    pow_avg = np.mean(np.array(pows), axis=0)

    diff = min(ip_avg) - max(pow_avg)
    print(diff)
    ip_avg_d = ip_avg - diff - 0.5
    per_avg = (toWatt(pow_avg)/toWatt(ip_avg_d))**2

    plot_spectrums(lam, pows, in_pows, pers, path, sn)
    plot_spectrum(lam, pow_avg, ip_avg_d, per_avg, path, sn)


plt.plot(lam, per_avg, linewidth=1, color='red', label='Anritsu')

#---Moglabs----
path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-18/"
dat = np.loadtxt(path+"moglabs_spectrum_0.dat")
dat = np.transpose(dat)
intensity = dat[1]

fwhm = 2.96 #nm
bragg = 1039.04 #nm
percent = intensity/max(intensity)
l, r, c = find_fwhm(intensity/max(intensity))

diff = r-l
res = fwhm/diff
start = bragg - c*res
end = bragg + (len(intensity) - c)*res
lam = np.linspace(start, end, len(intensity))

path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-28/"
fbg_num = ["T_33768_", "T_33768_33163_", "T_33966_33971_"]

for sn in fbg_num:
    i=0
    idx = True
    in_pows = []
    pows = []
    pers = []
    while(idx):
        try:
            dat = np.loadtxt(path+sn+"bare_"+str(i)+".dat")
            dat = np.transpose(dat)
            input_intensity = dat[1]
        except:
            print(sn+"bare_"+str(i)+" not existing")
            idx = False
        i += 1
        in_pows.append(np.array(input_intensity))
    ip_avg = np.mean(np.array(in_pows), axis=0)
    ip_avg = np.ones(len(intensity))*np.mean(ip_avg)

    i=0
    idx = True
    while(idx):
        try:
            dat = np.loadtxt(path+sn+str(i)+".dat")
            dat = np.transpose(dat)
            intensity = dat[1]
        except:
            print(sn+str(i)+" not existing")
            idx = False
        i += 1
        pows.append(np.array(intensity))
        diff = 0#min(ip_avg) - max(intensity)
        try:
            percent = (intensity/(ip_avg-diff))**2
        except:
            if len(intensity) > len(ip_avg):
                intensity = downsize(intensity, len(ip_avg))
                lam = downsize(lam, len(ip_avg))
            else:
                ip_avg = downsize(ip_avg, len(intensity))
            percent = (intensity/(ip_avg-diff))**2
        pers.append(percent)
    pow_avg = np.mean(np.array(pows), axis=0)

    plot_spectrum(lam, intensity, ip_avg, percent, path, "moglabs"+sn)
    plot_spectrums(lam, pows, in_pows, pers, path, "moglabs"+sn)
    plt.plot(lam, percent*100, linewidth=1,  color='blue', label='Moglabs')

#plt.xlim([1038, 1040])
plt.legend()
plt.show()




