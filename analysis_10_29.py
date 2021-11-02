from fbg_analysis import *

path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-29/"
#-----10-06 data-----
dfb = pd.read_csv(path+"T_cav2nm_bare_0")
lam = np.array(dfb["Wavelength(nm)"])
input_power = np.array(dfb["Power(dBm)"])

df = pd.read_csv(path+"T_cav2nm_0")
lam = np.array(df["Wavelength(nm)"])
power = np.array(df["Power(dBm)"])
percent = (toWatt(power)/max(toWatt(power)))**2
#plt.plot(lam, percent, linewidth=1, color='green') #label=fbg_num[0]+str(i))

dfb = pd.read_csv(path+"T_33966a33971_cav_2nm_bare_5.csv")
lam = np.array(dfb["Wavelength(nm)"])
input_power = np.array(dfb["Power(dBm)"])

df = pd.read_csv(path+"T_33966a33971_cav_2nm_4.csv")
lam = np.array(df["Wavelength(nm)"])
power = np.array(df["Power(dBm)"])
percent = (toWatt(power)/max(toWatt(power)))**2
#plt.plot(lam, percent, linewidth=1, color='purple') #label=fbg_num[0]+str(i))

path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/"
fbg_num = ["2021-10-12/T_33966_33971_cav_2nm_", "2021-10-29/T_33966_33971_"]
color = ['green', 'red']
label = ['L_Gap_53.33mm', 'L_Gap_10mm']
for sn in range(len(fbg_num)):
    i=0
    idx = True
    in_pows = []
    pows = []
    pers = []
    while(idx):
        try:
            dfb = pd.read_csv(path+fbg_num[sn]+"bare_"+str(i)+".csv")
        except:
            print("T_"+fbg_num[sn]+"bare_"+str(i)+" not existing")
            idx = False
            break
        lam = np.array(dfb["Wavelength(nm)"])
        input_power = np.array(dfb["Power(dBm)"])
        i += 1
        in_pows.append(np.array(input_power))
    ip_avg = np.mean(np.array(in_pows), axis=0)

    i=0
    idx = True
    while(idx):
        try:
            df = pd.read_csv(path+fbg_num[sn]+str(i)+".csv")
        except:
            print("T_"+fbg_num[sn]+str(i)+" not existing")
            idx = False
            break
        lam = np.array(df["Wavelength(nm)"])
        power = np.array(df["Power(dBm)"])
        pows.append(np.array(power))
        diff = min(ip_avg) - max(power)
        lin_power = toWatt(power)
        #percent = (lin_power/max(lin_power))**2
        percent = (lin_power/max(lin_power[int(len(power)/2):]))**2
        percent[:int(len(power)/2)] = (lin_power[:int(len(power)/2)]/max(lin_power[:int(len(power)/2)]))**2

        pers.append(percent)
        #plt.plot(lam, ip_avg-diff, linewidth=1, color='red') #label=fbg_num[0]+str(i))
        #plt.plot(lam, percent, linewidth=1, color='red') #label=fbg_num[0]+str(i))
        i += 1
    pow_avg = np.mean(np.array(pows), axis=0)

    diff = min(ip_avg) - max(pow_avg)
    ip_avg_d = ip_avg + diff
    per_avg = (toWatt(pow_avg)/toWatt(ip_avg_d))**2

    plt.plot(lam, np.max(pers, axis=0), linewidth=1, color= color[sn], label=label[sn])
    plot_spectrums(lam, pows, in_pows, pers, path+"2021-10-29/", label[sn])
    plot_spectrum(lam, pow_avg, ip_avg_d, per_avg, path+"2021-10-29/", label[sn])


plt.title("Anritsu OSA")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Transmission")
plt.xlim([1038, 1040])
plt.ylim([-0.1, 1.1])
plt.legend()
plt.show()

#---Moglabs----
path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-28/"
dat = np.loadtxt(path+"T_33768_0.dat")
dat = np.transpose(dat)
intensity = dat[1]

fwhm = 2.37 #nm
bragg = 1039.9 #nm
percent = intensity/max(intensity)
l, r, c = find_fwhm(intensity/max(intensity))

diff = r-l
res = 4.8e-3#fwhm/diff
start = bragg - c*res
end = bragg + (len(intensity) - c)*res
lam = np.linspace(start, end, len(intensity))
start = np.where(lam <= 1038)[-1][-1]
end = np.where(lam <= 1040)[-1][-1]

path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-29/"
fbg_num = ["T_33966_33971_"]
i=0
idx = True
in_pows = []
pows = []
pers = []
'''
while(idx):
    try:
        dat = np.loadtxt(path+fbg_num[0]+"bare_"+str(i)+".dat")
        dat = np.transpose(dat)
        input_intensity = dat[1]
    except:
        print(fbg_num[0]+"bare_"+str(i)+" not existing")
        idx = False
    i += 1
    in_pows.append(np.array(input_intensity))
ip_avg = np.mean(np.array(in_pows), axis=0)
ip_avg = np.ones(len(intensity))*np.mean(ip_avg)
'''

i=0
idx = True
while(idx):
    try:
        dat = np.loadtxt(path+fbg_num[0]+str(i)+".dat")
        dat = np.transpose(dat)
        intensity = dat[1]
    except:
        print(fbg_num[0]+str(i)+" not existing")
        idx = False
        break
    ip_avg = np.ones(len(intensity))*max(intensity[start:end])#1650
    pows.append(np.array(intensity))
    diff = min(ip_avg) - max(intensity)
    percent = ((np.array(intensity))/(ip_avg))**2
    pers.append(percent)

    #plt.plot(lam, percent, linewidth=1, color='blue') #label=fbg_num[0]+str(i))
    i += 1
pow_avg = np.mean(np.array(pows), axis=0)

plt.plot(lam, np.max(pers, axis=0), linewidth=1, color='blue') #label=fbg_num[0]+str(i))
plot_spectrum(lam, intensity, ip_avg, percent, path, "moglabs")
plot_spectrums(lam, pows, in_pows, pers, path, "moglabs")

plt.title("Moglabs WMW")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Transmission")
plt.xlim([1038, 1040])
plt.ylim([-0.1, 1.1])
#plt.legend()
plt.show()
