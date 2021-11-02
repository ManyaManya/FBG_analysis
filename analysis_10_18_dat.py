from fbg_analysis import *
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

i=0
idx = True
in_pows = []
pows = []
pers = []
while(idx):
    try:
        dat = np.loadtxt(path+"bare"+str(i)+".dat")
        dat = np.transpose(dat)
        input_intensity = dat[1]
    except:
        print("bare"+str(i)+" not existing")
        idx = False
    i += 1
    in_pows.append(np.array(input_intensity))
ip_avg = np.mean(np.array(in_pows), axis=0)
ip_avg = np.ones(len(intensity))*np.mean(ip_avg)

i=0
idx = True
while(idx):
    try:
        dat = np.loadtxt(path+"moglabs_spectrum_"+str(i)+".dat")
        dat = np.transpose(dat)
        intensity = dat[1]
    except:
        print("moglabs_spectrum_"+str(i)+" not existing")
        idx = False
    i += 1
    pows.append(np.array(intensity))
    diff = min(ip_avg) - max(intensity)
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

plot_spectrum(lam, intensity, ip_avg, percent, path, "moglabs")
plot_spectrums(lam, pows, in_pows, pers, path, "moglabs")

'''
plt.plot(lam, percent)
plt.plot(lam[[l, r, c]], percent[[l, r, c]], "ro")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Percent (nm)")
plt.show()
'''
