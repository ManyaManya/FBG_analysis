from fbg_analysis import *
path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-13_set2/"

dat = np.loadtxt(path+"5.dat")
dat = np.transpose(dat)
intensity = dat[1]

fwhm = 2.88 #nm
bragg = 1038.92 #nm
percent = intensity/max(intensity)
l, r, c = find_fwhm(intensity/max(intensity))

diff = r-l
res = fwhm/diff
start = bragg - c*res
end = bragg + (len(intensity) - c)*res
lam = np.linspace(start, end, len(intensity))

ip_avg = np.ones(len(intensity))*max(intensity)

plot_spectrum(lam, intensity, ip_avg, percent, path, "moglabs")
'''
plt.plot(lam, percent)
plt.plot(lam[[l, r, c]], percent[[l, r, c]], "ro")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Percent (nm)")
plt.show()
'''
