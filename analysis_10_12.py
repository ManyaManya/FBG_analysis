from fbg_analysis import *
path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-12/"
fbg_num = ["33966_33971_cav_1nm_", "33966_33971_cav_1nm_none_", "33966_33971_cav_1nm_rev_", \
           "33966_33971_cav_1nm_rev_dir_", "33966_33971_cav_1nm_rev_dir0_", "33966_33971_cav_2nm_"]

for sn in fbg_num:
    in_pows = []
    pows = []
    pers = []
    i=0
    idx = True
    while(idx):
        try:
            df = pd.read_csv(path+"T_"+sn+str(i)+".csv")
        except:
            print("T_"+sn+str(i)+" not existing")
            idx = False
        lam = np.array(df["Wavelength(nm)"])
        power = np.array(df["Power(dBm)"])
        i += 1
        power_watt = toWatt(power)
        percent = (power_watt/max(power_watt))**2
        pows.append(np.array(power))
        pers.append(percent)

    pow_avg = np.mean(np.array(pows), axis=0)
    per_avg = (toWatt(pow_avg)/max(toWatt(pow_avg)))**2

    plot_spectrums(lam, pows, in_pows, pers, path, sn)
    plot_spectrum(lam, pow_avg, np.ones(len(pow_avg))*max(pow_avg), per_avg, path, sn)

