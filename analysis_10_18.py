from fbg_analysis import *
path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-18/"
fbg_num = ["cav5nm_", "cav5nm_501_", "cav5nm_1001_1k_",  "cav5nm_1001_100k_", "cav5nm_", "cav5nm_5001_"]

for sn in fbg_num:
    i=0
    idx = True
    in_pows = []
    pows = []
    pers = []
    while(idx):
        try:
            dfb = pd.read_csv(path+"T_"+sn+"bare_"+str(i)+".csv")
        except:
            print("T_"+sn+"bare_"+str(i)+" not existing")
            idx = False
        lam = np.array(dfb["Wavelength(nm)"])
        input_power = np.array(dfb["Power(dBm)"])
        i += 1
        in_pows.append(np.array(input_power))
    ip_avg = np.mean(np.array(in_pows), axis=0)

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
        pows.append(np.array(power))
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
    ip_avg_d = ip_avg - diff
    per_avg = (toWatt(pow_avg)/toWatt(ip_avg_d))**2

    plot_spectrums(lam, pows, in_pows, pers, path, sn)
    plot_spectrum(lam, pow_avg, ip_avg_d, per_avg, path, sn)

