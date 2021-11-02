from fbg_analysis import *
path = "/run/media/anya/Seagate Backup Plus Drive/NPQO/Anritzu OSA/Data/2021-10-07/"
#path = "C:/Users/Anya/Desktop/Anritzu OSA/Data/2021-10-07/"
#files = glob.glob(path + '*.csv')
#fbg_num = list(set([f.split('_')[1].split('.')[0] for f in files]))
#print(fbg_num)
fbg_num = ["33212a33767_cav_", "33966a33971_cav_2nm_", "33212_", "33767_"]

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
        #plt.plot(lam, input_power, linewidth=1, label="bare"+str(i))
        in_pows.append(np.array(input_power))
    ip_avg = np.mean(np.array(in_pows), axis=0)
    #plt.show()

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
        #plt.plot(lam, power, linewidth=1, label="power"+str(i))
        pows.append(np.array(power))
        try:
            percent = (toWatt(power)/toWatt(ip_avg))**2
            #percent = (input_power/power)**2
        except:
            if len(power) > len(ip_avg):
                power = downsize(power, len(ip_avg))
                lam = downsize(lam, len(ip_avg))
            else:
                ip_avg_d = downsize(ip_avg, len(power))
            percent = (toWatt(power)/toWatt(ip_avg_d))**2
        pers.append(percent)
    pow_avg = np.mean(np.array(pows), axis=0)
    #plt.show()

    per_avg = (toWatt(pow_avg)/toWatt(ip_avg))**2

    plot_spectrums(lam, pows, in_pows, pers, path, sn)
    plot_spectrum(lam, pow_avg, ip_avg, per_avg, path, sn)


