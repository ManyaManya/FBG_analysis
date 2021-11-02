from fbg_analysis import *

path = "C:/Users/Anya/Desktop/Anritzu OSA/Data/2021-10-06/"
files = glob.glob(path + '*.csv')
fbg_num = list(set([f.split('_')[1].split('.')[0] for f in files]))

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

