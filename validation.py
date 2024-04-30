from Real_RC_Function import *
import datetime

starttime = datetime.datetime.now()

emissivity = pd.read_excel('Emissivity 0.3~20.xlsx')
emissivity_w = emissivity.loc[:, 'Wavelength']
emissivity_e = emissivity.loc[:, 'Emissivity']
emissivity_w = np.array(emissivity_w) * 10 ** -6
emissivity_e = np.array(emissivity_e)


SSSI = pd.read_excel('AM0AM1_5.xls')
SSSI = SSSI.drop(0)  
SSSI = SSSI.drop(['Unnamed: 1', 'Unnamed: 3', 'Unnamed: 4', 'ASTM E-490 AM0 Standard Spectra', 'Unnamed: 6'], axis=1)
SSSI.set_axis(['Wavelength', 'Irradiation'], axis=1, inplace=True)
s_w = SSSI.loc[:, 'Wavelength']
s_r = SSSI.loc[:, 'Irradiation']
SR = simps(s_r, s_w)


R_solar = pd.read_excel('R1_solar.xlsx')
R_solar_t = R_solar.loc[:, 'time']
R_solar_r = R_solar.loc[:, 'solar radiation']


solar_e = []
for i in range(0, len(R_solar_t)):
    solar_e.append(R_solar_r[i] / SR)


R_air = pd.read_excel('R1_air_t.xlsx')
R_air_time = R_air.loc[:, 'time']
R_air_t = R_air.loc[:, 'air temperature']
R_air_time_np = np.array(R_air_time)
R_air_t_np = np.array(R_air_t)


w_min = 0.3e-6
w_max = 20e-6
tcwv = 1000


matched_air_time = []
matched_air_t = []
for t in R_solar_t:
    matched_air_time.append(t)
    if t in R_air_time:
        matched_air_t.append(R_air_t[np.where(t == R_air_time_np)])
    else:
        t_md = min(R_air_time_np, key=lambda x: abs(x - t))
        if t > t_md:
            a = np.where(t_md == R_air_time_np)
            t_md_1 = R_air_t_np[a] + (R_air_t_np[a[0] + 1] - R_air_t_np[a]) / (R_air_time_np[a[0] + 1] - R_air_time_np[a]) * (
                    t - R_air_time_np[a])
            matched_air_t.append(t_md_1[0])
        else:
            a = np.where(t_md == R_air_time_np)
            t_md_1 = R_air_t_np[a[0] - 1] + (R_air_t_np[a[0] - 1] - R_air_t_np[a]) / (
                    R_air_time_np[a[0] - 1] - R_air_time_np[a]) * (t - R_air_time_np[a])
            matched_air_t.append(t_md_1[0])



P_br = []
P_RC = []
P_sun = []
P_atm = []
P_c = []
P_cool = []
t_cooler_min = []

for i in range(0, len(R_solar_t)):
    P_br_m, P_RC_m, P_atm_m, P_sun_m, P_c_m, t_cooler_min_m, P_cool_m = ClaRealPC(matched_air_t[i] + 273.15, 0, s_w, s_r, solar_e[i], tcwv,emissivity_w, emissivity_e, w_min, w_max)
    P_br.append(P_br_m)
    P_RC.append(P_RC_m)
    P_sun.append(P_sun_m)
    P_atm.append(P_atm_m)
    P_c.append(P_c_m)
    P_cool.append(P_cool_m)
    t_cooler_min.append(t_cooler_min_m)
    print('i = ', i)

P_RC = np.array(P_RC)
P_atm = np.array(P_atm)
P_br = np.array(P_br)
P_c = np.array(P_c)
P_cool = np.array(P_cool)
P_sun = np.array(P_sun)
t_cooler_min = np.array(t_cooler_min)


endtime = datetime.datetime.now()
print(endtime - starttime).seconds