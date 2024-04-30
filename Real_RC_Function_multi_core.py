import numpy as np
import pandas as pd
from scipy.integrate import simps
from scipy.integrate import quad

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
s = s_r
SR = simps(s_r, s_w)


w_min = 0.3e-6
w_max = 20e-6


class Radiation:


    def __init__(self, T):

        self.T = T

    def I(self, w):

        h = 6.62607004e-34  
        k_B = 1.38064852e-23  
        c = 3.0e8  
        return 2 * h * c ** 2 / w ** 5 * 1 / (np.exp(h * c / (w * k_B * self.T)) - 1)



def x(i, n, a, b):  
    return a + (b - a) * i / n


def y(i, n, a, b):  
    return a + (b - a) * i / n


def find_e(w, rc):

    if w in rc.e_w:
        e = rc.e[np.where(w == rc.e_w)]
        e_md = e[0]
    else:
        w_md = min(rc.e_w, key=lambda x: abs(x - w))
        if w > w_md:
            a = np.where(w_md == rc.e_w)
            e = rc.e[a] + (rc.e[a[0] + 1] - rc.e[a]) / (rc.e_w[a[0] + 1] - rc.e_w[a]) * (w - rc.e_w[a])
        else:
            a = np.where(w_md == rc.e_w)
            e = rc.e[a] + (rc.e[a[0] - 1] - rc.e[a]) / (rc.e_w[a[0] - 1] - rc.e_w[a]) * (w - rc.e_w[a])
        e_md = e[0]
    return e_md


def calculatePRC(rc, w_min, w_max):


    y = []
    for w in np.arange(w_min, w_max, 1e-6):
        e_md = find_e(w, rc)
        y_md = e_md * rc.I(w)
        y.append(y_md)
    
    w = np.arange(w_min, w_max, 1e-6)
    P_RC = simps(y, w)
    P_RC = P_RC * np.pi
    return P_RC


def int_cylinder(m, n, ax, bx, ay, by, T_air, rc, TR_np):

    ax = float(ax)
    bx = float(bx)
    ay = float(ay)
    by = float(by)
    Air = Radiation(T_air)
    
    
    f_sum = 0.0
    count = 0
    TR_w_list = TR_np[:, 0]
    TR_t_list = TR_np[:, 1]
    rc.e_w = rc.e_w * 10 ** 6  
    for i in range(m):
        xx = (x(i, m, ax, bx) + x(i + 1, m, ax, bx)) / 2
        for j in range(n):
            yy = (y(j, n, ay, by) + y(j + 1, n, ay, by)) / 2
            if yy in TR_w_list:
                TR_t = TR_t_list[np.where(TR_w_list == yy)]
            else:
                yy_min = min(TR_w_list, key=lambda w: abs(w - yy))
                if yy > yy_min:
                    a = np.where(TR_w_list == yy_min)
                    TR_t = TR_t_list[a] + (TR_t_list[a[0] + 1] - TR_t_list[a]) / (
                            TR_w_list[a[0] + 1] - TR_w_list[a]) * (yy - TR_w_list[a])
                else:
                    a = np.where(TR_w_list == yy_min)
                    TR_t = TR_t_list[a] + (TR_t_list[a[0] - 1] - TR_t_list[a]) / (
                            TR_w_list[a[0] - 1] - TR_w_list[a]) * (yy - TR_w_list[a])

            e_md = find_e(yy, rc)
            f = Air.I(yy * 10 ** -6) * e_md * (1 - TR_t ** (1 / np.cos(xx))) * np.sin(xx) * np.cos(xx)
            count = count + 1
            
            
            
            
            f_sum += f
            
            

    result = 2 * np.pi * f_sum * ((bx - ax) * 10 ** -6 / m) * ((by - ay) / n)
    
    rc.e_w = rc.e_w / 10 ** 6  
    return result


class Region:


    def __init__(self, sr, ws, t, twvc):

        self.sr = sr
        self.ws = ws
        self.t = t
        self.twvc = twvc


def dichotomy(left, right, eps, T_air, v_air, P_sun, P_atm, w_min, w_max, rc):

    middle = (left + right) / 2
    count = 0  

    def solve_function(x):

        h_c = 5.7 + 3.8 * v_air
        P_c = h_c * (T_air - x)
        rc_md = Radiation(x)
        rc_md.e = rc.e
        rc_md.e_w = rc.e_w
        D = calculatePRC(rc_md, w_min, w_max) - P_sun - P_atm - P_c
        return D

    while abs(solve_function(middle)) > eps:
        middle = (left + right) / 2
        if solve_function(left) * solve_function(middle) <= 0:
            right = middle
        else:
            left = middle
        count = count + 1
    return count, middle


def ClaRealPC(T_air, v_air, scoef, tcwv):
    T_cooler = T_air  

    
    s_r = s * scoef

    
    TR = pd.read_excel('%d' % (tcwv) + 'atm-cm.xlsx')  
    TR_np = np.array(TR)
    TR.set_axis(['Wavelength', 'Transmission'], axis=1, inplace=True)
    TR_w = TR.loc[:, 'Wavelength']
    TR_t = TR.loc[:, 'Transmission']

    
    rc = Radiation(T_cooler)
    P_br, err = quad(rc.I, w_min, w_max)
    P_br = P_br * np.pi  

    
    rc.e_w = emissivity_w
    rc.e = emissivity_e
    P_RC = calculatePRC(rc, w_min, w_max)

    
    
    P_atm_array = int_cylinder(50, 500, 0, np.pi / 2, min(TR_w), max(TR_w), T_air, rc, TR_np)
    P_atm = P_atm_array[0]

    
    P_sun_mid = []
    for index in range(len(s_w)):
        index = index + 1
        e = find_e(s_w[index] * 10 ** -9, rc)
        P_sun_mid.append(e * s_r[index])

    P_sun = simps(P_sun_mid, s_w)

    
    
    P_cool = P_RC - P_sun - P_atm

    
    
    P_RC_min = P_sun + P_atm  

    if P_RC_min < P_RC:
        count, T_cooler_min = dichotomy(T_air - 200, T_air, 10e-4, T_air, v_air, P_sun, P_atm, w_min, w_max, rc)
    elif P_RC_min > P_RC:
        count, T_cooler_min = dichotomy(T_air, T_air + 200, 10e-4, T_air, v_air, P_sun, P_atm, w_min, w_max, rc)
    else:
        T_cooler_min = T_air

    h_c = 5.7 + 3.8 * v_air  
    P_c = h_c * (T_air - T_cooler_min)
    print('P_RC = ', P_RC)
    print('P_atm = ', P_atm)
    print('P_sun = ', P_sun)
    print('P_c = ', P_c)
    print('P_cool = ', P_cool)
    print('T_cooler_min = ', T_cooler_min - 273.15)
    print('T_cooler_min - T_air = ', T_cooler_min - T_air)
    return P_br, P_RC, P_atm, P_sun, P_c, T_cooler_min, P_cool
