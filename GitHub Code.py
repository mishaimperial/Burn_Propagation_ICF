
""" 
Python Code written to analyse Chimera output data
Functions have been included, but for neatness, commands to execute the functions
have not

"""

"""
CONTENTS (code line in brackets):
    
Load modules (9)
Load data (18)
Colour plots (293)
Burn diagnostics (308)
Hotspot diagnostics (352)
Implosion plots (432)
Curve fit of burn pulse (543)
Hotspot-to-total DT ratio (595)
Power density plots for specified times (629)
Comparison of BUF and non-BUF S=1 yield (685)
Comparison of BUF yields for all scale factors (696)
Comparison of how much the yield is reduced for BUF vs non-BUF for all scale factors (715)
Changing coast time by changing frequency dependent spectrum (733)
Find burn front velocity after peak compression by using density spikes as the burn front (790)
Trajectory of a capsule in density-temperature space aroumd bang time (808)
Normalised integrated density vs normalised fusion plots (869)
Animation of integral plots (932)
Normalised integrated density and fusion for scale factor comparisons at bang time (952)
Fusion colour plots using the imshow and colormesh methods (994)
FWHM upper and lower bound 98% yield plots (1041)
Preheat off vs on for S=1 (1065)
Peak compression times for preheat on and off scale factors, and for different ice thicknesses (1075)
Scale factor comparison of diagnostics (1102)
Preheat off scale factor comparison of diagnostics (1138)
Scale Factor Comparison of Diagnostics for varying Ice Thickness (1173)
Mass ratio comparison for all scale factors (1240)
Burn-averaged temperature comparison for all scale factors (1259)
Burn pulse comparison for all scale factors (1273)
Preheat off mass ratio comparison for all scale factors (1290)
Preheat off burn-averaged temperature comparisons for all scale factors (1305)
Preheat off burn pulse comparisons for all scale factors (1319)
Ice thickness variation mass ratio comparisons (1336)
Ice thickness variation burn-averaged temperature comparisons (1351)
Ice thickness variation burn pulse comparisons (1365)
BUF mass ratio comparisons for all scale factors (1382)
Mass ablation rates for all scale factors (1397)
PVGamma calculation (1454)
Sound speed and burn front speed diagnostics (1477)
Initial toy problem density and temperature configuration (1801)
Density, temperature, pressure and power density animation (1817)
Temperature, density and pressure evolution plots (1860)
Temperature, density and pressure evolution plots (1918)
"""

#%% Load Modules 

import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import colors
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import time 

#%% Load Data 

""" Data for 'diagnostics' function- changed manually depending on which scale factor was being studied"""

work_elec = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04work_elec.dat')
work_ion = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04work_ion.dat')
data = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
tion = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04T_ion.dat')
telec = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04T_elec.dat')
rho = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04rho_DT.dat')
rhocc = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04rho_CC.dat')
qtot = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04alpha_qTot.dat')
elecwork = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04work_elec_cond.dat')
deltot = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04del_E_tot.dat')
pdv = work_elec + work_ion
yieldrate = data[2]
Yield = data[3]
dt = data[1]
burnavTi = data[4]
time = data[0]
scalefactor = 0.014696665446667919 #This is what to multiply radiation power by in implosion plots- it's the scale factor by which the source code is off 

""" Data for all scale factors for BUF and non-BUF """

index = 2
data_1 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_1 = data_1[index]
time_1 = data_1[0]
data_1_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy67\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_1_buf = data_1_buf[index]
time_1_buf = data_1_buf[0]
rho_1 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04rho_DT.dat')
rho_1_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy67\xy00\xy00rho_DT.dat')
data_09 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy11\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_09 = data_09[index]
time_09 = data_09[0]
data_09_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy72\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_09_buf = data_09_buf[index]
time_09_buf = data_09_buf[0]
rho_09 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy07\xy00\xy00rho_DT.dat')
rho_09_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy72\xy00\xy00rho_DT.dat')
data_11 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy07\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_11 = data_11[index]
time_11 = data_11[0]
data_11_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy68\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_11_buf = data_11_buf[index]
time_11_buf = data_11_buf[0]
rho_11 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy07\xy00\xy00rho_DT.dat')
rho_11_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy68\xy00\xy00rho_DT.dat')
data_12 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy08\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_12 = data_12[index]
time_12 = data_12[0]
data_12_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy69\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_12_buf = data_12_buf[index]
time_12_buf = data_12_buf[0]
rho_12 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy08\xy00\xy04rho_DT.dat')
rho_12_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy69\xy00\xy00rho_DT.dat')
data_13 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy09\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_13 = data_13[index]
time_13 = data_13[0]
data_13_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy70\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_13_buf = data_13_buf[index]
time_13_buf = data_13_buf[0]
rho_13 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy09\xy00\xy00rho_DT.dat')
rho_13_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy70\xy00\xy00rho_DT.dat')
data_14 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy10\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_14 = data_14[index]
time_14 = data_14[0]
data_14_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy71\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
yield_14_buf = data_14_buf[index]
time_14_buf = data_14_buf[0]
rho_14 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy10\xy00\xy00rho_DT.dat')
rho_14_buf = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy71\xy00\xy00rho_DT.dat')
tion_09 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy11\xy00\xy00T_ion.dat')
tion_1 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04T_ion.dat')
tion_11 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy07\xy00\xy00T_ion.dat')
tion_12 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy08\xy00\xy04T_ion.dat')
tion_13 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy09\xy00\xy00T_ion.dat')
tion_14 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy10\xy00\xy00T_ion.dat')
fusion_09 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy11\xy00\xy00Fusion.dat')
fusion_1 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04Fusion.dat')
fusion_11 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy07\xy00\xy00Fusion.dat')
fusion_12 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy08\xy00\xy04Fusion.dat')
fusion_13 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy09\xy00\xy00Fusion.dat')
fusion_14 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy10\xy00\xy00Fusion.dat')

""" Data for amending coast time from laser radiation driver """

hydra_energy = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\Source Code\Hydra_energy_N210808_difhyd.txt')
hydra_time = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\Source Code\Hydra_time_N210808_difhyd.txt')

""" Scale factor comparisons for preheat off, BUF and ice thickness variation """

data_09 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy11\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_09 = data_09[0]
yieldrate_09 = data_09[2]
burnav_09 = data_09[4]
data_1 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_1 = data_1[0]
yieldrate_1 = data_1[2]
burnav_1 = data_1[4]
rho_1 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04rho_DT.dat')
tion_1 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04T_ion.dat')
data_11 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy07\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_11 = data_11[0]
burnav_11 = data_11[4]
yieldrate_11 = data_11[2]
data_12 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy08\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_12 = data_12[0]
yieldrate_12 = data_12[2]
burnav_12 = data_12[4]
data_13 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy09\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_13 = data_13[0]
yieldrate_13 = data_13[2]
burnav_13 = data_13[4]
data_14 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy10\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_14 = data_14[0]
yieldrate_14 = data_14[2]
burnav_14 = data_14[4]
rho_09 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy11\xy00\xy00rho_DT.dat')
rho_1 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04rho_DT.dat')
rho_11 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy07\xy00\xy00rho_DT.dat')
rho_12 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy08\xy00\xy04rho_DT.dat')
rho_13 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy09\xy00\xy00rho_DT.dat')
rho_14 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy10\xy00\xy00rho_DT.dat')
tion_09 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy11\xy00\xy00T_ion.dat')
tion_1 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy06\xy00\xy04T_ion.dat')
tion_11 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy07\xy00\xy00T_ion.dat')
tion_12 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy08\xy00\xy04T_ion.dat')
tion_13 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy09\xy00\xy00T_ion.dat')
tion_14 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy10\xy00\xy00T_ion.dat')

#data_09 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy11\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
#time_09 = data_09[0]
#yieldrate_09 = data_09[2]
#burnav_09 = data_09[4]
data_1pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy12\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_1pre = data_1pre[0]
yieldrate_1pre = data_1pre[2]
burnav_1pre = data_1pre[4]
data_11pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy13\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_11pre = data_11pre[0]
burnav_11pre = data_11pre[4]
yieldrate_11pre = data_11pre[2]
data_12pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy14\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_12pre = data_12pre[0]
yieldrate_12pre = data_12pre[2]
burnav_12pre = data_12pre[4]
data_13pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy15\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_13pre = data_13pre[0]
yieldrate_13pre = data_13pre[2]
burnav_13pre = data_13pre[4]
data_14pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy16\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_14pre = data_14pre[0]
yieldrate_14pre = data_14pre[2]
burnav_14pre = data_14pre[4]
#rho_09 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy11\xy00\xy00rho_DT.dat')
rho_1pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy12\xy00\xy04rho_DT.dat')
rho_11pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy13\xy00\xy00rho_DT.dat')
rho_12pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy14\xy00\xy04rho_DT.dat')
rho_13pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy15\xy00\xy00rho_DT.dat')
rho_14pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy16\xy00\xy00rho_DT.dat')
#tion_09 = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy11\xy00\xy00T_ion.dat')
tion_1pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy12\xy00\xy04T_ion.dat')
tion_11pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy13\xy00\xy00T_ion.dat')
tion_12pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy14\xy00\xy04T_ion.dat')
tion_13pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy15\xy00\xy00T_ion.dat')
tion_14pre = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy16\xy00\xy00T_ion.dat')

data_1p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy52\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_1p = data_1p[0]
yieldrate_1p = data_1p[2]
burnav_1p = data_1p[4]
data_2p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy53\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_2p = data_2p[0]
yieldrate_2p = data_2p[2]
burnav_2p = data_2p[4]
data_3p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy54\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_3p = data_3p[0]
burnav_3p = data_3p[4]
yieldrate_3p = data_3p[2]
data_4p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy55\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_4p = data_4p[0]
yieldrate_4p = data_4p[2]
burnav_4p = data_4p[4]
data_5p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy56\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_5p = data_5p[0]
yieldrate_5p = data_5p[2]
burnav_5p = data_5p[4]
data_6p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy57\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_6p = data_6p[0]
yieldrate_6p = data_6p[2]
burnav_6p = data_6p[4]
data_7p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy58\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_7p = data_7p[0]
burnav_7p = data_7p[4]
yieldrate_7p = data_7p[2]
data_8p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy59\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_8p = data_4p[0]
yieldrate_8p = data_8p[2]
burnav_8p = data_8p[4]
data_9p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy60\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_9p = data_9p[0]
yieldrate_9p = data_9p[2]
burnav_9p = data_9p[4]
data_10p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy61\xy00\xy04.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_10p = data_10p[0]
yieldrate_10p = data_10p[2]
burnav_10p = data_10p[4]
rho_1p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy52\xy00\xy04rho_DT.dat')
rho_2p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy53\xy00\xy04rho_DT.dat')
rho_3p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy54\xy00\xy04rho_DT.dat')
rho_4p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy55\xy00\xy04rho_DT.dat')
rho_5p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy56\xy00\xy04rho_DT.dat')
rho_6p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy57\xy00\xy04rho_DT.dat')
rho_7p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy58\xy00\xy04rho_DT.dat')
rho_8p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy59\xy00\xy04rho_DT.dat')
rho_9p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy60\xy00\xy04rho_DT.dat')
rho_10p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy61\xy00\xy04rho_DT.dat')
tion_1p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy52\xy00\xy04T_ion.dat')
tion_2p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy53\xy00\xy04T_ion.dat')
tion_3p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy54\xy00\xy04T_ion.dat')
tion_4p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy55\xy00\xy04T_ion.dat')
tion_5p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy56\xy00\xy04T_ion.dat')
tion_6p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy57\xy00\xy04T_ion.dat')
tion_7p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy58\xy00\xy04T_ion.dat')
tion_8p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy59\xy00\xy04T_ion.dat')
tion_9p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy60\xy00\xy04T_ion.dat')
tion_10p = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy61\xy00\xy04T_ion.dat')

data_09BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy72\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_09BUF = data_09BUF[0]
yieldrate_09BUF = data_09BUF[2]
burnav_09BUF = data_09BUF[4]
data_1BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy67\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_1BUF = data_1BUF[0]
yieldrate_1BUF = data_1BUF[2]
burnav_1BUF = data_1BUF[4]
data_11BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy68\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_11BUF = data_11BUF[0]
burnav_11BUF = data_11BUF[4]
yieldrate_11BUF = data_11BUF[2]
data_12BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy69\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_12BUF = data_12BUF[0]
yieldrate_12BUF = data_12BUF[2]
burnav_12BUF = data_12BUF[4]
data_13BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy70\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_13BUF = data_13BUF[0]
yieldrate_13BUF = data_13BUF[2]
burnav_13BUF = data_13BUF[4]
data_14BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy71\xy00\xy00.dat', skiprows=1, unpack=True, usecols=[0, 15, 40, 44, 47])
time_14BUF = data_14BUF[0]
yieldrate_14BUF = data_14BUF[2]
burnav_14BUF = data_14BUF[4]
rho_09BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy72\xy00\xy00rho_DT.dat')
rho_1BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy67\xy00\xy00rho_DT.dat')
rho_11BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy68\xy00\xy00rho_DT.dat')
rho_12BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy69\xy00\xy00rho_DT.dat')
rho_13BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy70\xy00\xy00rho_DT.dat')
rho_14BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy71\xy00\xy00rho_DT.dat')
tion_09BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy72\xy00\xy00T_ion.dat')
tion_1BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy67\xy00\xy00T_ion.dat')
tion_11BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy68\xy00\xy00T_ion.dat')
tion_12BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy69\xy00\xy00T_ion.dat')
tion_13BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy70\xy00\xy00T_ion.dat')
tion_14BUF = np.loadtxt(r'C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\xy files\xy71\xy00\xy00T_ion.dat')

#%% Compute many diagnostics for a capsule 

def diagnostics(SCALEFACTOR):
    dist = np.arange(0, np.shape(tion)[1], 1)
    tarray = np.arange(0, np.shape(tion)[0], 1)
    density = rho/1e5
    ccdens = rhocc/1e5

    #Colour plots
    plot1 = plt.figure(1)
    plt.imshow(tion,norm=colors.LogNorm(1, 1e5), aspect='auto', origin='lower') 
    plt.colorbar(label= 'Temperature/eV')                                                  
    plt.xlabel('Radial Distance (μm)')
    plt.ylabel('Time (in intervals of 20ps)')
    plt.title('Variation of Ion Temperature and Capsule Radius with Time')
    
    plot2 = plt.figure(2)
    plt.imshow(rho,norm=colors.LogNorm(1, 1e5), aspect='auto', origin='lower') 
    plt.colorbar(label= 'Density (kg/m^3')                                                  
    plt.xlabel('Radial Distance (μm)')
    plt.ylabel('Time (in intervals of 20ps)')
    plt.title('Variation of DT Density and Capsule Radius with Time')
    
    #Burn Diagnostics
    def annot_max(x,y, ax=None):
            xmax = x[np.argmax(y)]
            ymax = y.max()
            text= "x={:.4}, y={:.4}".format(xmax, ymax)
            if not ax:
                ax=plt.gca()
                bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
                arrowprops=dict(arrowstyle="-",connectionstyle="angle,angleA=0,angleB=60")
                kw = dict(xycoords='data',textcoords="axes fraction",
                          arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")    
                ax.annotate(text, xy=(xmax, ymax), xytext=(1.01,0.98), **kw)

    def annot_maxinteg(x,y, ax=None):
        xmax = x[np.argmax(y)]
        ymax = y.max()
        text= "x={:.4}, y={:.4}".format(xmax, ymax)
        if not ax:
            ax=plt.gca()
            bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
            arrowprops=dict(arrowstyle="-",connectionstyle="angle,angleA=0,angleB=60")
            kw = dict(xycoords='data',textcoords="axes fraction",
                      arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")    
            ax.annotate(text, xy=(xmax, ymax), xytext=(0.47,0.98), **kw)
    
    plot3 = plt.figure(3)  
    
    annot_max(time, yieldrate)
    plt.plot(time, yieldrate, color='navy')
    #plt.xlim(0.8e-8, 0.83e-8)
    plt.xlabel('Time (s)')
    plt.ylabel('Yield Rate (Neutron/s)')
    plt.title('Burn Pulse')
    
    plot4 = plt.figure(4)
    
    annot_maxinteg(time, Yield)
    plt.plot(time, Yield, color='darkred')
    plt.xlim(0.8e-8, 0.84e-8)
    plt.xlabel('Time (s)')
    plt.ylabel('Total Yield (Neutrons)')
    plt.title('Integrated Burn Pulse')
    print('Maximum neutron yield', max(Yield))
    
    #Hotspot Diagnostics
    hotspot = []
    timesteps = []
    compressiontimes = []
    shapex = tion.shape[0]
    shapey = tion.shape[1]
    
    for i in range(0, shapex):
        if tion[i, 0]>= 2000:
            timesteps.append(i)
            for j in range(0, shapey):
                if tion[i, j]<2000:
                    hotspot.append(j)
                    compressiontimes.append(i)
                    break
                    
    mintimestep = min(timesteps)
    maxtimestep = max(timesteps)
    print('The range of timesteps over which the hotspot exceeds 2keV is', mintimestep, '-', maxtimestep)
    
    minimum = min(hotspot)
    print('The minimum hotspot radius is', minimum, 'um')
                    
    index = hotspot.index(minimum)
    peakcompressiontime = compressiontimes[index]
    print('The peak compression time index is', peakcompressiontime)

    peakdensity = rho[peakcompressiontime, minimum]
    print('The density at peak compression is', peakdensity, 'kg/m^3')

    hotspot = []
    timesteps = []
    shapex = tion.shape[0]
    shapey = tion.shape[1]
        
    for i in range(0, shapex):
        if tion[i, 0]>= 2000:
            timesteps.append(i)
            for j in range(0, shapey):
                if tion[i, j]<2000:
                    hotspot.append(j)
                    break
                
    mintimestep = min(timesteps)
    maxtimestep = max(timesteps)
    x = np.arange(mintimestep-1, maxtimestep, 1)
    plot5 = plt.figure(5)
   # plt.plot(x, hotspot, color= 'darkcyan')
    plt.xlabel('Time index')
    plt.ylabel('Hotspot Radius/um')
    plt.title('Hotspot Radius in the Regime where Hotspot Temperature exceeds 2keV')
                    
    shapey = rho.shape[1]
    r = np.arange(0, shapey, 1)
    dt_mass = np.trapz(4*np.pi*rho[0]*r*r*1e-12, dx=1e-6) #initial DT mass in kg
    
    Dmass = 3.3435e-27 #kg (i.e. atomic mass number (u)* kg conversion)
    Tmass = 5.0082711e-27 #kg
    
    initpairs = dt_mass/(Dmass + Tmass)
    yieldmax = Yield.max()
    
    burnupfraction = (yieldmax/initpairs)*100

    print('The percentage burn-up is', burnupfraction, '%')
    
    x = time.shape
    shapex = x[0]
    templist =[]
    for i in range(0, shapex):
        element = burnavTi[i] * yieldrate[i] * dt[i]
        templist.append(element)
        summation = sum(templist)
        totalyield = Yield.max()
        Teff = summation/totalyield 
        
    print('The effective temperature is', Teff, 'eV')
    
    plt.show()
    
    #Implosion Plots 
    plot5 = plt.figure(5)
    
    t1 = peakcompressiontime - 2
    t2 = peakcompressiontime - 1
    t3 = peakcompressiontime
    t4 = peakcompressiontime + 1
    t5 = peakcompressiontime + 2
    xrange = 100
    
    fig, axes= plt.subplots(nrows=2, ncols=5, sharex= 'row', sharey= 'row')
    ax1 = axes[0, 0]    #define axes 
    ax2 = ax1.twinx()
    ax3 = axes[0, 1]
    ax4 = ax3.twinx()
    ax5 = axes[0, 2]
    ax6 = ax5.twinx()
    ax7 = axes[0, 3]
    ax8 = ax7.twinx()
    ax9 = axes[0, 4]
    ax10 = ax9.twinx()
    ax11 = axes[1, 0]
    ax12 = axes[1, 1]
    ax13 = axes[1, 2]
    ax14 = axes[1, 3]
    ax15 = axes[1, 4]
    
    ax2.get_shared_y_axes().join(ax2, ax4, ax6, ax8, ax10)  #ensure axes share the same scale 
    
    ax1.set_ylabel('Density, ρ  x 10\N{SUPERSCRIPT FIVE} (kg/m\N{SUPERSCRIPT THREE})')
    ax1.plot(dist, density[t1],color='darkred')
    ax1.plot(dist, ccdens[t1], color='lightgreen')
    #ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))  #in case scientific labels are desired
    ax1.set_title('to-40ps')
    ax1.set_xlim([0, xrange])
    
    ax2.plot(dist, telec[t1], color='mediumseagreen', linestyle='dashed')
    ax2.plot(dist, tion[t1], color='teal')
    ax2.axes.yaxis.set_ticklabels([])
    ax2.axes.xaxis.set_ticklabels([])
    
    ax3.plot(dist, density[t2],color='darkred')
    ax3.plot(dist, ccdens[t2], color='lightgreen')
    ax3.set_title('to') 
    ax4.plot(dist, telec[t2], color='mediumseagreen', linestyle='dashed')
    ax4.plot(dist, tion[t2], color='teal')
    ax4.axes.yaxis.set_ticklabels([])
    
    ax5.plot(dist, density[t3],color='darkred')
    ax5.plot(dist, ccdens[t3], color='lightgreen')
    ax5.set_title('to+40ps')
    ax6.plot(dist, telec[t3], color='mediumseagreen', linestyle='dashed')
    ax6.plot(dist, tion[t3], color='teal')
    ax6.axes.yaxis.set_ticklabels([])
    
    ax7.plot(dist, density[t4],color='darkred')
    ax7.plot(dist, ccdens[t4], color='lightgreen')
    ax7.set_title('to+80ps')
    ax8.plot(dist, telec[t4], color='mediumseagreen', linestyle='dashed')
    ax8.plot(dist, tion[t4], color='teal')
    ax8.axes.yaxis.set_ticklabels([])
    
    ax9.plot(dist, density[t5],color='darkred', label='$ρ_{DT}$')
    ax9.plot(dist, ccdens[t5], color='lightgreen', label='$ρ_{CC}$')
    ax9.set_title('to+120ps')
    ax10.set_ylabel('Temperature (eV)') 
    ax10.plot(dist, telec[t5], color='mediumseagreen', linestyle='dashed', label='$T_{e}$')
    ax10.plot(dist, tion[t5], color='teal', label='$T_{i}$')
    #ax10.ticklabel_format(axis='y', style="sci", scilimits=(0,0))

    ax11.set_ylabel('Power Density (W/m\N{SUPERSCRIPT THREE})')
    ax11.plot(dist, qtot[t1],color='black')
    ax11.plot(dist, elecwork[t1], color='cornflowerblue')
    #ax11.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax11.plot(dist, deltot[t1]*scalefactor, color='sienna')
    ax11.plot(dist, pdv[t1], color='indigo')
    #ax11.axes.yaxis.set_ticklabels([])
    ax11.set_xlim([0, xrange])
    
    ax12.plot(dist, qtot[t2],color='black')
    ax12.plot(dist, elecwork[t2], color='cornflowerblue')
    ax12.plot(dist, deltot[t2]*scalefactor, color='sienna')
    ax12.plot(dist, pdv[t2], color='indigo')
    ax12.set_xlim([0, xrange])

    ax13.set_xlabel('Radial Distance (μm)')
    ax13.plot(dist, qtot[t3],color='black')
    ax13.plot(dist, elecwork[t3], color='cornflowerblue')
    ax13.plot(dist, deltot[t3]*scalefactor, color='sienna')
    ax13.plot(dist, pdv[t3], color='indigo')
    ax13.set_xlim([0, xrange])
    
    ax14.plot(dist, qtot[t4],color='black')
    ax14.plot(dist, elecwork[t4], color='cornflowerblue')
    ax14.plot(dist, deltot[t4]*scalefactor, color='sienna')
    ax14.plot(dist, pdv[t4], color='indigo')
    ax14.set_xlim([0, xrange])

    ax15.plot(dist, qtot[t5],color='black', label='$W_{α}$')
    ax15.plot(dist, elecwork[t5], color='cornflowerblue', label='$W_{e}$') 
    ax15.plot(dist, deltot[t5]*scalefactor, color='sienna', label='$W_{γ}$')
    ax15.plot(dist, pdv[t5], color='indigo', label='$W_{PdV}$')
    #ax20.ticklabel_format(axis='y', style="sci", scilimits=(0,0))
    ax15.set_xlim([0, xrange])
    
    plt.xlim([0,xrange])
    #fig.legend(loc=0)
    fig.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    fig.suptitle('Time evolution of 1D impolosion')
    plt.show()
    
    #Curve fit for Burn Pulse to find FWHM 
    plot6 = plt.figure(6)
    
    def gaussian(x, amp, mean, stddev):
        return ((amp / (np.sqrt(2*np.pi) * stddev))) * np.exp(-(x-mean)**2 / (2*stddev**2))
    
    ampguess = yieldrate.max()
    meanguess = time[np.argmax(yieldrate)]
    stdguess = 1e-9 #THIS CAN BE ADJUSTED IF THE CURVEFIT LOOKS INCORRECT
    xlim1 = 0.8e-8 #CHANGE ACCORDING TO INDIVIDUAL GRAPH
    xlim2 = 0.83e-8  #CHANGE ACCORDING TO INDIVIDUAL GRAPH
    
    ymax = yieldrate.max()
    halfymax = ymax/2
    plt.axhline(y=halfymax, color='black')
    plt.plot(time, yieldrate, color='navy')
    initialguess = [ampguess, meanguess, stdguess]
    popt, pcov = curve_fit(gaussian, time, yieldrate, initialguess)
    plt.xlim(xlim1, xlim2)
    
    xfit = np.arange(0, 2e-8, 2e-12) 
    plt.plot(xfit, gaussian(xfit, *popt), 'r', label = 'amp=%.3e, mean=%.3e, stddev=%.3e' % tuple(popt))
    plt.legend(loc='lower right', fontsize=8, handlelength = 0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Yield (neutrons/s)')
    plt.title('Curve Fit of Burn Pulse for FWHM')
    plt.show()
    print('Confinement Time (FWHM) =', tuple(popt)[2], 's')
    
    plot7 = plt.figure(7)
    
    sf = np.shape(Yield)[0]/np.shape(rho)[0]
    interval = (1e-8*SCALEFACTOR)/np.shape(Yield)[0]
    xmax = time[np.argmax(Yield)]
    int2 = xmax/interval
    t = int2/sf
    bangtime = round(t)
    print('Bang time index is', bangtime)
    plt.plot(dist , density[bangtime], color='black')
    plt.xlim(0, 300)
    plt.xlabel('Radius (μm)')
    plt.ylabel('Density, ρ x 10\N{SUPERSCRIPT FIVE} (kg/m\N{SUPERSCRIPT THREE})')
    plt.title('Density Profile at Bangtime')
    
    plot8 = plt.figure(8)
    
    plt.plot(time, burnavTi/1e3, color='navy')
    plt.xlabel('Time (s)')
    plt.ylabel('Burn-averaged Ion Temperature (keV)')
    plt.title('Burn-averaged Ion Temperature')
    plt.xlim(0.65e-8, 0.9e-8)
    
    #Hotspot-to-total DT mass ratio
    plot9 = plt.figure(9)
    
    def mass_ratio():    
        shapex = tion.shape[0] 
        shapey = tion.shape[1]
        r = np.linspace(0, shapey, shapey)
        ratio = []
        dt_mass = np.trapz(4*np.pi*rho[0, 0:shapey]*r*r*1e-12, dx=1e-6)
            
        for j in range(0, shapex):  #for loop to find upper limit of mass integral (where T_ion < 2keV)
            wherelist = np.where(tion>=2000) 
            zerolist = wherelist[0]
            zeroindex = zerolist[0] #lowest value for which Tion>=2000 i.e. hotspot radius
            if j<zeroindex:
                r_hs = 0
            else:
                r_hs = min(i for i in range(0, shapey) if tion[j, i] < 2000)
                    
            x = np.linspace(0, r_hs, r_hs)
            hs_mass = np.trapz(4*np.pi*rho[j, 0:r_hs]*x*x*1e-12, dx=1e-6) #mass integral
            mass_ratio = hs_mass/dt_mass
            ratio.append(mass_ratio)
                
        return ratio
        
    plt.plot(tarray, mass_ratio(), color='teal')
    plt.xlim(380, 430)
    plt.ylim(0, 1)
    plt.ylabel('Mhs/Mdt')
    plt.xlabel('Time Index')
    plt.title('Normalised Ratio of Hotspot and DT Mass')
    
        
#%% Power Density plots for five specified times- t1, t2, t3, t4 and t5

t1 = 495
t2 = 496
t3 = 497
t4 = 498
t5 = 499
xrange = 250
dist = np.arange(0, np.shape(rho)[1], 1)

def powerdens():
    fig = plt.figure(figsize = [6.5, 10])
    gs = fig.add_gridspec(5, 1, hspace = 0.35)
    (ax1), (ax2), (ax3), (ax4), (ax5) = gs.subplots()
    ax1.plot(dist, qtot[t1]/1e29, color='black')
    ax1.plot(dist, deltot[t1]*scalefactor/1e29, color='cornflowerblue')
    ax1.plot(dist, elecwork[t1]/1e29, color='sienna')
    ax1.plot(dist, pdv[t1]/1e29, color='indigo')
    ax1.set_xlim(0, xrange)
    ax1.grid()
    ax1.set_title('t = {}'.format(t1), x=0.5)
    ax1.axes.xaxis.set_ticklabels([])
    ax2.plot(dist, qtot[t2]/1e29, color='black')
    ax2.plot(dist, deltot[t2]*scalefactor/1e29, color='cornflowerblue')
    ax2.plot(dist, elecwork[t2]/1e29, color='sienna')
    ax2.plot(dist, pdv[t2]/1e29, color='indigo')
    ax2.set_xlim(0, xrange)
    ax2.grid()
    ax2.set_title('t = {}'.format(t2), x=0.5)
    ax2.axes.xaxis.set_ticklabels([])
    ax3.set_ylabel('Power Density x $10^{29}$ (W/m\N{SUPERSCRIPT THREE})')
    ax3.plot(dist, qtot[t3]/1e29, color='black')
    ax3.plot(dist, deltot[t3]*scalefactor/1e29, color='cornflowerblue')
    ax3.plot(dist, elecwork[t3]/1e29, color='sienna')
    ax3.plot(dist, pdv[t3]/1e29, color='indigo')
    ax3.set_xlim(0, xrange)
    ax3.grid()
    ax3.set_title('t = {}'.format(t3), x=0.5)
    ax3.axes.xaxis.set_ticklabels([])
    ax4.plot(dist, qtot[t4]/1e29, color='black')
    ax4.plot(dist, deltot[t4]*scalefactor/1e29, color='cornflowerblue')
    ax4.plot(dist, elecwork[t4]/1e29, color='sienna')
    ax4.plot(dist, pdv[t4]/1e29, color='indigo')
    ax4.set_xlim(0, xrange)
    ax4.grid()
    ax4.set_title('t = {}'.format(t4), x=0.5)
    ax4.axes.xaxis.set_ticklabels([])
    ax5.plot(dist, qtot[t5]/1e29, color='black', label='$W_{α}$')
    ax5.plot(dist, deltot[t5]*scalefactor/1e29, color='cornflowerblue', label='$W_{γ}$')
    ax5.plot(dist, elecwork[t5]/1e29, color='sienna', label='$W_{e}$')
    ax5.plot(dist, pdv[t5]/1e29, color='indigo', label='$W_{PdV}$')
    ax5.set_xlim(0, xrange)
    ax5.grid()
    ax5.set_title('t = {}'.format(t5), x=0.5)
    fig.legend(loc='lower right', bbox_to_anchor=(1.05, 0.76))
    
#%% Comparison of BUF and non-BUF S=1 yield 

fig, ax = plt.subplots()
ax.plot(time_1_buf, yield_1_buf/1e27, linestyle='--', label='BUF')
plt.plot(time_1, yield_1/1e27, linestyle = '--', label='No BUF')
plt.legend()
#plt.xlim(0.8e-8, 0.84e-8)
plt.xlabel('Time (s)')
plt.ylabel('Yield Rate (Neutron/s)')
plt.title('Burn Pulse with and without BUF')

#%% Comparison of BUF yields for all scale factors 

factor = 1e29

withbuf = [max(yield_09_buf)/factor, max(yield_1_buf)/factor, max(yield_11_buf)/factor, max(yield_12_buf)/factor, max(yield_13_buf)/factor, max(yield_14_buf)/factor]
withoutbuf = [max(yield_09)/factor, max(yield_1)/factor, max(yield_11)/factor, max(yield_12)/factor, max(yield_13)/factor, max(yield_14)/factor]

text_values = ["S=0.9", "S=1.0", "S=1.1", "S=1.2", "S=1.3", "S=1.4"]
x = np.arange(1, len(text_values) + 1, 1)
plt.scatter(x, withbuf, c='cornflowerblue', s=50, label='With BUF')
plt.grid(True)
plt.ylabel('Neutron Yield x $10^{29}$ /s')
plt.tick_params(labelcolor='black')
plt.scatter(x, withoutbuf, c='indianred', s=50, label='Without BUF')
plt.suptitle('Scale Factor Comparison of Neutron Yield with/without BUF', y=0.93)
plt.xticks(x, text_values)
#plt.ylim(4, 7)
plt.legend()

#%% Comparison of how much the yield is reduced for BUF vs non-BUF for all scale factors 

reductionfactor =[]

reductionfactor.append(max(yield_09_buf)/max(yield_09))
reductionfactor.append(max(yield_1_buf)/max(yield_1))
reductionfactor.append(max(yield_11_buf)/max(yield_11))
reductionfactor.append(max(yield_12_buf)/max(yield_12))
reductionfactor.append(max(yield_13_buf)/max(yield_13))
reductionfactor.append(max(yield_14_buf)/max(yield_14))

text_values = ["S=0.9", "S=1.0", "S=1.1", "S=1.2", "S=1.3", "S=1.4"]
x = np.arange(1, len(text_values) + 1, 1)
plt.scatter(x, reductionfactor, c='black', s=50)
plt.xticks(x, text_values, size=18)
plt.ylabel('Reduction factor with BUF')
plt.ylim(0.5, 1)

#%% Changing coast time by changing frequency dependent spectrum 

plt.plot(hydra_time, hydra_energy)
plt.xlabel('Time (μs)')
plt.ylabel('Energy (J)')
plt.title('Frequency Dependent Input Spectrum')

hydra_energy_trans = np.transpose(hydra_energy)

meanlist = []
for i in range(0, np.shape(hydra_energy_trans)[0], 1):
    peak = np.where(hydra_energy_trans[i] == max(hydra_energy_trans[i]))[0][0]
    meanlist.append(peak)
average = sum(meanlist)/len(meanlist)  
averagetime = hydra_time[round(average)]

array = []
for i in range(0, np.shape(hydra_time)[0]-1):
    diff = hydra_time[i+1] - hydra_time[i]
    array.append(diff)

x = np.arange(0, len(array), 1)
plt.plot(x, array)   
plt.xlabel('Radius (um)')
plt.ylabel('Difference in successive timesteps (us)') 

plt.plot(hydra_time, hydra_energy_trans[20])

y = hydra_energy_trans[20]
maximum = max(y)
plt.plot(hydra_time, hydra_energy_trans[20])
plt.axvline(x=0.0078256)

new_spectrum=[]
extension = 28

for i in range(0, np.shape(hydra_energy_trans)[0], 1):
    y = hydra_energy_trans[i]
    maximum = max(y)
    minrange = np.where(y==maximum)[0][0]
    maxrange = minrange + extension
    for j in range(minrange, maxrange):
        new_energy = y
        new_energy[j] = maximum
    new_spectrum.append(new_energy)

new_spectrum_trans = np.transpose(new_spectrum)
plt.plot(hydra_time, new_spectrum_trans)
plt.xlabel('Time (μs)')
plt.ylabel('Energy (J)')
plt.title('Frequency Dependent Input Spectrum')

textfile = open(r"C:\Users\misha\OneDrive\Desktop\YEAR 4\MSci Project\Coast Times\ext28.txt", "w")
for element in new_spectrum_trans:
    textfile.write(str(element))
textfile.close()

#%% Find burn front velocity after peak compression by using density spikes as the burn front 

peak_compression = 454
streak_freq = 20e-12

def velocity(rho):
    vel = []
    for i in range(peak_compression, np.shape(rho)[0], 1):
        peak = np.where(rho[i]==max(rho[i]))[0][0]
        vel.append(peak)
    print(np.shape(vel))
    difference = []
    for j in range(0, len(vel)-1, 1):
        difference.append((vel[j+1]-vel[j])/streak_freq)
    print(np.shape(difference))
    array = np.arange(0, len(difference), 1)
    plt.scatter(array, difference, s=5, c='black')

#%% Trajectory of a capsule in density-temperature space around bang time 

bang09 = 423
bang1 = 473
bang11 = 523
bang12 = 574
bang13 = 626
bang14 = 678

def trajectory():
    def averages(bangtime, rho, tion, col, lab):
        
        density = []
        temperature = []

        for i in range(0, np.shape(rho)[0]):
            density.append(sum(rho[i])/np.shape(rho)[1])

        for j in range(0, np.shape(tion)[0]):
            temperature.append(sum(tion[j])/np.shape(tion)[1])
        
        temp = np.array(temperature)/1e3
        dens = np.array(density)
        
        hotspot = []
        shapex = tion.shape[0]
        shapey = tion.shape[1]
        
        for k in range(0, shapex):
            if tion[k, 0]>= 2000:
                for l in range(0, shapey):
                    if tion[k, l]<2000:
                        hotspot.append(l)
                        break
            else:
                hotspot.append(0)
        
        rhor = np.multiply(dens, (hotspot))
        bang =  np.where(dens == max(dens))
       
        rangetemp = []
        rangetemp.append(temp[bang[0][0]-4:bang[0][0]+4])
        rangedens = []
        rangedens.append(rhor[bang[0][0]-4:bang[0][0]+4])
        
        plt.plot(rangedens[0], rangetemp[0], color = col, label= lab)
      #  plt.plot(dens,temp, color='yellow')
        plt.xlabel('Average Density (gm/cm\N{SUPERSCRIPT THREE})')
        plt.ylabel('Average Ion Temperature (keV)')
        plt.legend()
        plt.title('Areal Density vs Ion Temperature at Bangtime +/- 60ps')
        #plt.xlim(6.4, 7.5)
      #  plt.ylim(0.5, 0.8)
    
    averages(bang09, rho_09, tion_09, 'navy', 'S = 0.9')
    averages(bang1, rho_1, tion_1, 'black', 'S = 1')
    averages(bang11,  rho_11, tion_11, 'red', 'S = 1.1')
    averages(bang12, rho_12, tion_12, 'cornflowerblue', 'S = 1.2')
    averages(bang13, rho_13, tion_13, 'green', 'S = 1.3')
    averages(bang14, rho_14, tion_14, 'magenta', 'S = 1.4')

#%% Normalised integrated density vs normalised fusion plots 

def integralplots(DENSITY, FUSION, timestamp):
    rmax = np.shape(FUSION)[1]
    xrange = np.arange(0, rmax, 1)
    x = np.linspace(0, rmax, rmax)
    
    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111)
    lns1 = ax.plot(xrange, FUSION[timestamp], '--', color='navy', label = 'Fusion')
    ax2 = ax.twinx()
    lns2 = ax2.plot(xrange, DENSITY[timestamp]/1e3, '--', color='darkred', label = 'Density')
    plt.xlim(0, 600)
    ax.set_xlabel("Radius (μm)")
    ax.set_ylabel('Fusion Reaction Rate / s')
    ax2.set_ylabel('Density ($g/\,cm^{-3}$)')
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc ='upper left')
    
    fusionintegral = []
    
    for i in xrange:
        rfus = np.linspace(0, i, i)
        fintegral = np.trapz((4*np.pi*rfus*rfus*FUSION[timestamp,0:i]), dx=1e-6)
        fusionintegral.append(fintegral)
        
    totalfusion = np.trapz(4*np.pi*x*x*FUSION[timestamp, 0:rmax], dx=1e-6)
    normfusion = fusionintegral/totalfusion
    
    densityintegral = []
    
    for j in xrange:
        rdens = np.linspace(0, j, j)
        dintegral = np.trapz((4*np.pi*rdens*rdens*DENSITY[timestamp,0:j]), dx=1e-6)
        densityintegral.append(dintegral)
        
    totaldensity = np.trapz(4*np.pi*x*x*DENSITY[timestamp,0:rmax], dx=1e-6)
    normdensity = densityintegral/totaldensity
    
    fig2 = plt.figure(2)
    ax = fig2.add_subplot(111)
    ax.plot(xrange, normfusion, '--', color='navy', label = 'Fusion')
    ax2 = ax.twinx()
    ax2.plot(xrange, normdensity, '--', color='red', label = 'Density')
    plt.xlim(0, 600)
    ax.set_xlabel("Radius (μm)")
    ax.set_ylabel('Normalised Integrated Fusion Reaction Rate')
    ax2.set_ylabel('Normalised Integrated Density')
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)
    
    fig3 = plt.figure(3)
    plt.plot(normdensity, normfusion, color='black')
    plt.xlabel('Normalised Integrated Density')
    plt.ylabel('Normalised Integrated Fusion Reaction Rate')
    plt.axhline(y=0.98, color='rosybrown', linestyle='--') 
    plt.axvline(x=0.835, color='indianred', linestyle='--')
    #plt.title(timestamp)

    #plt.xlim(0.75, 0.86)

#%% Animation of integral plots

def integanim(DENSITY, FUSION, start, stop):
    rmax = np.shape(FUSION)[1]
    xrange = np.arange(0, rmax, 1)
    x = np.linspace(0, rmax, rmax)
    
    for k in range(start, stop):
        fig, ax = plt.subplots(1, 1)
        fusionintegral = []
        
        for i in xrange:
            rfus = np.linspace(0, i, i)
            fintegral = np.trapz((4*np.pi*rfus*rfus*FUSION[k,0:i]), dx=1e-6)
            fusionintegral.append(fintegral)
            
        totalfusion = np.trapz(4*np.pi*x*x*FUSION[k, 0:rmax], dx=1e-6)
        #plt.axvline(x=0.28, color='indianred', linestyle='--')
        plt.title(k)

#%% Normalised integrated density and fusion for scale factor comparisons at bang time 

def comp_values(DENSITY, FUSION, timestamp):
    rmax = np.shape(FUSION)[1]
    xrange = np.arange(0, rmax, 1)
    x = np.linspace(0, rmax, rmax)
    
    fusionintegral = []
    
    for i in xrange:
        rfus = np.linspace(0, i, i)
        fintegral = np.trapz((4*np.pi*rfus*rfus*FUSION[timestamp,0:i]), dx=1e-6)
        fusionintegral.append(fintegral)
        
    totalfusion = np.trapz(4*np.pi*x*x*FUSION[timestamp, 0:rmax], dx=1e-6)
    normfusion = fusionintegral/totalfusion
    
    densityintegral = []
    
    for j in xrange:
        rdens = np.linspace(0, j, j)
        dintegral = np.trapz((4*np.pi*rdens*rdens*DENSITY[timestamp,0:j]), dx=1e-6)
        densityintegral.append(dintegral)
        
    totaldensity = np.trapz(4*np.pi*x*x*DENSITY[timestamp,0:rmax], dx=1e-6)
    normdensity = densityintegral/totaldensity
    
    
    return normdensity, normfusion

plt.plot(comp_values(rho_09, fusion_09, 423)[0], comp_values(rho_09, fusion_09, 423)[1], color='brown', ls='dashdot', label='S=0.9')
plt.plot(comp_values(rho_1, fusion_1, 473)[0], comp_values(rho_1, fusion_1, 473)[1], color='midnightblue', ls='dashdot', label='S=1.0')
plt.plot(comp_values(rho_11, fusion_11, 523)[0], comp_values(rho_11, fusion_11, 523)[1], color='greenyellow', ls='dashdot', label='S=1.1')
plt.plot(comp_values(rho_12, fusion_12, 574)[0], comp_values(rho_12, fusion_12, 574)[1], color='mediumseagreen', ls='dashdot', label='S=1.2')
plt.plot(comp_values(rho_13, fusion_13, 626)[0], comp_values(rho_13, fusion_13, 626)[1], color='deeppink', ls='dashdot', label='S=1.3')
plt.plot(comp_values(rho_14, fusion_14, 678)[0], comp_values(rho_14, fusion_14, 678)[1], color='mediumorchid', ls='dashdot', label='S=1.4')
plt.xlabel('Normalised Capsule Mass')
plt.ylabel('Normalised Integrated Fusion Rate')
plt.title('Distribution of Fusion Reactions for Various Scale Factors at Bangtime')
plt.legend()
plt.axhline(y=0.98, color='black', linestyle='dotted') 

#%% Fusion colour plots using the imshow and colormesh methods 

plt.imshow(fusion_1,norm=colors.LogNorm(1, 1e16), aspect='auto', origin='lower') 
plt.colorbar(label= 'Temperature/eV')                                                  
plt.xlabel('Radial Distance (μm)')
plt.ylabel('Time (in intervals of 20ps)')
plt.title('Variation of Fusion Reactions and Capsule Radius with Time')

def intcolourplots(FUSION):
    rmax = np.shape(FUSION)[1]
    tmax = np.shape(FUSION)[0]
    xrange = np.arange(0, rmax, 1)
    x = np.linspace(0, rmax, rmax)
    
    l =  []

    for k in range(0, tmax):
        
        fusionintegral = []
        
        for i in xrange:
            rfus = np.linspace(0, i, i)
            fintegral = np.trapz((4*np.pi*rfus*rfus*FUSION[k,0:i]), dx=1e-6)
            fusionintegral.append(fintegral)
       
        totalfusion = np.trapz(4*np.pi*x*x*FUSION[k, 0:rmax], dx=1e-6)
        normfusion = fusionintegral/totalfusion
        where_are_NaNs = np.isnan(normfusion)
        normfusion[where_are_NaNs] = 0

        l.append(normfusion)
 
    fusionarray = np.array(l)
    print(np.shape(fusionarray)) 
    
    return fusionarray 

values = intcolourplots(fusion_1)    

x = np.linspace(0, 1, 1248)
y = np.linspace(0, 502*20e-12/1e-9, 502)
plt.pcolormesh(x, y, values)
plt.colorbar(label= 'Normalised Fusion Reaction Rate')
plt.ylabel('Time (ns)')
plt.xlabel('Normalised Capsule Mass')
#plt.ylim(454*20e-12/1e-9, )

#%% FWHM upper and lower bound 98% yield plots 

y09 = [26.5, 28, 29.5]
y1 = [78.5, 80.5, 82.5]
y11 = [84.4, 84, 83.5]
y12 = [84, 84, 83.3]
y13 = [83.5, 83.7, 83.9]
y14 = [83.5, 83.5, 83.7]

text_values = ["Lower FWHM of Burn Pulse", "Bangtime", "Upper FWHM of Burn Pulse"]
x = np.arange(1, len(text_values) + 1, 1)
plt.errorbar(x, y09, yerr = 0.1, color='brown', label = 'S=0.9', ls='dashdot', fmt='o')
plt.errorbar(x, y1, yerr = 0.1, color='midnightblue', label = 'S=1', ls='dashdot', fmt='o')
plt.errorbar(x, y11, yerr = 0.1, color='greenyellow', label = 'S=1.1', ls='dashdot', fmt='o')
plt.errorbar(x, y12, yerr = 0.05, color='mediumseagreen', label = 'S=1.2', ls='dashdot', fmt='o')
plt.errorbar(x, y13, yerr = 0.05, color='deeppink', label = 'S=1.3', ls='dashdot', fmt='o')
plt.errorbar(x, y14, yerr = 0.05, color='mediumorchid', label = 'S=1.4', ls='dashdot', fmt='o')
plt.xticks(x, text_values)
plt.ylabel('% of fuel in which 98% of fusion is occuring')
plt.title(' Scale Factors 1:1.4')
plt.legend(loc=(1.02, 0))
plt.ylim(78 , 86)
plt.show()

#%% Preheat off vs on for S=1

plt.plot(time_1pre/1e-9, yieldrate_1pre/1e27, label='Preheat off', color='navy')
plt.plot(time_1/1e-9, yieldrate_1/1e27, label='Preheat on', color='seagreen')
plt.xlim(8.8, 10)
plt.ylabel('Yield Rate x $10^{27}$/s')
plt.xlabel('Time (ns)')
plt.title('Burn Pulse for S=1')
plt.legend()

#%% Peak compression times for preheat on and off scale factors, and for different ice thicknesses 

peakcomp09 = 409
peakcomp1 = 454
peakcomp11 = 505
peakcomp12 = 556
peakcomp13 = 607
peakcomp14 = 659

#peakcomp09pre = 409
peakcomp1pre = 495
peakcomp11pre = 530
peakcomp12pre = 583
peakcomp13pre = 575 
peakcomp14pre = 623

peakcomp1p = 462
peakcomp2p = 465
peakcomp3p = 484
peakcomp4p = 488
peakcomp5p = 493
peakcomp6p = 492
peakcomp7p = 496
peakcomp8p = 500
peakcomp9p = 501
peakcomp10p = 501

#%% Scale factor comparison of diagnostics 

confinement = np.array([.36050696969726756, .3662869030383791, .22357293970475655, .19854275160871776, .1945160295619031,  .1990539295840084])
burnup = np.array([0.5262592677652878, 2.1741823572491166, 10.38399312434644, 18.262700348319928, 26.002076869709782, 31.831611894461236])
neutronyield = np.array([9.919e16, 5.656e17, 3.614e+18, 8.286e18, 1.473e19, 2.263e19])
efftemp = np.array([ 6332.029345188351, 8382.401010578855, 14352.959648039328, 19242.15542567321, 23524.528381317, 26751.703810987954])
peakcompdens = np.array([129832.0, 155532.0, 126022.0, 181660.0, 117180.0, 133576.0])

text_values = ["S=0.9", "S=1.0", "S=1.1", "S=1.2", "S=1.3", "S=1.4"]
x = np.arange(1, len(text_values) + 1, 1)
fig, ax = plt.subplots(3, 1, figsize=[14,14])
ax[0].scatter(x, confinement, c='cornflowerblue', s=100)
ax[0].grid(True)
ax[0].set_ylabel('Confinement Time (ns)', size=15)
ax[0].set_xticklabels([])
ax[0].tick_params(labelcolor='black', labelsize=14)
ax[1].scatter(x, burnup, c='indianred', s=100)
ax[1].grid(True)
ax[1].set_ylabel('Burnup Percentage (%)', size=15)
ax[1].set_ylim(0, 40)
ax[1].set_xticklabels([])
ax[1].tick_params(labelcolor='black', labelsize=14)
#ax[2].scatter(x, efftemp/1e3, c='seagreen', s=100)
#ax[2].grid(True)
#ax[2].set_ylabel('Effective Temperature (keV)', size=15)
#ax[2].set_xticklabels([])
#ax[2].tick_params(labelcolor='black', labelsize=14)
ax[2].scatter(x, neutronyield, c='indigo', s=100)
ax[2].set_ylabel('Neutron Yield', size=15)
ax[2].grid(True)
ax[2].set_yscale('log')
ax[2].tick_params(labelcolor='black', labelsize=14)
plt.suptitle('Scale Factor Comparison of Diagnostics', size=20, y=0.93)
plt.xticks(x, text_values, size=18)
fig.align_ylabels(ax[0:4])

#%% Preheat off scale factor comparison of diagnostics 

confinementpre = np.array([5.637916635825623e-11,  6.515678147733931e-11, 5.1960337539257616e-11, 1.5710501214354668e-10, 4.4983419453236167e-10])
burnuppre = np.array([0.6563580572130667, 1.0507758800517746, 3.772888544781614, 4.402214518311288e-07, 1.0745438118797126e-06])
neutronyieldpre = np.array([1.70761096e+17, 3.65690592e+17, 1.71178815e+18, 249381683000.0, 763791344000.0])
efftemppre = np.array([6107.163801392406, 6860.917044707083, 9100.571443242961, 2967.960994373824, 2528.632907527926 ])

text_values = ["S=1.0", "S=1.1", "S=1.2", "S=1.3", "S=1.4"]
x = np.arange(1, len(text_values) + 1, 1)
fig, ax = plt.subplots(4, 1, figsize=[14,14])
ax[0].scatter(x, confinementpre/1e-11, c='cornflowerblue', s=100)
ax[0].grid(True)
ax[0].set_ylabel('Confinement Time x $,10^{-11}$ s', size=15)
ax[0].set_xticklabels([])
ax[0].tick_params(labelcolor='black', labelsize=14)
ax[1].scatter(x, burnuppre, c='indianred', s=100)
ax[1].grid(True)
ax[1].set_ylabel('Burnup Percentage', size=15)
ax[1].set_ylim(0, 10)
ax[1].set_xticklabels([])
ax[1].tick_params(labelcolor='black', labelsize=14)
ax[2].scatter(x, efftemppre/1e3, c='seagreen', s=100)
ax[2].grid(True)
ax[2].set_ylabel('Effective Temperature (keV)', size=15)
ax[2].set_xticklabels([])
ax[2].tick_params(labelcolor='black', labelsize=14)
ax[3].scatter(x, neutronyieldpre, c='indigo', s=100)
ax[3].set_ylabel('Neutron Yield', size=15)
ax[3].grid(True)
ax[3].set_yscale('log')
ax[3].tick_params(labelcolor='black', labelsize=14)
plt.suptitle('Preheat-off Scale Factor Comparison of Diagnostics', size=20, y=0.93)
plt.xticks(x, text_values, size=18)
fig.align_ylabels(ax[0:4])

#%% Scale Factor Comparison of Diagnostics for varying Ice Thickness

confinementice = np.array([4.6552127103504e-11, 8.192941587470324e-11, 5.5595846203448083e-11 , 6.189628305097208e-11, 6.740014932406982e-11, 7.648607832814828e-11, 8.354158588989558e-11, 1.0433989946927491e-10, 1.8518648339998756e-10, 4.1311455485511927e-10])
burnupice = np.array([0.8297968815809181, 0.00037381392878975747, 0.09317051779774842,  0.056460083235717176, 0.0385045704273572, 0.027706145548766294, 0.021968017925823455, 0.016579549997634793,  0.012210755728153633, 0.007032217175844805])
neutronyieldice = np.array([2.50378881e+17, 8e16, 3.64922884e+16, 2.4853393e+16, 1.86809259e+16, 1.48393031e+16, 1.2792064e+16, 1.05222984e+16, 8341613540000000.0, 5186081740000000.0])
efftempice = np.array([6884.497293240667 , 2741.541636499615, 4674.691637812582,  4331.311135668235, 4044.224923082205, 3836.041316061816, 3708.733727519077, 3550.483857875987,  3383.6491403095515, 3340.20837980021])

text_values = ["1 ", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
x = np.arange(1, len(text_values) + 1, 1)
fig, ax = plt.subplots(4, 1, figsize=[14,14])
ax[0].scatter(x, confinementice/1e-11, c='cornflowerblue', s=100)
ax[0].grid(True)
ax[0].set_ylabel('Confinement Time x $,10^{-11}$ s', size=15)
ax[0].set_xticklabels([])
ax[0].tick_params(labelcolor='black', labelsize=14)
ax[1].scatter(x, burnupice, c='indianred', s=100)
ax[1].grid(True)
ax[1].set_ylabel('Burnup Percentage', size=15)
ax[1].set_ylim(0, 1)
ax[1].set_xticklabels([])
ax[1].tick_params(labelcolor='black', labelsize=14)
ax[2].scatter(x, efftempice/1e3, c='seagreen', s=100)
ax[2].grid(True)
ax[2].set_ylabel('Effective Temperature (keV)', size=15)
ax[2].set_xticklabels([])
ax[2].tick_params(labelcolor='black', labelsize=14)
ax[3].scatter(x, neutronyieldice, c='indigo', s=100)
ax[3].set_ylabel('Neutron Yield', size=15)
ax[3].grid(True)
ax[3].set_yscale('log')
ax[3].tick_params(labelcolor='black', labelsize=14)
ax[3].set_xlabel('Increase in Ice Thickness (%)', size=15)
plt.suptitle('Scale Factor Comparison of Diagnostics for varying Ice Thickness', size=20, y=0.93)
plt.xticks(x, text_values, size=18)
fig.align_ylabels(ax[0:4])


#%%

def timesteps(tion):
    tarray = np.arange(0, np.shape(tion)[0], 1)
    tarray1 = tarray*20e-12/1e-9
    return tarray1

def mass_ratio(tion, rho):    
    shapex = tion.shape[0] 
    shapey = tion.shape[1]
    r = np.linspace(0, shapey, shapey) 
    ratio = []
    dt_mass = np.trapz(4*np.pi*rho[0, 0:shapey]*r*r*1e-12, dx=1e-6)
    
    for j in range(0, shapex-1):  # ADDED TO FIX PREHEAT ERROR
        wherelist = np.where(tion>=2000) 
        zerolist = wherelist[0]
        zeroindex = zerolist[0] #lowest value for which Tion>=2000 i.e. hotspot radius
        if j<zeroindex:
            r_hs = 0
        else:
            r_hs = min(i for i in range(0, shapey) if tion[j, i] < 2000)
                    
        x = np.linspace(0, r_hs, r_hs)
        hs_mass = np.trapz(4*np.pi*rho[j, 0:r_hs]*x*x*1e-12, dx=1e-6) #mass integral
        mass_ratio = hs_mass/dt_mass
        ratio.append(mass_ratio)
    ratio.append(0)  #ADDED TO FIX PREHEAT ERROR
    return ratio

#%% Mass ratio comparison for all scale factors  

plt.plot(timesteps(tion_09)-peakcomp09, mass_ratio(tion_09, rho_09), color='brown', label='S=0.9', ls='dashdot')
plt.plot(timesteps(tion_1)-(peakcomp1*20e-12/1e-9), mass_ratio(tion_1, rho_1), color='black', label='S=1.0', ls='dashdot')
plt.plot(timesteps(tion_11)-peakcomp11, mass_ratio(tion_11, rho_11), color='orange', label='S=1.1', ls='dashdot')
plt.plot(timesteps(tion_12)-peakcomp12, mass_ratio(tion_12, rho_12), color='mediumseagreen', label='S=1.2', ls='dashdot')
plt.plot(timesteps(tion_13)-peakcomp13, mass_ratio(tion_13, rho_13), color='navy', label='S=1.3', ls='dashdot')
plt.plot(timesteps(tion_14)-peakcomp14, mass_ratio(tion_14, rho_14), color='mediumorchid', label='S=1.4', ls='dashdot')
plt.xlim(-0.5,0.6)
plt.ylim(0, 1)
plt.ylabel('Mhs/Mdt')
plt.xlabel('Time around Peak Compression t$_0$=0 (ns)')
plt.title(r'Ratio of Hotspot and DT Mass around Peak Compression t$_0$')
plt.legend(loc='upper right')
#print(np.where(mass_ratio(tion_1, rho_1) ==max(mass_ratio(tion_1, rho_1))))
#print((462*20e-12/1e-9)-(peakcomp1*20e-12/1e-9))
#print(max(mass_ratio(tion_09, rho_09)))
#print(np.where(mass_ratio(tion_09, rho_09) ==max(mass_ratio(tion_09, rho_09))))

#%% Burn-averaged temperature comparison for all scale factors 

plt.plot((time_09-(peakcomp09*20e-12))/1e-9, burnav_09/1e3, color='brown', label='S=0.9')
plt.plot((time_1-(peakcomp1*20e-12))/1e-9,burnav_1/1e3 , color='black', label='S=1.0')
plt.plot((time_11-(peakcomp11*20e-12))/1e-9, burnav_11/1e3, color='orange', label='S=1.1')
plt.plot((time_12-(peakcomp12*20e-12))/1e-9, burnav_12/1e3, color='mediumseagreen', label='S=1.2')
plt.plot((time_13-(peakcomp13*20e-12))/1e-9, burnav_13/1e3, color='navy', label='S=1.3')
plt.plot((time_14-(peakcomp14*20e-12))/1e-9, burnav_14/1e3, color='mediumorchid', label='S=1.4')
plt.xlim(-1.5,1)
plt.legend(loc='upper right')
plt.xlabel('Time around Peak Compression t$_0$=0 (ns)')
plt.ylabel('Burn-averaged Ion Temperature (keV)')
plt.title('Burn-averaged Ion Temperature for various Scale Factors')

#%% Burn pulse comparison for all scale factors 

plt.plot((time_09-(peakcomp09*20e-12))/1e-9, yieldrate_09, color='brown', label='S=0.9')
plt.plot((time_1-(peakcomp1*20e-12))/1e-9, yieldrate_1 , color='black', label='S=1.0')
plt.plot((time_11-(peakcomp11*20e-12))/1e-9, yieldrate_11, color='orange', label='S=1.1')
plt.plot((time_12-(peakcomp12*20e-12))/1e-9, yieldrate_12, color='mediumseagreen', label='S=1.2')
plt.plot((time_13-(peakcomp13*20e-12))/1e-9, yieldrate_13, color='navy', label='S=1.3')
plt.plot((time_14-(peakcomp14*20e-12))/1e-9, yieldrate_14, color='mediumorchid', label='S=1.4')
plt.xlim(-0.15,0.3)
#plt.xlim(-2, 1)
plt.yscale('log')
plt.ylim(5e25,1e30)
plt.legend(loc='upper right')
plt.xlabel('Time around Peak Compression t$_0$=0 (ns)')
plt.ylabel('Yield Rate (Neutron/s)')
plt.title('Burn Pulse for various Scale Factors')

#%% Preheat off mass ratio comparison for all scale factors 

#plt.plot(timesteps(tion_09pre)-peakcomp09pre, mass_ratio(tion_09, rho_09), color='brown', label='S=0.9', ls='dashdot')
plt.plot(timesteps(tion_1pre)-peakcomp1pre, mass_ratio(tion_1pre, rho_1pre), color='black', label='S=1.0', ls='dashdot')
plt.plot(timesteps(tion_11pre)-peakcomp11pre, mass_ratio(tion_11pre, rho_11pre), color='orange', label='S=1.1', ls='dashdot')
plt.plot(timesteps(tion_12pre)-peakcomp12pre, mass_ratio(tion_12pre, rho_12pre), color='mediumseagreen', label='S=1.2', ls='dashdot')
plt.plot(timesteps(tion_13pre)-peakcomp13pre, mass_ratio(tion_13pre, rho_13pre), color='navy', label='S=1.3', ls='dashdot')
plt.plot(timesteps(tion_14pre)-peakcomp14pre, mass_ratio(tion_14pre, rho_14pre), color='mediumorchid', label='S=1.4', ls='dashdot')
plt.xlim(-30, 30)
plt.ylim(0, 1)
plt.ylabel('Mhs/Mdt')
plt.xlabel('Time Index around Peak Compression t$_0$=0')
plt.title(r'Preheat-off Ratio of Hotspot and DT Mass around Peak Compression t$_0$')
plt.legend(loc='upper left')

#%% Preheat off burn-averaged temperature comparisons for all scale factors 

#plt.plot((time_09-(peakcomp09*20e-12))/1e-9, burnav_09/1e3, color='brown', label='S=0.9')
plt.plot((time_1pre-(peakcomp1pre*20e-12))/1e-9,burnav_1pre/1e3 , color='black', label='S=1.0')
plt.plot((time_11pre-(peakcomp11pre*20e-12))/1e-9, burnav_11pre/1e3, color='orange', label='S=1.1')
plt.plot((time_12pre-(peakcomp12pre*20e-12))/1e-9, burnav_12pre/1e3, color='mediumseagreen', label='S=1.2')
plt.plot((time_13pre-(peakcomp13pre*20e-12))/1e-9, burnav_13pre/1e3, color='navy', label='S=1.3')
plt.plot((time_14pre-(peakcomp14pre*20e-12))/1e-9, burnav_14pre/1e3, color='mediumorchid', label='S=1.4')
plt.xlim(-2,1)
plt.legend(loc='upper right')
plt.xlabel('Time around Peak Compression t$_0$=0 (ns)')
plt.ylabel('Burn-averaged Ion Temperature (keV)')
plt.title('Burn-averaged Ion Temperature for various Scale Factors')

#%% Preheat off burn pulse comparisons for all scale factors 

#plt.plot((time_09-(peakcomp09*20e-12))/1e-9, yieldrate_09, color='brown', label='S=0.9')
plt.plot((time_1pre-(peakcomp1pre*20e-12))/1e-9, yieldrate_1pre , color='black', label='S=1.0')
plt.plot((time_11pre-(peakcomp11pre*20e-12))/1e-9, yieldrate_11pre, color='orange', label='S=1.1')
plt.plot((time_12pre-(peakcomp12pre*20e-12))/1e-9, yieldrate_12pre, color='mediumseagreen', label='S=1.2')
plt.plot((time_13pre-(peakcomp13pre*20e-12))/1e-9, yieldrate_13pre, color='navy', label='S=1.3')
plt.plot((time_14pre-(peakcomp14pre*20e-12))/1e-9, yieldrate_14pre, color='mediumorchid', label='S=1.4')
#plt.xlim(-0.15,0.3)
plt.xlim(-2, 0.7)
plt.yscale('log')
plt.ylim(10e19,3e28)
plt.legend(loc='upper left')
plt.xlabel('Time around Peak Compression t$_0$=0 (ns)')
plt.ylabel('Yield Rate (Neutron/s)')
plt.title('Preheat-off Burn Pulse for various Scale Factors')

#%% Ice thickness variation mass ratio comparisons

plt.plot(timesteps(tion_1)-peakcomp1, mass_ratio(tion_1, rho_1), color='black', label='Original', ls='dashdot')
#plt.plot(timesteps(tion_2p)-peakcomp2p, mass_ratio(tion_2p, rho_2p), color='brown', label='2%', ls='dashdot')
plt.plot(timesteps(tion_4p)-peakcomp4p, mass_ratio(tion_4p, rho_4p), color='deeppink', label='4% increase', ls='dashdot')
plt.plot(timesteps(tion_6p)-peakcomp6p, mass_ratio(tion_6p, rho_6p), color='orange', label='6% increase', ls='dashdot')
plt.plot(timesteps(tion_8p)-peakcomp8p, mass_ratio(tion_8p, rho_8p), color='mediumseagreen', label='8% increase', ls='dashdot')
#plt.plot(timesteps(tion_10p)-peakcomp10p, mass_ratio(tion_10p, rho_10p), color='navy', label='10%', ls='dashdot')
plt.xlim(-30, 20)
plt.ylim(0, 1)
plt.ylabel('Mhs/Mdt')
plt.xlabel('Time Index around Peak Compression t$_0$=0')
plt.title(r'Ratio of Hotspot and DT Mass around Peak Compression t$_0$ for varying Ice Thickness')
plt.legend(loc='upper left')

#%% Ice thickness variation burn-averaged temperature comparisons

plt.plot((time_1-(peakcomp1*20e-12))/1e-9, burnav_1/1e3, color='black', label='0%')
plt.plot((time_2p-(peakcomp2p*20e-12))/1e-9, burnav_2p/1e3, color='brown', label='2%')
plt.plot((time_4p-(peakcomp4p*20e-12))/1e-9,burnav_4p/1e3 , color='deeppink', label='4%')
plt.plot((time_6p-(peakcomp6p*20e-12))/1e-9, burnav_6p/1e3, color='orange', label='6%')
#plt.plot((time_8p-(peakcomp8p*20e-12))/1e-9, burnav_8p/1e3, color='mediumseagreen', label='8%')
plt.plot((time_10p-(peakcomp10p*20e-12))/1e-9, burnav_10p/1e3, color='navy', label='10%')
plt.xlim(-2,1)
plt.legend(loc='upper right')
plt.xlabel('Time around Peak Compression t$_0$=0 (ns)')
plt.ylabel('Burn-averaged Ion Temperature (keV)')
plt.title('Burn-averaged Ion Temperature for various Ice Thicknesses')

#%% Ice thickness variation burn pulse comparisons

plt.plot((time_1-(peakcomp1*20e-12))/1e-9, yieldrate_1 , color='black', label='S=1.0')
plt.plot((time_2p-(peakcomp2p*20e-12))/1e-9, yieldrate_2p , color='brown', label='2%')
plt.plot((time_4p-(peakcomp4p*20e-12))/1e-9, yieldrate_4p, color='deeppink', label='4%')
plt.plot((time_6p-(peakcomp6p*20e-12))/1e-9, yieldrate_6p, color='orange', label='6%')
#plt.plot((time_8p-(peakcomp8p*20e-12))/1e-9, yieldrate_8p, color='mediumseagreen', label='8%')
plt.plot((time_10p-(peakcomp10p*20e-12))/1e-9, yieldrate_10p, color='navy', label='10%')
#plt.xlim(-0.15,0.3)
plt.xlim(-2, 0.7)
plt.yscale('log')
plt.ylim(10e19,3e28)
plt.legend(loc='upper left')
plt.xlabel('Time around Peak Compression t$_0$=0 (ns)')
plt.ylabel('Yield Rate (Neutron/s)')
plt.title('Burn Pulse for various Ice Thicknesses')

#%% BUF mass ratio comparisons for all scale factors 

plt.plot(timesteps(tion_09BUF)-peakcomp09, mass_ratio(tion_09BUF, rho_09BUF), color='brown', label='S=0.9', ls='dashdot')
plt.plot(timesteps(tion_1BUF)-peakcomp1, mass_ratio(tion_1BUF, rho_1BUF), color='black', label='S=1.0', ls='dashdot')
plt.plot(timesteps(tion_11BUF)-peakcomp11, mass_ratio(tion_11BUF, rho_11BUF), color='orange', label='S=1.1', ls='dashdot')
plt.plot(timesteps(tion_12BUF)-peakcomp12, mass_ratio(tion_12BUF, rho_12BUF), color='mediumseagreen', label='S=1.2', ls='dashdot')
plt.plot(timesteps(tion_13BUF)-peakcomp13, mass_ratio(tion_13BUF, rho_13BUF), color='navy', label='S=1.3', ls='dashdot')
plt.plot(timesteps(tion_14BUF)-peakcomp14, mass_ratio(tion_14BUF, rho_14BUF), color='mediumorchid', label='S=1.4', ls='dashdot')
plt.xlim(-20, 30)
plt.ylim(0, 1)
plt.ylabel('Mhs/Mdt')
plt.xlabel('Time Index around Peak Compression t$_0$=0')
plt.title(r'Ratio of Hotspot and DT Mass around Peak Compression t$_0$')
plt.legend(loc='upper right')

#%% Mass ablation rates for all scale factors 

def massplots():
    def mass_ratio(T_ion, rho):  
        shapex = T_ion.shape[0] 
        shapey = T_ion.shape[1]
        r = np.arange(0, shapey, 1)
        hsmass = []
        mass_diff = []
        xaxis = []
                
        for j in range(0, shapex):  #for loop to find upper limit of mass integral (where T_ion < 1keV)
            wherelist = np.where(T_ion>=1000) 
            zerolist = wherelist[0]
            zeroindex = zerolist[0] #lowest value for which Tion>=1000 i.e. hotspot radius
            if j<zeroindex:
                r_hs = 0
            else:
                r_hs = min(i for i in range(0, shapey) if T_ion[j, i] < 1000)
                        
            x = np.linspace(0, r_hs, r_hs)
            hs_mass = np.trapz(4*np.pi*rho[j, 0:r_hs]*x*x*1e-12, dx=1e-6) #mass integral
            hsmass.append(hs_mass)
        print(np.shape(hsmass))
    
        for k in range(1, np.shape(hsmass)[0]):
            massdiff = (hsmass[k] - hsmass[k-1])
            mass_diff.append(massdiff)
            xaxis.append(k)
        
        return xaxis, mass_diff
    
    x1 = mass_ratio(tion_1, rho_1)[0]
    mass1 = mass_ratio(tion_1, rho_1)[1]
    x11 = mass_ratio(tion_11, rho_11)[0]
    mass11 = mass_ratio(tion_11, rho_11)[1]
    x12 = mass_ratio(tion_12, rho_12)[0]
    mass12 = mass_ratio(tion_12, rho_12)[1]
    x13 = mass_ratio(tion_13, rho_13)[0]
    mass13 = mass_ratio(tion_13, rho_13)[1]
    x14 = mass_ratio(tion_14, rho_14)[0]
    mass14 = mass_ratio(tion_14, rho_14)[1]
    x09 = mass_ratio(tion_09, rho_09)[0]
    mass09 = mass_ratio(tion_09, rho_09)[1]
    
    plt.plot(x1, mass1, label='S=1')
    plt.plot(x11, mass11, label ='S=1.1')
    plt.plot(x12, mass12, label='S=1.2')
    plt.plot(x13, mass13, label='S=1.3')
    plt.plot(x14, mass14, label = 'S=1.4')
    plt.plot(x09, mass09, label='S=0.9')
    plt.xlim(390, 650)
    plt.title('Mass ablation rate for different scale factors')
    plt.xlabel('Timestep')
    plt.ylabel('Mass Ablation Rate (kg/ps')
    plt.legend()
    
#%% PVGamma calculation

hotspot = []
shapex = T_ion.shape[0]
shapey = T_ion.shape[1]
PVgamma = []
x = np.arange(0, shapex, 1)

for i in range(0, shapex):
    if T_ion[i, 0]>= 2000:
        for j in range(0, shapey):
            if T_ion[i, j]<2000:
                hotspot.append(j)
                break
    else:
        hotspot.append(0)

for i in range(0, shapex):  
    hotspotrad = hotspot[i]
    r = np.linspace(0, hotspotrad, hotspotrad)
    PVgam = np.trapz((((P_ion[i,0:hotspotrad])**3/5)*(4*np.pi*r*r*1e-12)), dx=1e-6)
    PVgamma.append(PVgam)
    
#%% Sound speed and burn front speed diagnostics 

instance = 490
peakline = 454
Time = data[0]
rmax=np.shape(tion)[1]
tmax=len(Time)
r=np.linspace(0, rmax, rmax)
peak_compression = 456
timestep = 1e-8/502
gamma = 5/3 #adiabatic ratio
k = 1.6e-19 #electron charge if T is in eV (otherwise use Boltzmann constant)
Z = 1 #no. of electrons per ion
Ad = 2 #atomic mass of D 
At = 3 #atomic mass of T
Mp = 1.67e-27 #mass of proton 
inweightav = (Ad+At)/(Ad*At) #inverse of weighted average
shapex = Time.shape[0]
kb = 1.38064852e-23

def soundspeedcheck():
    soundspeed = []
    for j in range(0, tion.shape[1]):
        c = np.sqrt((gamma * k * tion[instance, j] * inweightav)/Mp)
        soundspeed.append(c)
                
  #  print('Sound Speed List', soundspeed)
    plt.plot(r, soundspeed)
    plt.title(instance)
    plt.xlim(500, 610)
   # plt.vlines(peakline,0, 500000, color='black')
    speed = np.array(soundspeed)
    xneg = find_peaks(-speed)[0] #this used to be soundspeed i changed it to velocity
    xpos = find_peaks(speed)[0]
    print('i', instance)
    print('xpos', xpos)
    print('xneg', xneg)

def animatevelplot(start, end, interval, sleeptime):
    for i in range(start, end, interval):
        fig, ax = plt.subplots(1, 1)
        ax.set_ylabel('Velocity (μm/s)')
        plt.xlabel('Radial Distance (μm)')
        plt.title(i)
       # plt.plot(range(xlim))
        plt.xlim(0)
        ax.plot(r, velocity[i], 'r')
        time.sleep(sleeptime)
        
def velplot():
    #plt.xlim(517, 530)
    #plt.ylim(-300, 0)
    plt.plot(r,velocity[instance])
    plt.show()

def velocitycheck():
    print('i', instance)
    plt.plot(r,velocity[instance])
    plt.title(instance)
    xneg = find_peaks(-velocity[instance])[0]
    xneg1 = xneg.tolist()
    xpos = find_peaks(velocity[instance])[0]
    xpos2 = xpos.tolist()
    peaks = xneg1 + xpos2
    peaks.sort()
    print(peaks)
    #plt.vlines(peakline,-300000,100000, color='black')
    
def soundspeed(T):
    c = np.sqrt((gamma * k * T * inweightav)/Mp)
    return c

def Teff():
    templist =[]
    for i in range(0, shapex):
        element = burnavTi[i] * yieldrate[i] * dt[i]
        templist.append(element)
    summation = sum(templist)
    totalyield = yielddata.max()
    Teff = summation/totalyield 
    return Teff

def burnavTlist():
    burnav = []
    for i in range(0,shapex):
        element = soundspeed(burnavTi[i])
        burnav.append(element)
    return burnav

Tthresh = 1.723e3 #threshold temperature for fusion 
Teffect = Teff()

def shockfront():
    shockpoint = []
    speeds_mat = []
    for i in range(0, peak_compression):
        shockpoint.append(0)
    for i in range(peak_compression, 501): #tmax
        def lastelem():
            if len(shockpoint)==0:
                last_elem = 0
            else:
                last_elem = shockpoint[-1]
            return last_elem
        speed = []
        peaks = []
        corr_speed = []
        for j in range(0, tion.shape[1]):
            c = np.sqrt((gamma * k * tion[i,j] * inweightav)/Mp)
            speed.append(c)
        soundspeed = np.array(speed)
        xneg = find_peaks(-soundspeed)[0] 
        xpos = find_peaks(soundspeed)[0]
        for l in range(0, len(xpos)):
            if xpos[l] >= lastelem():
                peaks.append(xpos[l])
        for m in range(0, len(xneg)):   #is this even needed?
            if xneg[m] >= lastelem():
                peaks.append(xneg[m])
        peaks.sort()
        for n in range(0, len(peaks)):
            corr_speed.append(speed[peaks[n]])
        maximum = max(corr_speed)
        speeds_mat.append(maximum)
        index = corr_speed.index(maximum)
        shockpoint.append(peaks[index])
    return shockpoint 

def materialspeed():
    shockpoint = []
    speeds_mat = []
    t=np.arange(0,502,1)
    for i in range(0, peak_compression):
        shockpoint.append(0)
        speeds_mat.append(0)
    for i in range(peak_compression, 501):
        def lastelem():
            if len(shockpoint)==0:
                last_elem = 0
            else:
                last_elem = shockpoint[-1]
            return last_elem
        speed = []
        peaks = []
        corr_speed = []
        for j in range(0, tion.shape[1]):
            c = np.sqrt((gamma * k * tion[i,j] * inweightav)/Mp)
            speed.append(c)
        soundspeed = np.array(speed)
        xneg = find_peaks(-soundspeed)[0]
        xpos = find_peaks(soundspeed)[0]
        for l in range(0, len(xpos)):
            if xpos[l] >= lastelem():
                peaks.append(xpos[l])
        for m in range(0, len(xneg)):
            if xneg[m] >= lastelem():
                peaks.append(xneg[m])
        peaks.sort()
        for n in range(0, len(peaks)):
            corr_speed.append(speed[peaks[n]])
        maximum = max(corr_speed)
        speeds_mat.append(maximum)
        index = corr_speed.index(maximum)
        shockpoint.append(peaks[index])
    return speeds_mat

def plotshockfront():
    plt.ylim(0.9e-8,1e-8)
    t=np.arange(0,tmax,1)
    sf=shockfront(tion) 
    plt.hlines((peak_compression*1.99e-11),0,700, color='grey', label='Peak Compression')
    plt.plot(sf,Time, 'navy', label='Shock Radius')
    plt.xlabel('Shock Radius/um')
    plt.ylabel('Timestep') 
    plt.legend()
    plt.title('Shock front after peak compression')
    
def velocities():
    fig, ax = plt.subplots()
    shockfrontlist = shockfront()
    materialspeeds = materialspeed()
    burnavspeed = burnavTlist()
    shockvelocity = []
    t=np.arange(0,(1e-8-timestep),timestep)
    t_mat = np.linspace(0, 1e-8, 501)
    for i in range(0, len(shockfrontlist)):
        if i<len(shockfrontlist)-1:
            difference = (shockfrontlist[i+1]-shockfrontlist[i])*1e-6
            velocity = difference/timestep 
            shockvelocity.append(velocity)
    shockvelocity.append(shockvelocity[-1])
    t_adj = t - (peak_compression* timestep)
    t_mat_adj = t_mat - (peak_compression * timestep)
    plt.xlim(0, 0.9)
    plt.ylim(0, 1.5e6)
    plt.scatter(t_adj/1e-9, shockvelocity, color='black', marker='.', label='Burn Front Velocity')
    plt.plot(t_mat_adj/1e-9, materialspeeds, color='cornflowerblue', label='Sound Speed at Burn Front')
   # plt.hlines(soundspeed(Teffect),0,9e-10, color='grey',linestyles='dotted', label='Average Sound Speed')
    #plt.plot(tarray, burnavspeed, color='cyan', label='BurnAvTi Sound Speed')
    plt.xlabel('Time after Peak Compression (ns)')
    plt.ylabel('Velocity (m/s)')
   # plt.title('Shock Velocity after Peak Compression')
    plt.legend()
    #print(shockvelocity)

def velocities2():
    fig, ax = plt.subplots()
    shockfrontlist = shockfront()[0]
    materialindex = shockfront()[1]
    burnavspeed = burnavTlist()
    shockvelocity = []
    materialspeeds = []
    t=np.arange(0,(1e-8-timestep),timestep)
    t_mat = np.linspace(0, 1e-8, 502)
    for i in range(0, len(shockfrontlist)):
        if i<len(shockfrontlist)-1:
            difference = (shockfrontlist[i+1]-shockfrontlist[i])*1e-6
            vel = difference/timestep 
            shockvelocity.append(vel)
    shockvelocity.append(shockvelocity[-1])
    for j in range(1, len(materialindex)):
        materialspeeds.append(velocity[j-1, shockfrontlist[j]]) #for every timestep, find the velocity at one timestep before the peak
    difference = []
    zip_object = zip(shockvelocity, materialspeeds)
    for shockvelocity_i,materialspeeds_i in zip_object:
        difference.append(shockvelocity_i-materialspeeds_i)
    t_adj = t - (peak_compression* timestep)
    t_mat_adj = t_mat - (peak_compression * timestep)
    index = [0, -1]
    new_t_mat = np.delete(t_mat_adj, index)
    
    plot1 = plt.figure(1)
    ax.ticklabel_format(axis='y', style="sci", scilimits=(6,6))
    ax.ticklabel_format(axis='x', style="sci", scilimits=(-9,-9))
    plt.scatter(t_adj, shockvelocity, color='teal', marker='.', label='Unadjusted Burn Front Velocity')
    plt.scatter(new_t_mat, difference, color = 'midnightblue', marker ='.', label='Adjusted Burn Front Velocity')
    plt.plot(new_t_mat, materialspeeds, color = 'dimgrey',label='Background Velocity ')
    plt.xlim(0, 0.9e-9)
    plt.ylim(-0.5e6, 2e6)
    plt.xlabel('Time after Peak Compression (s)')
    plt.ylabel('Velocity (m/s)')
   # plt.title('Shock Velocity after Peak Compression')
    plt.grid()
    plt.legend()
    
    plot2 = plt.figure(2)
    materialindex.pop(-1)
    vellist = []
    for i in range(0, 500):
        speed = []
        for j in range(1, tion.shape[1]): #might need to start from 1
            c = np.sqrt(((kb*telec[i,j]) + (3*kb*tion[i,j]))*inweightav/Mp)
           # c = np.sqrt((gamma * k * tion[i,j] * inweightav)/Mp)
           #c = np.sqrt(gamma * pressure[i, j]/rho[i,j])
           # c = np.sqrt(gamma * k * tion[i,j] * inweightav/Mp)
            speed.append(c) 
       # print(speed)
        soundspeed = np.array(speed)
        vellist.append(soundspeed[materialindex[i]])
    #print(vellist[456:502])
    ax.ticklabel_format(axis='y', style="sci", scilimits=(6,6))
    ax.ticklabel_format(axis='x', style="sci", scilimits=(-9,-9))
   # plt.scatter(new_t_mat, difference, color = 'midnightblue', marker ='.', label='Adjusted Shock Velocity')
    plt.plot(new_t_mat, vellist, label='Sound Speed ahead of Shock Front')
    plt.xlabel('Time after Peak Compression (s)')
    plt.ylabel('Velocity (m/s)')
    plt.title('Shock Velocity after Peak Compression')
    plt.xlim(0, 0.9e-9)
   # plt.ylim(0, 1.2e12)
   # plt.ylim(0, 2e6)
    plt.grid()
    plt.legend()

#Animate plots of the soundspeed and temperature to see if the soundspeed peaks are accurate
def check1():
    for i in range(peak_compression, tmax-1):
        speed = []
        for j in range(0, tion.shape[1]):
            c = np.sqrt((gamma * k * tion[i,j] * inweightav)/Mp)
            speed.append(c) 
        soundspeed = np.array(speed)
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Radius (um)')
        ax1.set_ylabel('Sound Speed x 10\N{SUPERSCRIPT FIVE} (m/s)', color='black')
        ax1.plot(r, soundspeed/1e5, color='darkslategrey')
        ax1.tick_params(axis='y', labelcolor='darkslategrey')
        ax2 = ax1.twinx()
        ax2.set_ylabel('Ion Temperature (eV)', color='black') 
        ax2.plot(r, tion[i], color='indianred', dashes=(1, 2, 1, 2) , dash_joinstyle='round')
        ax2.tick_params(axis='y', labelcolor='indianred')
        fig.tight_layout()  
        plt.title(i)
        plt.show()
        
#Plot the soundspeed itself (so that we can zoom in to see behaviour of peak at a certain time)
def check2():
    instance = 480
    speed = []
    for j in range(0, tion.shape[1]):
        c = np.sqrt((gamma * k * tion[instance,j] * inweightav)/Mp)
        speed.append(c) 
    soundspeed = np.array(speed)
    plt.plot(r, soundspeed)
    plt.xlim(400, 440)

#Animate the soundspeed vs the material speed 
def check3():
    for i in range(peak_compression, tmax-1):
        speed = []
        for j in range(0, tion.shape[1]):
            c = np.sqrt((gamma * k * tion[i,j] * inweightav)/Mp)
            speed.append(c) 
        soundspeed = np.array(speed)
        plt.xlabel('Radius (um)')
        plt.ylabel('Speed x 10\N{SUPERSCRIPT FIVE} (m/s)', color='black')
        plt.plot(r, soundspeed/1e5, color='darkslategrey', label='Sound Speed')
        plt.plot(r, velocity[i]/1e5, color='indianred', dashes=(1, 2, 1, 2) , dash_joinstyle='round', label='Material Speed') 
        toround = (i-peak_compression)*20e-12/1e-9
        plt.title('t = $t_0$ + {} ns'.format(round(toround, 2), x=0.5))
        #plt.title(i)
        plt.legend()
       # plt.xlim(400, 420)
        plt.show()

#%% Initial toy problem density and temperature configuration

fig1 = plt.figure(1)
plt.plot(r, T_ion[0]/1e3, 'black', label='0ps')
plt.ylabel('Temperature (keV)')
plt.xlabel('Radius (μm)')
plt.xlim(0, 250)
#plt.axvline(x=75)

fig2 = plt.figure(2)
plt.plot(r, rho_DT[0]/1e3, 'black', label='0ps')
plt.xlabel('Radius (μm)')
plt.ylabel('Density ρ (g/cm\N{SUPERSCRIPT THREE})')
plt.xlim(0, 250)
#plt.axvline(x=75)

#%% Density, temperature, pressure and power density animation

start = 450
end = 500
interval = 1
xmin = 0
xrange = 800
dist = np.arange(0, np.shape(rho_DT)[1], 1)
ps_per_index = 20

def animatetoy_I():
    for i in range(start, end, interval):
        fig = plt.figure(figsize = [6.7252, 10])
        gs = fig.add_gridspec(4, 1, hspace = 0.35)
        (ax1), (ax2), (ax3), (ax4) = gs.subplots()
        ax1.set_xlim(xmin, xrange)
        ax1.set_ylabel('Density ρ (g/cm\N{SUPERSCRIPT THREE})')
        ax1.plot(dist, rho_DT[i]/1e3,color='darkred') #1e3 kgm^-3 in one gcm^-3
        ax1.axes.xaxis.set_ticklabels([])
        actualtime = i * 1
        ax1.set_title('t = {} ps'.format(actualtime), x=0.49)
        ax2.set_xlim(xmin, xrange)
        ax2.set_ylabel('Temperature T (keV)')
        ax2.plot(dist, T_ion[i]/1e3, color='rosybrown')
        ax2.axes.xaxis.set_ticklabels([])
        ax3.set_xlim(xmin, xrange)
        ax3.set_ylabel('Pressure P x $10^{18}$ (Gbar)') #1e5 N/m^2 in a bar, 1e9 bars in a Gbar
        ax3.plot(dist, pres[i]/1e14, color='navy')
        ax3.axes.xaxis.set_ticklabels([])
        ax4.set_xlim(xmin, xrange)
        ax4.set_ylabel('Power Density x $10^{29}$ (W/m\N{SUPERSCRIPT THREE})')
        ax4.plot(dist, alpha[i]/1e29, color='black', label='$W_{α}$')
        ax4.plot(dist, rad[i]*scalefactor/1e29, color='cornflowerblue', label='$W_{γ}$')
        ax4.plot(dist, elec_cond[i]/1e29, color='sienna', label='$W_{e}$')
        ax4.plot(dist, pdv[i]/1e29, color='indigo', label='$W_{PdV}$')
        ax4.set_xlabel('Radius (μm)')
        fig.legend(loc='lower right', bbox_to_anchor=(1.08, 0.13))
        fig.suptitle('Index = {}'.format(i), y=0.93)
        fig.align_ylabels([ax1, ax2, ax3, ax4])
       # plt.ylim(-10, 22)
        
        plt.show()

#%% Temperature, density and pressure evolution plots 

t1 = 0
t2 = 50
t3 = 60
t4 = 70
t5 = 80
xrange = 400
yrange= 10000
dist = np.arange(0, np.shape(rho_DT)[1], 1)
ps_per_index = 1

def tpd(): #tpd = temperature, density and pressure
    fig = plt.figure(figsize = [6.5, 10])
    gs = fig.add_gridspec(3, 1, hspace = 0.35)
    (ax1), (ax2), (ax3) = gs.subplots()
    ax1.set_ylabel('Density ρ (g/cm\N{SUPERSCRIPT THREE})')
    ax1.plot(dist, rho_DT[t1]/1e3,color='black', label='t={} ps'.format(t1*ps_per_index))
    ax1.plot(dist, rho_DT[t2]/1e3,color='navy',linestyle='dashed', markersize=3 ,label='t={} ps'.format(t2*ps_per_index))
    ax1.plot(dist, rho_DT[t3]/1e3,color='cornflowerblue',linestyle='dashed', markersize=3, label='t={} ps'.format(t3*ps_per_index))
    ax1.plot(dist, rho_DT[t4]/1e3,color='darkolivegreen', linestyle='dashed', markersize=3, label='t={} ps'.format(t4*ps_per_index))
    ax1.plot(dist, rho_DT[t5]/1e3,color='indigo', linestyle='dashed', markersize=3, label='t={} ps'.format(t5*ps_per_index))
   # ax1.plot(dist, rho_DT[t1]/1e3,color='black', label='t= $t_0$')
   # ax1.plot(dist, rho_DT[t2]/1e3,color='navy',linestyle='dashed', markersize=3 ,label='t=$t_0$ + 40ps')
   # ax1.plot(dist, rho_DT[t3]/1e3,color='cornflowerblue',linestyle='dashed', markersize=3, label='t=$t_0$ + 80ps')
   # ax1.plot(dist, rho_DT[t4]/1e3,color='darkolivegreen', linestyle='dashed', markersize=3, label='t=$t_0$ + 120ps')
   # ax1.plot(dist, rho_DT[t5]/1e3,color='indigo', linestyle='dashed', markersize=3, label='t=$t_0$ + 160ps')
    ax1.set_yscale('log')
    ax1.set_xlim(0, xrange)
    ax1.set_ylim(1,)
    ax1.axes.xaxis.set_ticklabels([])
    ax2.set_ylabel('Temperature T (keV)')
    ax2.plot(dist, T_ion[t1]/1e3, color='black')
    ax2.plot(dist, T_ion[t2]/1e3, color='navy',linestyle='dashed', markersize=3 )
    ax2.plot(dist, T_ion[t3]/1e3, color='cornflowerblue',linestyle='dashed', markersize=3)
    ax2.plot(dist, T_ion[t4]/1e3, color='darkolivegreen', linestyle='dashed', markersize=3)
    ax2.plot(dist, T_ion[t5]/1e3, color='indigo', linestyle='dashed', markersize=3)
    ax2.set_yscale('log')
    ax2.set_xlim(0, xrange)
    ax2.set_ylim(1e-1,)
    ax2.axes.xaxis.set_ticklabels([])
    ax3.set_ylabel('Pressure P (Gbar)') #1e5 N/m^2 in a bar, 1e9 bars in a Gbar
    ax3.plot(dist, pres[t1]/1e14, color='black')
    ax3.plot(dist, pres[t2]/1e14, color='navy',linestyle='dashed', markersize=3)
    ax3.plot(dist, pres[t3]/1e14, color='cornflowerblue',linestyle='dashed', markersize=3)
    ax3.plot(dist, pres[t4]/1e14, color='darkolivegreen', linestyle='dashed', markersize=3)
    ax3.plot(dist, pres[t5]/1e14,color='indigo', linestyle='dashed', markersize=3)
    ax3.set_yscale('log')
    ax3.set_xlim(0, xrange)
    ax3.set_ylim(1,)
    ax3.set_xlabel('Radius (μm)')
    fig.align_ylabels([ax1, ax2, ax3])
    #fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.13))
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.88))
    #fig.legend(bbox_to_anchor=(0.85, 0.52))
    plt.show()


#%% Power density evolution plots 

t1 = 454
t2 = 456
t3 = 458
t4 = 460
t5 = 461
xmin = 0
xrange = 100
#y1 = -10
#y2 = 25
dist = np.arange(0, np.shape(rho_DT)[1], 1)

def powerdens():
    fig = plt.figure(figsize = [6.5, 10])
    gs = fig.add_gridspec(5, 1, hspace = 0.35)
    (ax1), (ax2), (ax3), (ax4), (ax5) = gs.subplots()
    ax1.plot(dist, alpha[t1]/1e29, color='black')
    ax1.plot(dist, rad[t1]*scalefactor/1e29, color='cornflowerblue')
    ax1.plot(dist, elec_cond[t1]/1e29, color='sienna')
    ax1.plot(dist, pdv[t1]/1e29, color='indigo')
    ax1.set_xlim(xmin, xrange)
    #ax1.set_ylim(y1, y2)
    ax1.grid()
    ax1.set_title('t = {} ns'.format(t1*ps_per_index/1e3), x=0.5)
    #ax1.set_title('t = $t_0$', x=0.5)
    ax1.axes.xaxis.set_ticklabels([])
    ax2.plot(dist, alpha[t2]/1e29, color='black')
    ax2.plot(dist, rad[t2]*scalefactor/1e29, color='cornflowerblue')
    ax2.plot(dist, elec_cond[t2]/1e29, color='sienna')
    ax2.plot(dist, pdv[t2]/1e29, color='indigo')
    ax2.set_xlim(xmin, xrange)
    #ax2.set_ylim(y1, y2)
    ax2.grid()
    ax2.set_title('t = {} ns'.format(t2*ps_per_index/1e3, x=0.5))
    #ax2.set_title('t = $t_0$ + 40ps', x=0.5)
    ax2.axes.xaxis.set_ticklabels([])
    ax3.set_ylabel('Power Density x $10^{29}$ W/m\N{SUPERSCRIPT THREE}')
    ax3.plot(dist, alpha[t3]/1e29, color='black')
    ax3.plot(dist, rad[t3]*scalefactor/1e29, color='cornflowerblue')
    ax3.plot(dist, elec_cond[t3]/1e29, color='sienna')
    ax3.plot(dist, pdv[t3]/1e29, color='indigo')
    ax3.set_xlim(xmin, xrange)
    #ax3.set_ylim(y1, y2)    
    ax3.grid()
    ax3.set_title('t = {} ns'.format(t3*ps_per_index/1e3), x=0.5)
    #ax3.set_title('t = $t_0$ + 80ps', x=0.5)
    ax3.axes.xaxis.set_ticklabels([])
    ax4.plot(dist, alpha[t4]/1e29, color='black')
    ax4.plot(dist, rad[t4]*scalefactor/1e29, color='cornflowerblue')
    ax4.plot(dist, elec_cond[t4]/1e29, color='sienna')
    ax4.plot(dist, pdv[t4]/1e29, color='indigo')
    ax4.set_xlim(xmin, xrange)
    #ax4.set_ylim(y1, y2)
    ax4.grid()
    ax4.set_title('t = {} ns'.format(t4*ps_per_index/1e3), x=0.5)
    #ax4.set_title('t = $t_0$ + 100ps', x=0.5)
    ax4.axes.xaxis.set_ticklabels([])
    ax5.plot(dist, alpha[t5]/1e29, color='black', label='$W_{α}$')
    ax5.plot(dist, rad[t5]*scalefactor/1e29, color='cornflowerblue', label='$W_{γ}$')
    ax5.plot(dist, elec_cond[t5]/1e29, color='sienna', label='$W_{e}$')
    ax5.plot(dist, pdv[t5]/1e29, color='indigo', label='$W_{PdV}$')
    ax5.set_xlabel('Capsule Radius (μm) ')
    ax5.set_xlim(xmin, xrange)
   # ax5.set_ylim(y1, y2)
    ax5.grid()
    ax5.set_title('t = {} ns'.format(t5*ps_per_index/1e3), x=0.5)
    #ax5.set_title('t = $t_0$ + 120ps', x=0.5)
    fig.legend(loc='lower right', bbox_to_anchor=(1.05, 0.76))
    
#%%