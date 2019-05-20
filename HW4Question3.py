import numpy as np
import matplotlib.pyplot as plt
PbH2O = [5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0]
TbH2O = [32.90,36.18,39.03,41.54,43.79,45.83,49.45,52.58,55.34,57.83,60.09]
hvap=[2423.8,2416,2409.2,2403.2,2397.9,2392.9,2384.2,2376.7,2370,2363.9,2358.4]

mw_solute = 39.97
mw_water = 18.015
P = 821 #kPa
U = 1232 #W/m^2K
Tfeed = 22 + 273.15 #Kelvin
Tsteam = 171.5 + 273.15 #Kelvin
F = 5890 #kg/hr
mfrac_sol_feed = 2.01/100
mfrac_water_feed = 1-mfrac_sol_feed
mfrac_sol_liquid = 50/100
m1 = F*mfrac_sol_feed/100

xf_s = m1/mw_solute/(m1/mw_solute + (F-m1)/mw_water) #mol fraction
L = mfrac_sol_feed*F/mfrac_sol_liquid
V = F-L
P_Tank = 6.89 #kPa
P_sat = P_Tank/(1-xf_s)
Tboil2 = np.interp(P_sat,PbH2O,TbH2O)+273.15
Tboil = np.interp(P_sat, PbH2O, TbH2O) + 273.15 #K
latent = np.interp(P_sat, PbH2O, hvap) #latent heat for liquor side
#hl - hf #get Cp from thermo textbook
delta_h = (1/mw_water)*(7.243e1*(Tboil-Tfeed)+1.039e-2*(Tboil**2-Tfeed**2)/2- 1.497e-6/3*(Tboil**3-Tfeed**3))
hv_minus_hf = latent+delta_h

#heat generated within system
Q = V*(hv_minus_hf)+ L*delta_h

#heat of vaporization for steam from steam table
hs = 2043+(821-820)*(2039.6-2043.0)/(840-820)

#%% PART A:
"""
Calc following quantities assuming neg enthalpy of mixing and adiabatic tanks
"""
#i: Mass flow of steam needed for this process
msteam = Q/hs
#ii: surface area of evaporator
Q = Q*1000/3600
area= Q/(U*(Tsteam-Tfeed))
#iii: economy of evaporator
economy = V/msteam
print('Mass flow rate in kg/h  = %0.2f' %(msteam))
print('Surface area of the evaporator in m2  = %0.2f' %(area))
print('Economy of evaporator  = %0.2f' %(economy))
#%% PART B:
'''
Plot mass flow rate of steam needed to adjust concentraiton from 2.01wt% to 
60wt% under same assumptions 
'''
mfrac_sol_liq_plot = np.linspace(mfrac_sol_feed, 0.6, 100)
xf_s_plot = m1/mw_solute/(mw_solute + (F-m1)/mw_water)
L_plot = mfrac_sol_feed*F/(mfrac_sol_liq_plot)
V_plot = F-L_plot
Psat_plot = P_Tank/(1-xf_s_plot)
Tboil_plot = np.interp(Psat_plot, PbH2O, TbH2O) + 273.15
delta_h_plot = 1/(mw_water)*(7.243e1*(Tboil_plot-Tfeed)+1.039e-2*(Tboil_plot**2-Tfeed**2)/2- 1.497e-6/3*(Tboil**3-Tfeed**3))
hv_minus_hf_plot = latent+delta_h_plot
Q_plot = V_plot*(hv_minus_hf_plot)+ L_plot*delta_h_plot
hs = 2043+(821-820)*(2039.6-2043.0)/(840-820)
Q_plot = Q_plot*1000/3600
ms_plot = Q_plot/hs
plt.plot(mfrac_sol_liq_plot, ms_plot, 'k-')
plt.xlabel('MW % of NaOH in Liquor Outlet')
plt.ylabel('Mass flow of steam (kg/hr)')
plt.title('Mass Flow Rate of Steam needed vs. Mass Fraction NaOH')
plt.show()