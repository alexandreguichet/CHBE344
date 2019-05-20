import numpy as np
import matplotlib.pyplot as plt

#continuous cycling centrifuge
rho_liquid = 1205.00 #kg/m^3
rho_solid = 1610.00 #kg/m^3
mu_liquid = 1.95e-3 #Pa/s
l = 0.025 #thickness of liquid layer in meters
ri = 0.575/2-l #inner wall diameter in meters
ro = (0.575)/2 #outer wall in meters. Inner wall diameter+wall thickness
d = 0.375 #depth in meters
w = 1200.00*2*np.pi/60 #rotation rate in rpm
pc = 35*10**(-6) #particle cut size in meters
A = np.pi*(pc/2)**2 #area of cylindrical centrifuge in m^2
vp = 4/3*np.pi*(pc/2)**3 #volume of particles
mp = (pc/2)**3 *4/3*np.pi * rho_solid #mass of particles

#Q4a
#calc separation capacity in m^3/h for a particle cut size of 35 micrometers. Assume spherical particles

Q = (np.pi*(ro**2-ri**2)*d*mp*w**2*(1-rho_liquid/rho_solid)*pc)/(12*mu_liquid*A*np.log(2*ro/(ro+ri)))

print(Q*3600)

d_p = np.linspace(5e-6, 1e-3, 1000)
A_p = np.pi*(d_p/2)**2
m_p = (d_p/2)**3 *4/3*np.pi * rho_solid  # kg
Q = 60*(np.pi*(ro**2 - ri**2)* m_p * w**2 * (1-rho_liquid/rho_solid)*d_p)/(12*mu_liquid*A_p*np.log(2*ro/(ro+ri)))

plt.figure(3)
plt.plot(d_p, Q, 'brown')
plt.xlabel('Particle Diameter [m]')
plt.ylabel(r'Separation Capacity [m$^3$/h]')
plt.title('Diameter of particle vs Separation Capacity')
plt.show()