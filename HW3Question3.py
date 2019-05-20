import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

'''
Question 3
'''
'''
Part A
'''


def up_t(m, w, r, rho_f, rho_p, d, mu, A):

    u_p = (m * w**2 * r * d * (1-(rho_f/rho_p)))/(12*mu*A)
    return u_p


d_p = 157e-6  # particle diameter
rho_p = 2800
rho = 998  # density of water
mu = 1e-3
d = 0.25
r = d/2
g = 50
rpm = np.sqrt(g/(1.118e-3*r))

print('Rotational Speed: ', rpm, 'rpm\n')

m = rho_p * 4/3 * np.pi * (d_p/2)**3
w = rpm*2*np.pi/60  # angular velocity in rads^-1
A = np.pi*(d_p/2)**2

v = up_t(m, w, r, rho, rho_p, d_p, mu, A)

print('Terminal Particle Velocity = ', v, 'm/s')


'''
Part B
'''

r = 0.2  # radius
rpm = np.linspace(60, 4000, 1000)
v = np.empty_like(rpm)

i = 0
for n in rpm:
    w = n * 2 * np.pi / 60  # angular velocity in rads^-1
    v[i] = up_t(m, w, r, rho, rho_p, d_p, mu, A)
    i += 1

plt.figure(1)
plt.plot(rpm, v, '-')
plt.xlabel('RPM')
plt.ylabel('Terminal Settling Velocity (m/s)')
plt.title('Effect of RPM on Settling Velocity at r = 20cm')
plt.show()