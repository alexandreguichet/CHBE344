import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve

# Input Parameters
N = 100 # Number of samples

Q = 4.7e-7 # volume flow rate [m^3/s]
A = (0.09/2)**2*np.pi # superficial area [m^2]
L = 1 # depth of porous media [m]
mu = 1.13e-3 # viscosity of fluid [Pa-s]
rho = 998 # kg/m^3
sphericity = 0.730 #[-]

#Part a
# Dp = 0.002
epsilon = 0.498 # porosity [-]

#Part 2
Dp = np.linspace(50e-6,5e-3,N) # Sauter-mean particle diameter [m]

#Part 3

def ergun_k(epsilon, sphericity, Dp, Re):
   kInv = 150*((1-epsilon)**2/(epsilon**3*sphericity**2*Dp**2)*(mu*L*Q/A))
   return kInv


#Calculate the pressure drop based on the input parameters

#Question 1

# Re = rho*Q/A*Dp/mu
# P = ergun_k(epsilon,sphericity,Dp,Re)
#
# print('The pressure Drop is %0.5g' %P)

#Question 2

P = np.zeros(N)

for i in range(N):
    Re = rho*Q/A*Dp[i]/mu
    P[i] = ergun_k(epsilon,sphericity,Dp[i],Re)

plt.loglog(Dp,P)
plt.xlabel("Diameter of Particle [m]")
plt.ylabel("Pressure drop [kPa]")
plt.title("Pressure Drop of sand particle with increasing size from 50e-6m to 5e-3m")
plt.show()

#Question 3
# P = 511.5 # [Pa]
# Dp = 203e-6
# Re = rho * Q / A * Dp / mu
# kInv = P * A / (mu * L * Q)
# sphericity = 1
#
# def solve(epsilon):
#
#    fun = P - 150*((1-epsilon)**2/(epsilon**3*sphericity**2*Dp**2)*(mu*L*Q/A))
#    return fun
#
# epsilon1 = fsolve(solve,0.5)
# print(epsilon1)
