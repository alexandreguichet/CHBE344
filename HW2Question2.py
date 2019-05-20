import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve

#Question 2 - a - i

d = 105*10**-6 #m
rho_p = 2650.0 #kg/m^3
volume = (4.0/3.0)*np.pi*(d/2)**3
m = volume*rho_p #kg
rho = 935.0 #kg/m^3
muList = np.linspace(1e-3, 10.0, 1000)
g = 9.81 #kg/s

#m = 10.0
# g = 9.81
#rho = 1000.0
#rho_p = 2036.0
# muList = np.linspace(1e-3, 14.0, 100)
y0 = [0,1e-15]
tMax = 10
nTime = 1000

#volume = m/rho_p
#d = (6.0*volume/np.pi)**(1.0/3.0)
A = np.pi/4.0*d**2
uc = np.sqrt((2*m*g*(1-rho/rho_p)/rho/A))
tc = uc/(1-rho/rho_p)/g

def drag_coeff(Re):

    def fun1(Re):
        return 24.0/Re
    def fun2(Re):
        return 24.0/Re*(1+0.14*Re**(0.7))
    def fun3(Re):
        return 0.445

    c = np.piecewise(Re,[(Re<0.1),(Re>=0.1) & (Re<1e3),(Re>=1e3)],[fun1,fun2,fun3])
    return c


def settling(y,t,rho,d,mu,uc):

    z, v = y
    dy_dt = np.zeros_like(y)
    Re = rho*d*v*uc/mu
    C = drag_coeff(Re)
    print(C)
    dy_dt[0] = v
    dy_dt[1] = 1 - C*v**2
    return dy_dt


nMu = len(muList)
tEq = np.zeros(nMu)

# for i in range(nMu):
#
#     mu = muList[i]
#     t = np.linspace(0, tMax, nTime)
#     y = odeint(settling, y0, t, args=(rho, d, mu, uc))
#     z = y[:,0]
#     v = y[:,1]
#
#     tEq[i] = v[-1]*uc
#
# tEq[-1]=tEq[-2]

#Question 2 - a - ii

newMu = 44.5e-3 #Pa.s

newT = np.linspace(0, tMax, nTime)
newY = odeint(settling, y0, newT, args=(rho, d, newMu, uc))
newZ = newY[:, 0]
newV = newY[:, 1]

ind = np.argwhere(newV >= 0.9 * newV[-1])

v2 = newV[ind[0]]*uc
z2 = newZ[ind[0]]*uc*tc

print("The settling velocity for the particle in vegetable oil (44.5Cp) is:", v2)
print("t max =", newT[ind[0]]*tc)
print("z max = ", newZ[ind[0]]*tc*uc)


plt.loglog(muList, tEq)
plt.xlabel('Viscosity [Pa/s]')
plt.ylabel('Terminal Velocity [m/s]')
plt.title('Terminal velocity as a function of viscosity on a logarithm scale')
plt.show()

