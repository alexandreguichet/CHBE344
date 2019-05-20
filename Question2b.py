import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Question 2 - b - i

sampleSize = 100

d = np.linspace(1e-6, 1e-1, sampleSize) #m
rho_p = 2650.0 #kg/m^3
rho = 827.0 #kg/m^3
muList = 4.91e-3
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

def drag_coeff(Re):

    def fun1(Re):
        return 24.0/Re
    def fun2(Re):
        return 24.0/Re*(1+0.14*Re**(0.7))
    def fun3(Re):
        return 0.445

    c = np.piecewise(Re,[(Re<0.1),(Re>=0.1) & (Re<1e3),(Re>=1e3)],[fun1,fun2,fun3])
    return c


def settling(y,t,rho,da,mu,uc):

    z, v = y
    dy_dt = np.zeros_like(y)
    Re = rho*da*v*uc/mu
    C = drag_coeff(Re)
    dy_dt[0] = v
    dy_dt[1] = 1 - C*v**2
    return dy_dt


nMu = sampleSize
tEq = np.zeros(nMu)

for i in range(nMu):

    volume = (4 / 3) * np.pi * (d[i] / 2) ** 3
    m = volume * rho_p  # kg
    A = np.pi / 4.0 * d[i] ** 2
    uc = np.sqrt((2 * m * g * (1 - rho / rho_p) / rho / A))
    tc = uc / (1 - rho / rho_p) / g

    t = np.linspace(0, tMax, nTime)
    y = odeint(settling, y0, t, args=(rho, d[i], muList, uc))
    z = y[:,0]
    v = y[:,1]

    ind = np.argwhere(v >= 0.9*v[-1])
    tEq[i] = v[-1]*uc

tEq[-1] = tEq[-2]

max(tEq)

#Question 2 - b - ii

newpSize = 1e-3

newVolume = (4 / 3) * np.pi * (newpSize / 2) ** 3
newM = newVolume * rho_p  # kg
newA = np.pi / 4.0 * newpSize ** 2
newuc = np.sqrt((2 * newM * g * (1 - rho / rho_p) / rho / newA))
newtc = newuc / (1 - rho / rho_p) / g

newT = np.linspace(0, tMax, nTime)
newY = odeint(settling, y0, newT, args=(rho, newpSize, muList, newuc))
newZ = newY[:, 0]
newV = newY[:, 1]

newInd = np.argwhere(newV >= 0.9 * newV[-1])

print("The velocity of a 1mm particle is: ", newV[newInd[-1]]*newuc)
print("t max =", newT[ind[-1]]*newtc)
print("z max = ", newZ[ind[-1]]*newtc*newuc)

plt.loglog(d, tEq)
plt.xlabel('Diameter [m]')
plt.ylabel('Terminal Velocity [m/s]')
plt.title('Terminal velocity as a function of diameter in Texas Crude Oil on a logarithm scale')
plt.show()