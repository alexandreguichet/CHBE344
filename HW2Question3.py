import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
import csv
from scipy.integrate import odeint

g = 9.81 #m/s**2
rho = 1594
rho_p = 4.0*1000.0
mu = 0.913e-3

Area = 9.75
q = 20.2/60*0.001
uTerminalMin = q/Area
phi = 0.198

N = 100

dmin = 1e-6
dmax = 5e-3
d = np.linspace(dmin,dmax,N)

vEq = np.zeros(N)
tEq = np.zeros(N)
x = np.zeros(N)

volume = (4 / 3) * np.pi * (d / 2) ** 3
m = volume * rho_p  # kg
A = np.pi / 4.0 * d ** 2
uc = np.sqrt((2 * m * g * (1 - rho / rho_p) / rho / A))
tc = uc / (1 - rho / rho_p) / g

uMax = np.zeros(N)
Revec = np.zeros(N)
n = np.zeros(N)
uHindered = np.zeros(N)
x_removed = np.zeros(N)

y0 = [0.0, 1e-15]
tMax = 10
nTime = 1000

Sieve_Aperture_1 = [] #[microns]

with open('HW2-1 Data.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        Sieve_Aperture_1.append(row[0])

Sieve_Aperture = np.array(Sieve_Aperture_1, dtype=float)
Sieve_Aperture = np.delete(Sieve_Aperture, slice(0,2)).astype(np.float)

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

for i in range(N):

    t = np.linspace(0, tMax, nTime)
    y = odeint(settling, y0, t, args=(rho, d[i], mu, uc[i]))
    z = y[:,0]
    v = y[:,1]

    ind = np.argwhere(v >= 0.9 * v[-1])
    tEq[i] = t[ind[0]] * tc[i]

    u = v*uc[i]
    x = z*uc[i]*tc[i]
    uMax[i] = u[nTime - 1]
    Revec[i] = rho*uMax[i]*d[i]/mu
    n[i] = 3.575 - 1.075*np.tanh(0.6516*(np.log(Revec[i]))-1.153)
    uHindered[i] = uMax[i]*((1-phi)**n[i])
    if uHindered[i]/uTerminalMin >1:
        x_removed[i] = 1
    else:
        x_removed[i] = uHindered[i]/uTerminalMin


davg = np.zeros(len(Sieve_Aperture))
davg[0] = Sieve_Aperture[0]
for i in range(len(Sieve_Aperture)-1):
    davg[i+1] = (Sieve_Aperture[i] + Sieve_Aperture[i+1])/2

davg*=1e-6

print("Diameter [m] of the particle is:", d)
print("Fraction of the particle removed is:", x_removed)
print("Average particle diameter[m] is:", davg)

if davg[0]>d[0]:
    print("The average diameter is bigger than the maximum particle diameter, all particles will settle down")

plt.loglog(d,x_removed)
plt.title("Sphalerite Settling in Sedimentation Unit")
plt.xlabel("Particle Diameter, d(m)")
plt.ylabel("Fraction of Particles Removed, x")
plt.show()