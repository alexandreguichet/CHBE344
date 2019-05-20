import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#calc terminal velocity and subtract velocity of flow up to get elutriation velocity
rho = .10 #kh/m^3, density of water
rho_g = 7.50 #density of ga
rho_s = 2.65 #density of silica

dmax_g = 0.0725 #in cm
dmax_s = 0.0325
dmin_g = 0.0055
dmin_s = 0.00075

Sph = 0.806 #sphericity, the trident looking thing

darray = np.linspace(dmin_g, dmax_g, 10000)

mu = 14
y0 = [0,1e-15]
tMax = 10
nTime = 1000
g = 9.81

#for sphericity 0.806:
A = 0.2734
B = 0.5510
C = 1.406
D = 762.39

volume = 4/3*np.pi*(dmax_s/2)**3 #calculate volume of largest diameter silica particle
m = (rho_s*volume) #calculate mass of largest diameter silica particle
A = np.pi/4.0*dmax_s**2
uc = np.sqrt(2*m*g/rho*(1-rho/rho_s)/A)
tc = uc*(1-rho/rho_s)/g


volume2 = 4/3*np.pi*(dmin_g/2)**3 #calculate volume of smallest diameter galena particle
m2 = (rho_g*volume2) #calculate mass of smallest diameter galena particle
A2 = np.pi/4.0*dmin_g**2
uc2 = np.sqrt(2*m2*g/rho*(1-rho/rho_g)/A2)
tc2 = uc2*(1-rho/rho_g)/g

def drag_coeff(Re):
    Cd = (24.0/Re)*(1.0+A*(Re**B))+C/(1.0+D/Re)
    return Cd

def settling(y,t,rho,d,mu,uc):
    z,v = y
    dy_dt = np.zeros_like(y)
    Re = rho*d*v*uc/mu
    DragCoeff = drag_coeff(Re)
    dy_dt[0] = v
    dy_dt[1] = 1- DragCoeff*v**2
    return dy_dt


t = np.linspace(0, tMax, nTime)
y = odeint(settling,y0,t,args=(rho,dmax_s,mu,uc))
z = y[:,0]
v = y[:,1]
t2 = np.linspace(0, tMax, nTime)
y2 = odeint(settling,y0,t,args=(rho,dmin_g,mu,uc2))
z2 = y2[:,0]
v2 = y2[:,1]

ind = np.argwhere(v >= 0.9 * v[-1])
ind2 = np.argwhere(v2 >= 0.9 * v2[-1])

#Question 4 - a/b

tEq = np.zeros(len(darray))

#inp = np.zeros(2)
inp = []

a = 0

for i in range(len(darray)):

    r3 = darray[i]/2
    volume3 = 4/3*np.pi*r3**3 #calculate volume of largest diameter silica particle
    m3 = (rho_g*volume3) #calculate mass of largest diameter silica particle
    A3 = np.pi/4.0*darray[i]**2
    uc3 = np.sqrt(2.0*m3*g/rho*(1-rho/rho_g)/A3)
    tc3 = uc3*(1.0-rho/rho_g)/g
    t3 = np.linspace(0, tMax, nTime)
    y3 = odeint(settling,y0,t3,args=(rho,darray[i],mu,uc3))
    z3 = y3[:,0]
    v3 = y3[:,1]

    ind3 = np.argwhere(v3>=0.9*v3[-1])
    tEq[i] = v3[ind3[0]]*uc3

    if (v[ind[0]]*uc - tEq[i]) > 0.1e-7: #This if loop finds the last diameter smaller than the maximum velocity of the silica particle.
        a += 1 #it means that the diameter which is bigger by 1 in the array won't be carried over.

ind = np.argwhere(v>=0.9*v[-1])

print("v silicat max", v[ind[0]]*uc)

ind2 = np.argwhere(v2>=0.9*v2[-1])

print("v galena min", v2[ind2[0]]*uc2)

#Question 4 - c

print(a)
r4 = darray[a]/2
volume4 = 4/3*np.pi*r4**3
m4 = (rho_g*volume4)
A4 = np.pi/4.0*darray[a]**2
uc4 = np.sqrt(2.0*m4*g/rho*(1-rho/rho_g)/A4)
tc4 = uc4*(1.0-rho/rho_g)/g
t4 = np.linspace(0, tMax, nTime)
y4 = odeint(settling,y0,t4,args=(rho,darray[a],mu,uc4))
z4 = y4[:,0]
v4 = y4[:,1]

ind4 = np.argwhere(v4>=0.9*v4[-1])
print("v galena at a speed of", v4[ind4[0]]*uc4, "At diameter ", darray[a])

