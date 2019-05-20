import numpy as np
import matplotlib.pyplot as plt

np.seterr(divide='ignore', invalid='ignore')

# Input parameter
V_flux = 2.45*0.001/0.092903/60 #[m/sec]
p_min = 0e-4 #[Pa]
p_max = 472e3 #[Pa]
N = 3611
V_total = 1400*0.00378541 #[m3]
mu = 1.01e-3 #[Pa-s]
Cs = 1.25/0.0283168 #[kg/m3]
H = 1/39.37 #[m]
A = 71.7*0.092903 #[m2]
K = 6.94e-12*0.092903 #[m2]
Q = V_flux*A
p_filter = 1/K*mu*H*Q/A
p_cake = p_max
alpha = 7.95e10*(1+4.15e-3*(p_cake/1000)**0.74)


time2 = np.linspace(0,N,N)
time = np.linspace(0,N,N)

P_newP = time2*mu*Cs*V_flux**2*alpha

Volume = np.linspace(V_flux, 0, 3611)
P = np.zeros(N)

t1 = (p_max - mu * H * Q / K / A) / (alpha * mu * Cs * Q ** 2 / A ** 2)
V1 = t1 * Q

V2 = V_total - V1
t2 = alpha * mu * Cs / (2 * p_max * A ** 2) * (V_total ** 2 - V1 ** 2) + mu * H / (p_max * K * A) * (V_total - V1)

for i in range(N):
    if (i < 476):
        P[i] = P_newP[i]
    elif (i >= 476):
        P[i] = p_max

t_total = (t1+t2)#Total time in minutes

print(t_total)
print(t1)
#Question 2b.


#Constant pressure filtration
#Pfilter =

Kc = alpha*mu*Cs/A**2/P
Q0 = K*A*P/mu/H
Q2 = np.zeros(len(Volume))

#calculate the flow rate as a function of time

Q1 = np.zeros(N)
V3 = np.linspace(0, V_flux, 3611-476)
a = 0

for i in range(N):
    if(i < 476):
        Q1[i] = V_flux
    elif(i >= 476):
        Q1[i] = V_flux - (V3[a]/A)
        a +=1

for i in range(len(Volume)):
    if(i < int(t1)):
        Q2[i] = 1 / (1 / Q0[i] + Kc[i] * Volume[i])
    elif(i >= int(t1)):
        Q2[i] = 1 / (1 / Q0[i] + Kc[i] * Volume[i]) + (Q2[i-1] - 1 / (1 / Q0[i] + Kc[i] * Volume[i]))

print(len(Q2))
print(len(Kc))
print(len(Volume))
#Plot Pressure drop vs. time
plt.figure(1)
plt.loglog(time/60, P*np.ones_like(time)/1000, 'k.')
plt.title("Pressure Drop as function of time")
plt.xlabel("Time [min]")
plt.ylabel("Pressure Drop [Pa]")

#Plot Volume flux of filtrate
plt.figure(2)
plt.loglog(time/60,Q1,'r.')
plt.title("Volume Flux as a function of time")
plt.xlabel("Time [min]")
plt.ylabel("Volume Flux [m/s]")

#Plot cake mass vs time
plt.figure(3)
plt.plot(time/60, Q2*time*Cs, 'b.')
plt.yscale('log')
plt.xscale('log')
plt.title("Cake Mass as a function of time")
plt.xlabel("Time [min]")
plt.ylabel("Cake Mass [kg]")
plt.show()

