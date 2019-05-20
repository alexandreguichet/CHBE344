import numpy as np
import matplotlib.pyplot as plt
t = 388.15
pabs = 1.15
g = 9.81
atm = 101325
rho = 786
floss1 = 7500
NPSHa = 2.5
mu = 0.278
pstar = 10**(4.07827-1343.943/(388.15-53.773))*10**5
hvap = pstar/(rho*g)
hatm = pabs*atm/(rho*g)
hf = floss1/(rho*g)
ht = hvap + NPSHa + hf - hatm
htt = ht + mu
print("pstar",pstar)
print("hvap",hvap)
print("hatm",hatm)
print("hf",hf)
print("ht is answer to part a",ht)
print("htt",htt)
x = np.linspace(0.5,1.5,9)
y = [htt,htt,htt,htt,htt,htt,htt,htt,htt]
print("y",y)
plt.plot(x,y)
plt.xlabel("Pressure(atm)")
plt.ylabel("Minimum Height(m)")
plt.title("Pressure vs Minimum Height")
plt.show()