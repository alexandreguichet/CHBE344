import numpy as np
import matplotlib.pyplot as plt

#Question 4


T1 = 365 # K
T2 = 273 # K

a = list(range(365,273,-1))
TemperatureRange = np.array(a, dtype=float)

# part a
H20 = 100 - 31.2
H20T = H20 - ((1/10*(365-TemperatureRange))/100)*H20 #H2O mass left with temperature
BaNO3 = 31.2

BaNO3T = 9.64*10**-4*(TemperatureRange)**2-.3284*TemperatureRange+22.718 #BaNO3 mass left with temperature in the solution


PercentageRemoved = ((BaNO3 - BaNO3T/100*H20T)/BaNO3)*100 #Calculations for the new weight mass with both mass water and mass Bano3
BaNO3290 = 9.64*10**-4*(290)**2-.3284*290+22.718
H20290 = H20 - (1/10*(365-290))*H20*0.01
PercentageRemoved290 = ((31.2 - BaNO3290/100*H20290)/31.2)*100

print(PercentageRemoved290)



plt.plot(TemperatureRange,PercentageRemoved)
plt.gca().invert_xaxis()
plt.title("Percentage of Ba(NO3)2 removed from the initial solution as the temperature decreases from 365K to 273K")
plt.ylabel("Percentage of Ba(NO3)2 removed of the solution")
plt.xlabel("Temperature in K")
plt.show()