import matplotlib.pyplot as plt
import numpy as np

ParticleInitialSize1 = 60 # mm
ParticleFinalSize1 = 10 # mm

EnergyComsuption = 13.7 # kW/(Kg/s)

ParticleInitialSize2 = 75 # mm
ParticleFinalSize2 = 5
FinalRange = [50,5]


CKicks = EnergyComsuption/np.log(ParticleInitialSize1/ParticleFinalSize1)
CRittinger = EnergyComsuption/(1/ParticleFinalSize1 - 1/ParticleInitialSize1)

print(CKicks)
print(CRittinger)

Range = list(range(50,4,-1))

#Rittinger

EnergyRittinger = CRittinger*(1/np.array(Range, dtype = float) - 1/ParticleInitialSize2)

Energyfor5mmRittinger = CRittinger*(1/ParticleFinalSize2-1/ParticleInitialSize2)

print("The Particle Final Energy for 5mm ", Energyfor5mmRittinger)
plt.plot(Range,EnergyRittinger)
plt.gca().invert_xaxis()
plt.xlabel("Finale size of the Particle from 50mm to 5mm")
plt.ylabel("Energy Consumption in kW*s/kg")
plt.title("Energy Consumption to crush a particle from 75 mm to a range of 50mm to 5mm assuming Rittinger's law ")

#Kicks

EnergyKicks = CKicks*np.log(ParticleInitialSize2/np.array(Range, dtype = float))

Energyfor5mmKicks = CKicks*np.log(ParticleInitialSize2/ParticleFinalSize2)
print("The Particle Final Energy for 5mm ", Energyfor5mmKicks)

plt.figure()

plt.plot(Range,EnergyKicks)
plt.gca().invert_xaxis()
plt.xlabel("Finale size of the Particle from 50mm to 5mm")
plt.ylabel("Energy Consumption in kW*s/kg")
plt.title("Energy Consumption to crush a particle from 75 mm to a range of 50mm to 5mm assuming Kick's law ")
plt.show()