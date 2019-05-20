import csv
import numpy as np
import matplotlib.pyplot as plt
import statistics


Sieve_Aperture_1 = [] #[microns]
Feed_1 = [] #[mass fraction]
Accepts_1 = [] #[mass fraction]
Rejects_1 = [] #[mass fraction]

with open('HW2-1 Data.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        Sieve_Aperture_1.append(row[0])
        Feed_1.append(row[1])
        Accepts_1.append(row[2])
        Rejects_1.append(row[3])

Sieve_Aperture = np.array(Sieve_Aperture_1, dtype=float)
Feed = np.array(Feed_1, dtype=float) #It is xf
Accepts = np.array(Accepts_1, dtype=float) #It is xa
Rejects = np.array(Rejects_1, dtype=float) #It is xr

E_a = []
E_r = []
E_f = []

IndustrialAperture = (1 - 10*0.0401)/10*25400

counter = 0 #Counter that will help to sum all the values bellow the Sieve Diameter of 1397

for i in range(0, len(Feed)):
    if Sieve_Aperture[i] <= 1387:
        E_a.append(Accepts[i])
        E_r.append(Rejects[i])
        E_f.append(Feed[i])
        counter += 1

Xa = np.sum(E_a) + (IndustrialAperture - 1397)/(1651-1397) * Accepts[counter]
Xr = np.sum(E_r) + (IndustrialAperture - 1397)/(1651-1397) * Rejects[counter]
Xf = np.sum(E_f) + (IndustrialAperture - 1397)/(1651-1397) * Feed[counter]

print(Xa,Xr,Xf)

OverallEfficiency = (Xa/Xf)*((Xf-Xr)/(Xa-Xr))*((Xa-Xf)/(Xa-Xr))*((1-Xr)/(1-Xf))
print(OverallEfficiency)

