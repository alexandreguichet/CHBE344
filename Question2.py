import csv
import numpy as np
import matplotlib.pyplot as plt
import statistics

# Open csv File

with open('HW1-2 Data.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        line_count += 1
    print(f'Processed {line_count} lines.')

# Variables/Lists

Volume = []
Mass = []
FeretD_1mm = []
FeretD_2mm = []

Density = []
LargestFeret = []
MinimumF = []
MaximumF = []
FeretRatio = []
Dsv = []

bins = 250
binsDensity = 100

with open('HW1-2 Data.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        Volume.append(row[0])
        Mass.append(row[1])
        FeretD_1mm.append(row[2])
        FeretD_2mm.append(row[3])


# Define Calculations for moments, mean, kurtosis, skewness, standard variation

def moments(y, w, y0=0, m=1):
    numerator = 0
    denominator = 0
    N = len(y)
    for i in range(N):
        numerator += w[i] * (y[i] - y0) ** m
        denominator += w[i]

    moment = numerator / denominator
    return moment


def calculations(
        particles):  # Will give an array of (in order): mean, sauter, debrouker, variance, std, skewness, kurtosis

    w = np.ones_like(particles)
    mean = moments(particles, w)
    variance = moments(particles, w, y0=mean, m=2)
    std = np.sqrt(variance)
    skewness = moments(particles, w, y0=mean, m=3) / std ** 3
    kurtosis = moments(particles, w, y0=mean, m=4) / std ** 4


    return mean, std, skewness, kurtosis


# Particle Density and compute mean, kurtosis, etc...

a = np.array(Mass, dtype=np.float)
b = np.array(Volume, dtype=np.float)

Density = a / b

print("length of density", len(Density))

DensityMean, DensityStd, DensitySkewness, DensityKurtosis = calculations(Density)

print("")
print("Mean of Density is", DensityMean)
print("Standard Deviation of Density is", DensityStd)
print("Skewness of Density is", DensitySkewness)
print("kurtosis of Density is", DensityKurtosis)


# Equivalent Spherical particle diameter based on volume

D_e = (6*np.array(Volume, dtype = float)/np.pi)**(1/3)

D_e_toPlot = np.array(D_e, dtype = np.float)

D_e_toPlotMean, D_e_toPlotStd, D_e_toPlotSkewness, D_e_toPlotKurtosis = calculations(D_e_toPlot)

print("")
print("Mean of Equivalent Volume Diameter is", D_e_toPlotMean)
print("Standard Deviation of Equivalent Volume Diameter is", D_e_toPlotStd)
print("Skewness of Equivalent Volume Diameter is", D_e_toPlotSkewness)
print("kurtosis of Equivalent Volume Diameter is", D_e_toPlotKurtosis)

# Largest Feret Diameter
for i in range(0, 1000):
    if FeretD_1mm[i] > FeretD_2mm[i]:
        LargestFeret.append(FeretD_1mm[i])
    else:
        LargestFeret.append(FeretD_2mm[i])

WantedFeret = np.array(LargestFeret, dtype=np.float)

LargestFeretMean, LargestFeretStd, LargestFeretSkewness, LargestFeretKurtosis = calculations(WantedFeret)

print("")
print("Mean of Feret is", LargestFeretMean)
print("Standard Deviation of Feret is", LargestFeretStd)
print("Skewness of Feret is", LargestFeretSkewness)
print("kurtosis of Feret is", LargestFeretKurtosis)

# Particle Aspect Ratio : Ratio of the minimum feret to the maximum feret diameter.

for i in range(0,1000):
    if FeretD_1mm[i] > FeretD_2mm[i]:
        MinimumF.append(FeretD_2mm[i])
        MaximumF.append(FeretD_1mm[i])
    else:
        MinimumF.append(FeretD_1mm[i])
        MaximumF.append(FeretD_2mm[i])

c = np.array(MinimumF, dtype = np.float)
d = np.array(MaximumF, dtype = np.float)

FeretRatio = c/d

FeretRatioMean, FeretRatioStd, FeretRatioSkewness, FeretRatioKurtosis = calculations(FeretRatio)

print("")
print("Mean of Particle Aspect Ratio is", FeretRatioMean)
print("Standard Deviation of Particle Aspect Ratio is", FeretRatioStd)
print("Skewness of Particle Aspect Ratio is", FeretRatioSkewness)
print("kurtosis of Particle Aspect Ratio is", FeretRatioKurtosis)

# Calculate Sauter : dsv = 6V/S

Area_Sauter = 2.0*(c/2.0)**2*np.pi + c*np.pi*d
Dsv = 6*np.array(Volume, dtype = float)/Area_Sauter

DsvMean, DsvStd, DsvSkewness, DsvKurtosis = calculations(Dsv)

print("")
print("Mean of D[3,2] is", DsvMean)
print("Standard Deviation of D[3,2] is", DsvStd)
print("Skewness of D[3,2] is", DsvSkewness)
print("Kurtosis of D[3,2] is", DsvKurtosis)

#Calculate D[4,3]

Ds = []

Ds = np.sqrt(Area_Sauter/np.pi)
Dv = (6*np.array(Volume, dtype = float)/np.pi)**(1/3)

DeBrouckere = Dv**4/Ds**3

DeBrouckereMean, DeBrouckereStd, DeBrouckereSkewness, DeBrouckereKurtosis = calculations(DeBrouckere)

print("")
print("Mean of D[4,3] is", DeBrouckereMean)
print("Standard Deviation of D[4,3] is", DeBrouckereStd)
print("Skewness of D[4,3] is", DeBrouckereSkewness)
print("Kurtosis of D[4,3] is", DeBrouckereKurtosis)

print("")
print("The D50 for smallest feret is ",statistics.median(c))
print("")
print("The D10 and D90 of mass weighted is ", np.percentile(DeBrouckere, 10), " ", np.percentile(DeBrouckere, 90))

plt.hist(Density,bins)
plt.yscale('log')
plt.title("Density plotted with a logarythmic scale")
plt.xlabel("density in g/mm^3")
plt.ylabel("frequency")


plt.figure()

plt.hist(WantedFeret, bins)
plt.title('Largest Feret Diameter')
plt.xlabel('feret diameter in mm')
plt.ylabel("frequency")


plt.figure()

plt.hist(D_e_toPlot, bins)
plt.title("Equivalent Spherical Particle Diameter based on Volume")
plt.xlabel("diameter in mm")
plt.ylabel("frequency")


plt.figure()

plt.hist(FeretRatio, bins)
plt.title("Particle Aspect Ratio")
plt.xlabel("Feret Diameters ratio")
plt.ylabel("frequency")


plt.figure()

plt.hist(Dsv, bins)
plt.yscale('log')
plt.title("Sauter Mean Diameter for particles")
plt.ylabel("frequency")

plt.figure()

plt.hist(DeBrouckere, bins)
plt.yscale('log')
plt.title("De Brouckere Mean Diameter for particles")
plt.ylabel("frequency")

print(FeretRatio)

plt.show()
