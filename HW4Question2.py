import csv
import numpy as np
from matplotlib.ticker import FormatStrFormatter

import matplotlib.pyplot as plt
from scipy.optimize import fmin

bulkP11 = []
adsorbP11 = []
bulkP22 = []
adsorbP22 = []
N = 30 # number of sample
gGuess = 20.0
KeqGuess = 4.0

with open("HW4 Data.csv", "r") as csv_file:
    # csv_reader = csv_file.readlines()
    # csv_file.close()

    csv_reader = csv.DictReader(csv_file, delimiter=',')
    line_count = 0

    for row in csv_reader:
        bulkP11.append(row['Bulk Concentration Protein1 [mg/mL]'])
        adsorbP11.append(row['Adsorbed Concentration Protein1 [mg/mL]'])
        bulkP22.append(row['Bulk Concentration Protein2 [mg/mL]'])
        adsorbP22.append(row['Adsorbed Concentration Protein2 [mg/mL]'])

        #print(bulkP1, '\n\n', adsorpP1, '\n\n', bulkP2, '\n\n', adsorpP2, '\n')

bulkP1 = []
bulkP1 = np.array(bulkP11, dtype = float)

adsorbP1 = []
adsorbP1 = np.array(adsorbP11, dtype = float)

bulkP2 = []
bulkP2 = np.array(bulkP22, dtype = float)

adsorbP2 = []
adsorbP2 = np.array(adsorbP22, dtype = float)

def langmuir(C, G, Keq):
    return G * C / (1.0 / Keq + C)


def fr(C, m, n):
    return m * C ** (1.0 / n)


def linear(C, K):
    return C * K


def obj_funF(p, cData, qData):
    print(p, " here")
    m = p[0]
    n = p[1]
    qtheory = fr(cData, m, n)
    error = np.sum((qtheory - qData) ** 2)
    return error


def obj_funL(p, cData, qData):
    G = p[0]
    Keq = p[1]
    qtheory = langmuir(cData, G, Keq)
    error = np.sum((qtheory - qData) ** 2)
    return error


def obj_funLin(p, cData, qData):
    K = p
    qtheory = linear(cData, K)
    error = np.sum((np.array(qtheory) - np.array(qData)) ** 2)
    return error


cData = bulkP1
qData = adsorbP1
cData2 = bulkP2
qData2 = adsorbP2
#Plot Guess for particle 1

guess = [gGuess, KeqGuess]
sol = fmin(obj_funL, guess, args=(cData, qData))
GFit = sol[0]
KeqFit = sol[1]

sol = fmin(obj_funF, guess, args=(cData, qData))
mFit = sol[0]
nFit = sol[1]

sol = fmin(obj_funLin, 6.0, args=(cData, qData))
kFit = sol

#Plot Guess for particle 2

guess = [gGuess, KeqGuess]
sol = fmin(obj_funL, guess, args=(cData2, qData2))
GFit2 = sol[0]
KeqFit2 = sol[1]

sol = fmin(obj_funF, guess, args=(cData2, qData2))
mFit2 = sol[0]
nFit2 = sol[1]

sol = fmin(obj_funLin, 6.0, args=(cData2, qData2))
kFit2 = sol

# #Plot Guess
# qGuess = langmuir(cData, guess[0], KeqGuess)
# plt.plot(cData, qGuess, 'r')
#print(kFit)

#Plot Best Fit for particle 1
plt.figure(1)
plt.plot(bulkP1,adsorbP1,'k-',label = "Actual Data",lw=2)
plt.xlabel(r'Bulk Concentration [mol/m$^3$]', fontsize=14)
plt.ylabel(r'Adsorbed Concentration [mol/m$^3$]', fontsize=14)
qFit = langmuir(cData,GFit, KeqFit)
plt.plot(cData, qFit, 'y:',label="Langmuir Isotherm",lw=2)
qFit2 = fr(cData, mFit, nFit)
plt.plot(cData,qFit2, 'c.',label="Freundlich Isotherm",lw=2)
qFit3 = linear(cData, kFit)
plt.plot(cData, qFit3, 'm--',label="Linear Isotherm", lw=2)
plt.title("Best Fit parameters for particle 1 with the Langmuir, Freundlich, and Linea Isotherms")
plt.legend()

#Plot Best Fit for particle 2
plt.figure(2)
plt.plot(bulkP2,adsorbP2,'k-',label = "Actual Data",lw=2)
plt.xlabel(r'Bulk Concentration [mol/m$^3$]', fontsize=14)
plt.ylabel(r'Adsorbed Concentration [mol/m$^3$]', fontsize=14)
qFit2 = langmuir(cData2,GFit2, KeqFit2)
plt.plot(cData2, qFit2, 'y:',label="Langmuir Isotherm",lw=2)
qFit22 = fr(cData2, mFit2, nFit2)
plt.plot(cData2,qFit22, 'c.',label="Freundlich Isotherm",lw=2)
qFit32 = linear(cData2, kFit2)
plt.plot(cData2, qFit32, 'm--',label="Linear Isotherm", lw=2)
plt.title("Best Fit parameters for particle 2 with the Langmuir, Freundlich, and Linea Isotherms")
plt.legend()
plt.show()
