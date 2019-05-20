import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def round(x):
    localround = "%#.3g" % (x)
    return localround


H_1 = 5.0  # height of pump below resevoir
H_2 = 2.75  # height of pump below heat exchanger
pGauge = 325000  # required pressure
L_1 = 125.0  # resevoir
D_1 = 4.0 * .0254  # resevoir
L_2 = 6.5  # heat exchanger
D_2 = 1.0 * .0254  # pipe
elbow_num1 = 7.0
roughness = 0.001
A_1 = np.pi / 4 * D_1 ** 2
A_2 = np.pi / 4 * D_2 ** 2

rho = 998
mu = 8.9 * (10 ** -4)
g = 9.81
print('Question 1')
print('a)')
print('\n')
powerCons = 6250
eff = .635
power = powerCons * eff


def solve(mDot, rho, D_1, D_2, A_1, A_2, g, H_1, H_2, pGauge, mu, power, L_1, L_2, elbow_num1):
    U_in = mDot / rho / A_1
    U_out = mDot / rho / A_2
    head_V = (U_in ** 2 - U_out ** 2) / 2.0
    head_Z = g * (H_1 - H_2)
    head_P = -pGauge / rho
    Re1 = rho * U_in * D_1 / mu
    Re2 = rho * U_out * D_2 / mu
    f1 = f(Re1)
    f2 = f(Re2)
    if (A_2 / A_1) <= 0.715:
        E_exp = (2 / 5) * (5 / 4 - A_2 / A_1)
    else:
        E_exp = (3 / 4) * (1 - A_2 / A_1)

    E_elbows = (elbow_num1 * 1.5 * U_in ** 2 / 2) + ((0.7 + 0.2) * U_out ** 2 / 2)
    E_major1 = (f1 * 2.0 * L_1 * U_in ** 2 / D_1)
    E_major2 = (f2 * 2.0 * L_2 * U_out ** 2 / D_2)
    E_contraction = 0.5 * U_in ** 2 / 2
    E_expansion = E_exp * U_out ** 2 / 2

    Ef = E_elbows + E_major1 + E_major2 + E_contraction + E_expansion
    zero = power + mDot * (head_V + head_Z + head_P - Ef)
    return zero


def f(Re):
    if Re > 4000:
        f = (-3.6 * np.log10((6.9 / Re) + (roughness / 3.7) ** 1.11)) ** (-2)
    else:
        f = 16.0 / Re
    return f


mDot = fsolve(solve, 0.5, args=(rho, D_1, D_2, A_1, A_2, g, H_1, H_2, pGauge, mu, power, L_1, L_2, elbow_num1))
U_out = mDot / rho / (np.pi / 4 * D_2 ** 2)
Re2 = rho * U_out * D_2 / mu
print('Mass flow rate is: %s kg/s' % (round(mDot)))
print('Reynolds number is: %s' % (round(Re2)))
print('b)')
print('\n')

n_elbow = np.linspace(0, 50)
L_1_v = np.zeros(len(n_elbow))
L = L_1 - 7 * 3

for i in range(len(n_elbow)):
    L_1_v[i] = L + 3 * i

eff = np.zeros(len(n_elbow))
work_v = np.zeros(len(n_elbow))

U_in = mDot / rho / A_1
U_out = mDot / rho / A_2
head_V = (U_in ** 2 - U_out ** 2) / 2.0
head_Z = g * (H_1 - H_2)
head_P = -pGauge / rho
Re1 = rho * U_in * D_1 / mu
Re2 = rho * U_out * D_2 / mu
f1 = f(Re1)
f2 = f(Re2)

if (A_2 / A_1) <= 0.715:
    E_exp = (2 / 5) * (5 / 4 - A_2 / A_1)
else:
    E_exp = (3 / 4) * (1 - A_2 / A_1)


def pump_eff(eff_pump, w):
    zero = -0.48 * (w / (eff_pump / 100)) ** 2 + 8.2 * (w / (eff_pump / 100)) + 31 - eff_pump
    return zero


for i in range(len(n_elbow)):
    E_elbows = (n_elbow[i] * 1.5 * U_in ** 2 / 2) + ((0.7 + 0.2) * U_out ** 2 / 2)
    E_major1 = (f1 * 2.0 * L_1_v[i] * (U_in ** 2 / D_1))
    E_major2 = (f2 * 2.0 * L_2 * U_out ** 2 / D_2)
    E_contraction = 0.5 * U_in ** 2 / 2
    E_expansion = E_exp * U_out ** 2 / 2

    Ef = E_elbows + E_major1 + E_major2 + E_contraction + E_expansion

    work_v[i] = -mDot * (-U_out ** 2 / 2 + head_Z + head_P - Ef) / 1000
    eff[i] = fsolve(pump_eff, 60, args=(work_v[i]))

P_brake = work_v / (eff / 100)

plt.title('''Brake Horsepower vs. Number of 
Elbows Between Resevoir and Pump''')
plt.xlabel('Number of Regular 90 degree threaded elbows')
plt.ylabel('Brake Horsepower hp')
plt.plot(n_elbow, P_brake * 1.34102)
plt.show()
# c
print('c)')

print('''Brake horsepower only seems to increas by a small amount when additional 
elbows are added and the friction factor increases. Because of this, I would
rather reroute the pipes instead of reducing the pumping energy.''')