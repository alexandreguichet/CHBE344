import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

print('Problem 1 a)')
m_dot_guess = 1  # [kg/0]
m_dot = m_dot_guess
precision1 = 3


def friction_factor(Re, relative_roughness, D):
    def fun1(Re):
        return 16 / Re

    def fun2(Re):
        return (-3.6 * np.log10(6.9 / Re + (relative_roughness / (3.7 * D)) ** 1.11)) ** (-2)

    f = np.piecewise(Re, [(Re < 4000), (Re >= 4000) & (Re < 10**7)], [fun1, fun2])
    return f

Efficiency = 63.5 / 100  # efficiency of the pump

def Energy_Balance(m_dot):
    Power_cons = 6250  # [W] #power consumption
    z_resevoir = 5  # [m]Reservoir is 5m up from the pump
    z_exchanger = 2.75  # [m]exchanger is 2.75m above pump
    z_pump = 0  # [m]
    P_e = 325000  # [Pa] Gauge pressure required at entrance of exchanger
    P_r = 0  # gauge pressure of resevoir is zero
    D_rp = 4 * 0.0254  # [m] 4 inches converted to meters
    A_rp = (D_rp / 2) ** 2 * np.pi
    L_rp = 125  # [m]
    D_pe = 1 * 0.0254  # [m] 1 inch to meters
    A_pe = (D_pe / 2) ** 2 * np.pi
    L_pe = 6.5  # [m]
    rho = 1000  # [kg/m^3]
    g = 9.81  # [m/s^2]
    relative_roughness = 0.001
    mu = 1.00e-3  # [Pa*s]
    Power_req = Power_cons * Efficiency  # The mechanical power that the pump actually outputs
    u_1 = m_dot / (A_rp * rho)
    u_2 = m_dot / (A_pe * rho)
    Re1 = rho * u_1 * D_rp / mu
    Re2 = rho * u_2 * D_pe / mu
    n_elbow1 = 7  # there are 7 threaded elbows
    Cf_elbow1 = 1.5
    Cf_elbow2threaded = 0.7
    Cf_elbow2flanged = 0.2

    E_f_minor1 = u_1 ** 2 / 2 * (n_elbow1 * Cf_elbow1)
    E_f_minor2 = u_2 ** 2 / 2 * (Cf_elbow2threaded + Cf_elbow2flanged)
    E_f_minor = E_f_minor1 + E_f_minor2

    E_f_major1 = friction_factor(Re1, relative_roughness * D_rp, D_rp) * 2 * L_rp * u_1 ** 2 / D_rp
    E_f_major2 = friction_factor(Re2, relative_roughness * D_pe, D_pe) * 2 * L_pe * u_2 ** 2 / D_pe
    E_f_major = E_f_major1 + E_f_major2

    E_f = E_f_minor + E_f_major
    fun = Power_req + m_dot * (
                1 / 2 * (u_1 ** 2 - (u_2 ** 2)) + g * (z_resevoir - z_exchanger) + 1 / rho * (P_r - P_e) - E_f)
    return fun


m_dot = fsolve(Energy_Balance, m_dot)

rho = 1000  # [kg/m^3]
mu = 1.00e-3  # [Pa*s]
D_pe = 1 * 0.0254  # [m] 1 inch to meters
A_pe = (D_pe / 2) ** 2 * np.pi
u_2 = m_dot / (A_pe * rho)
Re2 = rho * u_2 * D_pe / mu

print('Mass flow rate = ', m_dot, 'kg/s')
print('Reynolds number for water entering the exchanger = ', Re2)


#Part b


#Hydrolic power
HyP = 6250/(9.81*m_dot)
brakePower = HyP/Efficiency



