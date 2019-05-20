import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#Input variables

Hp = 5.0 #Pump location below liquid level [m]
Hhe = 2.75 #Heat exchanger location above the pump [m]
Phe = 325000 #Required Heat Exchanger inupt [Pa]

LpPR = 125.0 #Pipe length from pump to reservoir [m]
DpPR = 0.1016 #Pipe diameter from pump to reservoir [m]
ApPR = np.pi*(DpPR/2)**2 #Area of pipe from pump to reservoir [m**2]
LpPHe = 1.0 #Pipe length from pump to heat exchanger [m]
DpPhe = 0.0254 #Pipe diameter from pump to heat exchanger [m]
ApPhe = np.pi*(DpPhe/2)**2

Gravity = 9.81
Density = 998

ElbthNumberPR = 7 #Number of threaded elbows that connects the pump to the reservoir [-]
ElbthNumberPHe = 1 #Number of threaded elbows that connects the pump to the heat exchanger [-]
ElbflNumberPhe = 1 #Number of flanged elbows connecting the pump to the heat exchanger [-]

CfthLR = 0.7 #Friction Losses coefficient for Long Radius 90degrees threaded elbows
CfthRr = 1.5 #Friction Losses coefficient for Regular Radius 90Degree threaded elbows
CfflLr = 0.2 #Friction Losses coefficient for Long Radius 90Degree flanged elbows

Rhp = 0.001 #Relative roughtness of all pipes

mu = 1.0e-3

#Part a)

#TotalDynamicHead = Total Discharge Head - Total Suction Head
#Bernouilly : 0 = Total Dynamic Head + Velocity head + Elevation Head + Pressure Head - Friction Head
#Frictional Losses : Ef = Ef,major + Ef,minor
#Ef,minor = Cf*u**2/2 -> Elbow fricion loss * velocity square
#Ef,major = f*(2Lu**2)/d
#f = [-3.6log(6.9/Re+[E/3.7D]**1.11)] for 4000<Re<10**7; 16/Re for Re<2300;


ConsumptionP = 6250.0 #Pump consumption [W]
EfficiencyP = 0.635 #Pump efficiency [-]

def PowerConsumption(Q):
    fun = ConsumptionP + \
           Q/Density*(1/2*((Q/ApPR)**2 - (Q/ApPhe)**2))+ Gravity *(Hp-Hhe) + 1/Density * (0-Phe) - hl(Q)
    return fun

def F1PR(Q):
    return (FanningFrictionFactor(Q,(Rhp*DpPR),DpPR,ApPR)*2*LpPR*(Q/ApPR)**2/DpPR)

def F2PHe(Q):
    return (FanningFrictionFactor(Q,(Rhp*DpPhe),DpPhe,ApPhe)*2*LpPHe*(Q/ApPhe)**2/DpPhe)

def F3(Q):
    return ElbthNumberPR*(CfthRr*(Q/ApPR)**2)/2 + 1*((CfflLr*(Q/ApPhe)**2)/2) + 1*((CfthLR*(Q/ApPhe)**2)/2)

def hl(Q):
    return F1PR(Q) + F2PHe(Q) + F3(Q)


# def FanningFrictionFactor(Re ,E, d):
#     def fun1(Re):
#         return 16.0 / Re
#
#     def fun2(Re):
#         return [-3.6*np.log(6.9/Re+(E/(3.7*d))**1.11)]**-2
#
#     c = np.piecewise(Re, [(Re < 2300), (Re >= 4000) & (Re < 10**7)], [fun1, fun2])
#     return c


def FanningFrictionFactor(Q,E, d, A):
    def fun1(Q):
        return 16.0 / ((Density*Q*d)/(8.9*1**-4*A))

    def fun2(Q):
        return (-3.6*np.log(6.9/((Density*Q*d)/(8.9*1**-4*A))+(E/(3.7*d))**1.11))**-2


    Minimum = 2300/Density/d*mu*np.pi*(d/2)**2
    Medium = 4000/Density/d*mu*np.pi*(d/2)**2
    High = 10e7/Density/d*mu*np.pi*(d/2)**2
    c = np.piecewise(Q, [(Q < Minimum), (Q >= Medium) & (Q < High)], [fun1, fun2])
    return c

# def ExpansionFrictionalLosses(A1,A2):
#     return (1-A1/A2)

# def ContractionFrictionalLosses(A1,A2):
#     def fun1(A1,A2):
#         return 2/5*(5/4-A2/A1)
#
#     def fun2(A1,A2):
#         return 3/4*(1-A2/A1)
#
#     c = np.piecewise(A1,A2, [(A2/A1 < 0.715), (A2.A1 > 0.715)], [fun1, fun2])
#     return c
#
# def TotalDynamicHead(E, Power): #From Efficiency and Power Input
#     return E*Power

MassFlowRate = fsolve(PowerConsumption,1)*Density/6.5

print(MassFlowRate)


