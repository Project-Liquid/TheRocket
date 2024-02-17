#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 23:50:58 2024

@author: ryansmithers
"""

#Important Libraries
import math
import matplotlib.pyplot as plt
import numpy as np
import CoolProp
from CoolProp.CoolProp import PropsSI



def Iterator(m, area, K, prop, R, h, nodes):
    n = 0 #Start the iteration at count 0
    
    while n<nodes-1:
        #Solve for velocity and mach number at node
        VelocitySolver(n, m, area, prop, R)
    
        #Solve for new Pressure
        Darcy(K, n)
        
        #Solve for other necessary parameters
        rho_ox_SI[n+1] = PropsSI('D', 'P', P_ox_SI[n], 'H', h, prop)
        x_ox[n+1] = PropsSI('Q', 'P', P_ox_SI[n], 'H', h, prop)
        T_ox_SI[n+1] = PropsSI('T', 'P', P_ox_SI[n], 'H', h, prop)
        
        n+=1
    
    #Filling out final element
    v_ox_SI[n] = v_ox_SI[n-1]
    M_ox[n] = M_ox[n-1]


def VelocitySolver (n, m, area, prop, R):
    v_ox_SI[n] = m/(rho_ox_SI[0]*area) #Velocity at this node
    Cp = PropsSI('CPMASS', 'P', P_ox_SI[n], 'Q', x_ox[n], prop)
    Cv = PropsSI('CVMASS', 'P', P_ox_SI[n], 'Q', x_ox[n], prop)
    gamma = Cp/Cv
    c = np.sqrt(gamma * R * T_ox_SI[n]) #Speed of sound at this node
    M_ox[n] = v_ox_SI[n]/c #Inputs mach number at that node


def Darcy (K, n):
    dP = K * 0.5 * rho_ox_SI[n] * v_ox_SI[n]**2
    P_ox_SI[n+1] = P_ox_SI[n] - dP

#Numerical Solver Tweaks
nodes = 100 #How many times we discretize the tube

#Unit Conversions
psi_to_Pa = 6894.76
m2_to_in2 = 1550 
kg_to_lb = 2.205
lb_per_ft3_to_kg_per_m3 = 16.0185
m_to_ft = 3.28

#Starting Propellant Characteristics: Ox
P1_ox_init = 750 #PSI... assuming thermodynamic equilibrium
x_ox_init = 0 #Quality
P1_ox_init_SI = P1_ox_init * psi_to_Pa


#Starting Propellant Characteristics: Fuel
P1_fuel_init = 400 #PSI... assuming thermodynamic equilibrium
x_fuel = 0
P1_fuel_SI = P1_fuel_init * psi_to_Pa

#Assumptions
m_ox = 3.602 #lb/s
m_ox_SI = m_ox / kg_to_lb
m_fuel = 0.721 #lb/s

#Line Characteristics: Ox
L_ox = 10 #ft
OD_ox = 0.75 #in
Thick_ox = 0.0375 #in
ID_ox = OD_ox - 2*Thick_ox
Area_ox = (math.pi * ID_ox**2)/4
Area_ox_SI = Area_ox / m2_to_in2

#Line Characteristics: Fuel
L_fuel = 1 #ft
OD_fuel = 0.75 #in
Thick_fuel = 0.0375 #in
ID_fuel = OD_ox - 2*Thick_fuel

#Line Characteristics: General
f = 0.02

#Discretized
l_ox = L_ox / nodes #ft
l_fuel = L_fuel / nodes #ft

#Descretized K factor for each line
K_ox = (f * l_ox) / (ID_ox/12)
K_fuel = (f * l_fuel) / (ID_fuel/12)

#Propellant Characteristics
nitrous = 'NITROUSOXIDE'
MOL_nitrous = 0.044013 #kg/mol
R_nitrous = 8.314/MOL_nitrous #J/kg-K

rho_ox_SI_init = PropsSI('D', 'P', P1_ox_init_SI, 'Q', x_ox_init, nitrous)
T_ox_SI_init = PropsSI('T', 'P', P1_ox_init_SI, 'Q', x_ox_init, nitrous)
h_ox_SI_init = PropsSI('H', 'P', P1_ox_init_SI, 'Q', x_ox_init, nitrous)

#All the lists we are keeping track of
P_ox_SI = np.empty(nodes)
x_ox = np.empty(nodes)
T_ox_SI = np.empty(nodes)
rho_ox_SI = np.empty(nodes)
v_ox_SI = np.empty(nodes)
M_ox = np.empty(nodes)


P_ox_SI[0] = P1_ox_init_SI
x_ox[0] = x_ox_init
T_ox_SI[0] = T_ox_SI_init
rho_ox_SI[0] = rho_ox_SI_init


dist_ox = np.linspace(0,L_ox,nodes)


# Solver
Iterator(m_ox_SI, Area_ox_SI, K_ox, nitrous, R_nitrous, h_ox_SI_init, nodes)

#Imperial Units
P_ox = P_ox_SI / psi_to_Pa
T_ox = (T_ox_SI - 273.15) * (9./5) + 32
rho_ox = rho_ox_SI / lb_per_ft3_to_kg_per_m3
v_ox = v_ox_SI * m_to_ft


# Graphing
# Sample data
x = dist_ox
y1 = P_ox
y2 = T_ox
y3 = rho_ox
y4 = x_ox
y5 = v_ox
y6 = M_ox


# Create a figure with two subplots
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, dpi = 1000)

# Plot the first graph in the first subplot
ax1.plot(x, y1)
ax1.set_title('Pressure')

# Plot the second graph in the second subplot
ax2.plot(x, y2)
ax2.set_title('Temperature')

# Plot the second graph in the second subplot
ax3.plot(x, y3)
ax3.set_title('Density')

# Plot the second graph in the second subplot
ax4.plot(x, y4)
ax4.set_title('Quality')

# Plot the second graph in the second subplot
ax5.plot(x, y5)
ax5.set_title('Velocity')

# Plot the second graph in the second subplot
ax6.plot(x, y6)
ax6.set_title('Mach #')

# Show the figure
plt.tight_layout()
plt.show()
    

    


