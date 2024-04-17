#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 12:44:53 2024

@author: ryansmithers
"""

# Converters
lb2kg = 0.4536
in2m = 0.0254


# Clear all variables
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
plt.close('all')


def f_to_Kelvin(f):
    return ((f-32) * (5/9)) + 273

def Kelvin_to_f(K):
    return (9/5) * (K - 273.15) + 32

def R_Cond_Cylinder (r2, r1, k, L):
    return math.log(r2/r1) / (2*math.pi*k*L)

def R_Conv_Cylinder (d, h, L):
    return 1/(h*math.pi*d*L)

def Temperature_Iteration(Q, R_air, R_prop, T_air, T_prop):
    # Thermal model at node n
    T_s = ((R_prop * R_air * Q) + (R_prop * T_air) + (R_air * T_prop)) / (R_air + R_prop) #K
    
    
    # q_air = (T_s - T_air) / R_air
    q_prop = (T_s - T_prop) / R_prop

    # Change in temperature
    T_prop2 = (((q_prop*t) / (rho*V_inner*cp)) * T_prop) + T_prop
    return T_prop2


# Solver numerical parameters
time_step = 2 #s
max_iterations = 1000


# Transport Properties
h_air = 5 #W/m^2-K
h_prop = 33 #W/m^2-K

# Tank Dimensions
d_outer_imperial = 4 #in
L_imperial = 16.75 #in
m_bot_imperial = 9 #lb
cp_bottle = 490 #J/kg-C
d_outer = d_outer_imperial * in2m #m
L = L_imperial * in2m #m
m_bottle = m_bot_imperial * lb2kg #kg

# Environmental Conditions
T_set_imperial = 55 #F
T_air_imperial = T_set_imperial #f
T_prop_initial_imperial = T_set_imperial #F
T_air = f_to_Kelvin(T_air_imperial) #K
T_prop_initial = f_to_Kelvin(T_prop_initial_imperial) #K
T_bottle_initial = T_prop_initial #K... keep it simple

# Heat Input
Q = 240 #W

# Propellant State Properties
cp_prop = 1800 #J/kg-K 
m_prop_imperial = 2 #lb
m_prop = m_prop_imperial * lb2kg #kg


# Define resistances
R_air = R_Conv_Cylinder(d_outer, h_air, L) #K/W
R_prop = R_Conv_Cylinder(d_outer, h_prop, L) #K/W.... using outer diameter


# Solving end state steady state heat transfer
T_prop_final = Q * R_air + T_air

# Numerical Thresholds
error = 100
threshold = 1 # Percentage error
n = 0

# Arrays: Initialization
T_prop_history = np.zeros(max_iterations)
T_bottle_history = np.zeros(max_iterations)
q_inf = np.zeros(max_iterations)
q_prop = np.zeros(max_iterations) 
q_bottle = np.zeros(max_iterations) 

# Arrays: Initial Conditions
T_prop_history[0] = T_prop_initial
T_bottle_history[0] = T_bottle_initial

while error > threshold and n+1 < max_iterations:
    
    #Solve for heat flow values
    q_inf[n] = (T_bottle_history[n] - T_air) / R_air
    q_prop[n] = (T_bottle_history[n] - T_prop_history[n]) / R_prop
    q_bottle[n] = Q - q_inf[n] - q_prop[n]
    
    #Solve for temperature change
    dT_prop = (time_step * q_prop[n]) / (m_prop * cp_prop)
    dT_bottle = (time_step * q_bottle[n]) / (m_bottle *cp_bottle)
    
    # New Temperature
    T_prop_history[n+1] = T_prop_history[n] + dT_prop
    T_bottle_history[n+1] = T_bottle_history[n] + dT_bottle
    
    # Check Error: Ends on n+1 if condition met
    error = ((T_prop_final - T_prop_history[n+1]) / T_prop_final) * 100
    
    # Iterate Count
    n += 1


# Time array
time = np.arange(0, time_step * n+1, time_step)

# Remove the trailing zeros in all the arrays
T_prop_history = T_prop_history[:n+1]
T_bottle_history = T_bottle_history[:n+1]
q_inf = q_inf[:n+1]
q_prop = q_prop[:n+1]
q_bottle = q_bottle[:n+1]

# Add another set of values at end of heat flow so arrays match (could also remove one off time and temp)
q_inf[n] = q_inf[n-1]
q_prop[n] = q_prop[n-1]
q_bottle[n] = q_bottle[n-1]

# In the units we want
time_minute = time/60
T_prop_history_Fahrenheit = Kelvin_to_f(T_prop_history)
T_prop_final_Fahrenheit = Kelvin_to_f(T_prop_final)
T_bottle_history_Fahrenheit = Kelvin_to_f(T_bottle_history)

# Generate 300 points from min(x) to max(x)
x_smooth = np.linspace(time_minute.min(), time_minute.max(), 300)

# Create spline of order 3 (cubic)
y_smooth_prop = make_interp_spline(time_minute, T_prop_history_Fahrenheit)(x_smooth)
y_smooth_bottle = make_interp_spline(time_minute, T_bottle_history_Fahrenheit)(x_smooth)
y_smooth_q_inf = make_interp_spline(time_minute, q_inf)(x_smooth)
y_smooth_q_prop = make_interp_spline(time_minute, q_prop)(x_smooth)
y_smooth_q_bottle = make_interp_spline(time_minute, q_bottle)(x_smooth)


plt.figure(dpi=1000)

# Plot Temperature Curves
plt.plot(x_smooth, y_smooth_bottle, color='red', label="Bottle Temperature")  # color is optional
plt.plot(x_smooth, y_smooth_prop, color='blue', label="Prop Temperature")  # color is optional
#plt.hlines(T_prop_final_Fahrenheit, time_minute.min(), time_minute.max(), linestyles="--", colors="orange", label = "Steady-State Temperature")
plt.xlabel("Time (min)")
plt.ylabel("Tempearture (F)")
plt.title("Bottle Heating Time")
plt.legend(loc="lower right")


# Plot Heat Flow Curves
plt.figure(dpi=1000)
plt.plot(x_smooth, y_smooth_q_bottle, color='red', label="Q into Bottle")  # color is optional
plt.plot(x_smooth, y_smooth_q_prop, color='blue', label="Q into Prop")  # color is optional
plt.plot(x_smooth, y_smooth_q_inf, color='orange', label="Q into Air")  # color is optional
#plt.hlines(T_prop_final_Fahrenheit, time_minute.min(), time_minute.max(), linestyles="--", colors="orange", label = "Steady-State Temperature")
plt.xlabel("Time (min)")
plt.ylabel("Heat Flow Rate (W)")
plt.title("Where Heat Flows")
plt.legend(loc="upper right")

