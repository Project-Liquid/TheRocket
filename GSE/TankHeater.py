#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 12:44:53 2024

@author: ryansmithers
"""


# Clear all variables
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
plt.close('all')


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
t = 10 #s
max_iterations = 1000


# Transport Properties
h_air = 25 #W/m^2-K
h_prop = 2.5 #W/m^2-K
k_alum = 237 #K/W-K

# Tank Dimensions
d_outer = 0.1524 #m
d_inner = 0.1397 #m
L = 0.25 #m
V_inner = math.pi * d_inner * L #m^3

# Environmental Conditions
T_air = 266 #K
T_prop_initial = 266 #K (Starting condition of propellant)

# Heat Input
Q = 400 #W

# Propellant State Properties
cp = 1800 #J/kg-K ... come back to this later
rho = 1000 #kg/m^3 ... come back to this later


# Define resistances
R_air = R_Conv_Cylinder(d_outer, h_air, L) #K/W
R_prop = R_Conv_Cylinder(d_inner, h_prop, L) #K/W


# Solving end state steady state heat transfer
T_prop_final = Q * R_air + T_air

# Numerical Thresholds
error = 100
threshold = 1 # Percentage error
n = 0

# Arrays
T_prop_history = np.array([T_prop_initial]) # start of temperature array
time = np.array([0]) # start of time array

while error > threshold and n < max_iterations:
    T_prop_history = np.append(T_prop_history, Temperature_Iteration(Q, R_air, R_prop, T_air, T_prop_history[n])) #Next temperature step
    time = np.append(time, time[n]+t) #Timestep saved
    error = ((T_prop_final - T_prop_history[n+1]) / T_prop_final) * 100
    n += 1

# In the units we want
time_minute = time/60
T_prop_history_Fahrenheit = (9/5) * (T_prop_history - 273.15) + 32
T_prop_final_Fahrenheit = (9/5) * (T_prop_final - 273.15) + 32

# Generate 300 points from min(x) to max(x)
x_smooth = np.linspace(time_minute.min(), time_minute.max(), 300)

# Create spline of order 3 (cubic)
y_smooth = make_interp_spline(time_minute, T_prop_history_Fahrenheit)(x_smooth)

plt.figure(dpi=1000)
# Plot the smooth curve
plt.plot(x_smooth, y_smooth, color='red', label="Bottle Temperature")  # color is optional
plt.hlines(T_prop_final_Fahrenheit, time_minute.min(), time_minute.max(), linestyles="--", colors="orange", label = "Steady-State Temperature")
# Plot Graph
#plt.scatter(time_minute, T_prop_history)
plt.xlabel("Time (min)")
plt.ylabel("Tempearture (F)")
plt.title("Bottle Heating Time")
plt.legend(loc="lower right")

