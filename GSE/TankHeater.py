#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 12:44:53 2024

@author: ryansmithers
"""

# System Conditions
h_air = 25 #W/m^2-K
h_prop = 72 #W/m^2-K
T_set_imperial = 20 #F
t_total = 5 #hrs
T_sustained_imperial = 85 #f

# Numerical Thresholds
threshold = 1 # Percentage error
time_step = 2 #s
max_iterations = 1000

# Converters
lb2kg = 0.4536
in2m = 0.0254


# Clear all variables
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
plt.close('all')

class Bottle:
    def __init__(self, name, d_outer_imperial, L_imperial, m_prop_imperial, m_bottle_imperial, cp_bottle, cp_prop, Q):
        self.name = name
        self.d_outer = d_outer_imperial * in2m
        self.L = L_imperial * in2m
        self.m_bottle = m_bottle_imperial * lb2kg
        self.m_prop = m_prop_imperial *lb2kg
        self.cp_bottle = cp_bottle
        self.cp_prop = cp_prop
        self.Q = Q
        self.R_air = R_Conv_Cylinder(self.d_outer, h_air, self.L) #K/W
        self.R_prop = R_Conv_Cylinder(self.d_outer, h_prop, self.L) #K/W.... using outer diameter
        self.T_prop_final = self.Q * self.R_air + T_air
        self.Q_sustained = (T_sustained - T_air) / self.R_air
        self.time = 0 #s... time the graph will run to, solved later
        
        # Array Initialization
        self.T_prop_history = np.zeros(max_iterations)
        self.T_bottle_history = np.zeros(max_iterations)
        self.q_inf = np.zeros(max_iterations)
        self.q_prop = np.zeros(max_iterations) 
        self.q_bottle = np.zeros(max_iterations) 
        
        # Solver
        self.solve_for_temperature()
        
        
    def solve_for_temperature(self):
        # Arrays: Initial Conditions
        self.T_prop_history[0] = T_prop_initial
        self.T_bottle_history[0] = T_bottle_initial
        
        # Preparing Iteration
        n = 0
        error = 100
        
        while error > threshold and n+1 < max_iterations:
            #Solve for heat flow values
            self.q_inf[n] = (self.T_bottle_history[n] - T_air) / self.R_air
            self.q_prop[n] = (self.T_bottle_history[n] - self.T_prop_history[n]) / self.R_prop
            self.q_bottle[n] = self.Q - self.q_inf[n] - self.q_prop[n]
            
            #Solve for temperature change
            dT_prop = (time_step * self.q_prop[n]) / (self.m_prop * self.cp_prop)
            dT_bottle = (time_step * self.q_bottle[n]) / (self.m_bottle * self.cp_bottle)
            
            # New Temperature
            self.T_prop_history[n+1] = self.T_prop_history[n] + dT_prop
            self.T_bottle_history[n+1] = self.T_bottle_history[n] + dT_bottle
            
            # Check Error: Ends on n+1 if condition met
            error = ((self.T_prop_final - self.T_prop_history[n+1]) / self.T_prop_final) * 100
            
            # Iterate Count
            n += 1
        
        # Process Data
        self.time = np.arange(0, time_step * (n+1), time_step)
        
        # Remove the trailing zeros in all the arrays
        self.T_prop_history = self.T_prop_history[:n+1]
        self.T_bottle_history = self.T_bottle_history[:n+1]
        self.q_inf = self.q_inf[:n+1]
        self.q_prop = self.q_prop[:n+1]
        self.q_bottle = self.q_bottle[:n+1]
        
        # Add another set of values at end of heat flow so arrays match (could also remove one off time and temp)
        self.q_inf[n] = self.q_inf[n-1]
        self.q_prop[n] = self.q_prop[n-1]
        self.q_bottle[n] = self.q_bottle[n-1]
    
    def Grapher(self):
        # In the units we want
        time_minute = self.time/60
        T_prop_history_Fahrenheit = Kelvin_to_f(self.T_prop_history)
        T_prop_final_Fahrenheit = Kelvin_to_f(self.T_prop_final)
        T_bottle_history_Fahrenheit = Kelvin_to_f(self.T_bottle_history)

        # Generate 300 points from min(x) to max(x)
        x_smooth = np.linspace(time_minute.min(), time_minute.max(), 300)

        # Create spline of order 3 (cubic)
        y_smooth_prop = make_interp_spline(time_minute, T_prop_history_Fahrenheit)(x_smooth)
        y_smooth_bottle = make_interp_spline(time_minute, T_bottle_history_Fahrenheit)(x_smooth)
        y_smooth_q_inf = make_interp_spline(time_minute, self.q_inf)(x_smooth)
        y_smooth_q_prop = make_interp_spline(time_minute, self.q_prop)(x_smooth)
        y_smooth_q_bottle = make_interp_spline(time_minute, self.q_bottle)(x_smooth)


        plt.figure(dpi=1000)

        # Plot Temperature Curves
        plt.plot(x_smooth, y_smooth_bottle, color='red', label="Bottle Temperature")  # color is optional
        plt.plot(x_smooth, y_smooth_prop, color='blue', label="Prop Temperature")  # color is optional
        #plt.hlines(T_prop_final_Fahrenheit, time_minute.min(), time_minute.max(), linestyles="--", colors="orange", label = "Steady-State Temperature")
        plt.xlabel("Time (min)")
        plt.ylabel("Tempearture (F)")
        plt.title(f"{self.name} Heating Time")
        plt.legend(loc="lower right")


        # Plot Heat Flow Curves
        plt.figure(dpi=1000)
        plt.plot(x_smooth, y_smooth_q_bottle, color='red', label="Q into Bottle")  # color is optional
        plt.plot(x_smooth, y_smooth_q_prop, color='blue', label="Q into Prop")  # color is optional
        plt.plot(x_smooth, y_smooth_q_inf, color='orange', label="Q into Air")  # color is optional
        #plt.hlines(T_prop_final_Fahrenheit, time_minute.min(), time_minute.max(), linestyles="--", colors="orange", label = "Steady-State Temperature")
        plt.xlabel("Time (min)")
        plt.ylabel("Heat Flow Rate (W)")
        plt.title(f"{self.name} Heat Flows")
        plt.legend(loc="upper right")

def f_to_Kelvin(f):
    return ((f-32) * (5/9)) + 273

def Kelvin_to_f(K):
    return (9/5) * (K - 273.15) + 32

def R_Cond_Cylinder (r2, r1, k, L):
    return math.log(r2/r1) / (2*math.pi*k*L)

def R_Conv_Cylinder (d, h, L):
    return 1/(h*math.pi*d*L)



# Conversions
T_set = f_to_Kelvin(T_set_imperial) #K
T_sustained = f_to_Kelvin(T_sustained_imperial) #K

# Set starting temperature
T_air = T_prop_initial = T_bottle_initial = T_set #K

# Tank Dimensions
Ethane = Bottle("Ethane", 4, 16.75, 2, 9, 490, 1800, 200)
Nitrous = Bottle("Nitrous", 5.25, 16.75, 5, 8.25, 910, 1700, 240)

# Graphs
Ethane.Grapher()
Nitrous.Grapher()