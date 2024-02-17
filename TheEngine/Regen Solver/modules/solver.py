# general
import os
import numpy as np
import pandas as pd

# thermo and prop
import cantera as ct
from CoolProp.CoolProp import PropsSI

# numerical methods
import scipy.optimize

# visualization / other
import matplotlib.pyplot as plt

plotfolder = "solverplots/"
solvedplotfolder = plotfolder + "solved/"
analysisplotfolder = plotfolder + "analysis/"

''' Equations and Functions Used in Solver '''
# define exhaust gas transport properties - using Cantera
def cp_exhaust(T, p):
    gas.TP = T, p
    return gas.cp # [J/kg K]

def mu_exhaust(T, p):
    gas.TP = T, p
    return gas.viscosity # [Pa s]

def Pr_exhaust(T, p):
    gas.TP = T, p
    return gas.cp * gas.viscosity / gas.thermal_conductivity # Definition of Prandtl number

def rho_exhaust(T, p):
    gas.TP = T, p
    return gas.density # [kg/m^3]

def free_vel_exhaust(T, p, M):
    gas.TP = T, p
    return gas.sound_speed * M # [m/s]


# defining function to get P_sat and T_sat
coolant_name = 'NITROUSOXIDE'

def P_sat(T):
    return PropsSI('P','T', T, 'Q', 0, coolant_name)

def T_sat(P):
    return PropsSI('T','P', P, 'Q', 0, coolant_name)


def rho_f(P):
    return PropsSI('DMASS','P', P, 'Q', 0, coolant_name)

def rho_g(P):
    return PropsSI('DMASS','P', P, 'Q', 1, coolant_name)

def rho_coolant_phase(P, X):
    return X * rho_g(P) + (1 - X) * rho_f(P)


def cp_f(P):
    return PropsSI('CPMASS','P', P, 'Q', 0, coolant_name)

def cp_g(P):
    return PropsSI('CPMASS','P', P, 'Q', 1, coolant_name)

def cp_coolant_phase(P, X):
    return X * cp_g(P) + (1 - X) * cp_f(P)


# first we need to define the data that will produce our functional relation
T = np.array([260, 270, 280, 290, 300, 307])
P = np.array([P_sat(260), P_sat(270), P_sat(280), P_sat(290), P_sat(300), P_sat(307)])
k_f_data = np.array([0.12, 0.115, 0.095, 0.08, 0.075, 0.07])
k_g_data = np.array([0.02, 0.023, 0.025, 0.027, 0.03, 0.065])

coeff_k_f = np.polyfit(P, k_f_data, 1)
coeff_k_g = np.polyfit(P, k_g_data, 1)

# define the thermal conductivity as a function of temperature
def k_f(P):
    return np.polyval(coeff_k_f, P)

def k_g(P):
    return np.polyval(coeff_k_g, P)

def k_coolant_phase(P, X):
    return (1-X)*k_f(P) + X*k_g(P)


# first we need to define the data that will produce our functional relation
T = np.array([260, 270, 280, 290, 300, 307])
P = np.array([P_sat(260), P_sat(270), P_sat(280), P_sat(290), P_sat(300), P_sat(307)])
mu_f_data = np.array([1.05e-4, 9.5e-5, 8.3e-5, 7e-5, 6e-5, 4e-5])
mu_g_data = np.array([1.5e-5, 1.6e-5, 1.7e-5, 1.8e-5, 2e-5, 2.3e-5])

coeff_mu_f = np.polyfit(P, mu_f_data, 1)
coeff_mu_g = np.polyfit(P, mu_g_data, 1)

# define the dynamic viscosity as a function of temperature
def mu_f(P):
    return np.polyval(coeff_mu_f, P)

def mu_g(P):
    return np.polyval(coeff_mu_g, P)

def mu_coolant_phase(P, X):
    return (1-X)*mu_f(P) + X*mu_g(P)


def Pr_coolant_phase(P, X):
    return cp_coolant_phase(P, X) * mu_coolant_phase(P, X) / k_coolant_phase(P, X)


# define function to calculate vapor quality X (in CoolProp, Q) of coolant (given enthalpy h and pressure P)
def X_coolant(h, P):
    h_f = PropsSI('HMASS','P', P, 'Q', 0, coolant_name)
    h_g = PropsSI('HMASS','P', P, 'Q', 1, coolant_name)
    return (h - h_f) / (h_g - h_f)


# define functions to get coolant properties as functions of T and P
def cp_coolant_super(T, P):
    # note CoolProp might not like retrieving values very close to saturation - then use saturation values if get ValueError
    try:
        cp = PropsSI('CPMASS','P', P, 'T', T, coolant_name)
    except ValueError:
        cp = PropsSI('CPMASS','P', P, 'Q', 1, coolant_name)
        print("Warning: CoolProp failed to find cp at T = {} K and P = {} Pa. Using saturated gas estimate.".format(T, P))
    return cp
def rho_coolant_super(T, P):
    try:
        rho = PropsSI('DMASS','P', P, 'T', T, coolant_name)
    except ValueError:
        rho = PropsSI('DMASS','P', P, 'Q', 1, coolant_name)
        print("Warning: CoolProp failed to find rho at T = {} K and P = {} Pa. Using saturated gas estimate.".format(T, P))
    return rho

set_mu = 0 
set_k = 0 
def mu_coolant_super(T, P):
    global set_mu
    if set_mu == 0: # the first time this function is called, set set_mu to the value of mu_g at that pressure
        set_mu = mu_g(P)
        #print("Viscosity set at", set_mu, "Pa s")
    return set_mu # not defined - get estimates as if saturated gas state
def k_coolant_super(T, P):
    global set_k
    if set_k == 0: # the first time this function is called, set set_mu to the value of mu_g at that pressure
        set_k = k_g(P)
        #print("Conductivity set at", set_k, "W/m K")
    return set_k # not defined - get estimates as if saturated gas state

def Pr_coolant_super(T, P):
    return cp_coolant_super(T, P) * mu_coolant_super(T, P) / k_coolant_super(T, P)


# film coefficient functions
def h_gas_bartz(D, cp, mu, Pr, rho, v, rho_am, mu_am, mu0):
    """
    Bartz equation, using Equation (8-23) from page 312 of RPE 7th edition (Reference [1]). 'am' refers to the gas being at the 'arithmetic mean' of the wall and freestream temperatures.

    Args:
        D (float): Gas flow diameter (m)
        cp_inf (float): Specific heat capacity at constant pressure for the gas, in the freestream
        mu_inf (float): Absolute viscosity in the freestream
        Pr_inf (float): Prandtl number in the freestream
        rho_inf (float): Density of the gas in the freestream
        v_inf (float): Velocity of the gas in in the freestream
        rho_am (float): Density of the gas, at T = (T_wall + T_freestream)/2, P - exhaust static pressure
        mu_am (float): Absolute viscosity of the gas, at T = (T_wall + T_freestream)/2, P - exhaust static pressure
        mu0 (float): Absolute viscosity of the gas under stagnation conditions.
        
    Returns:
        float: Convective heat transfer coefficient (W/m2/K), h, for the exhaust gas side (where q = h(T - T_inf)).
    """

    return (0.026/D**0.2) * (cp*mu**0.2)/(Pr**0.6) * (rho * v)**0.8 * (rho_am/rho) * (mu_am/mu0)**0.2

def h_coolant_dittus_boelter(rho, v, Dh, mu, Pr, k):
    """Dittus-Boelter equation for convective heat transfer coefficient.

    Args:
        rho (float): Coolant bulk density (kg/m^3).
        V (float): Coolant bulk velocity (m/s)
        D (float): Hydraulic diameter of pipe (m)
        mu (float): Coolant bulk viscosity (Pa s)
        Pr (float): Coolant bulk Prandtl number
        k (float): Coolant thermal conductivity

    Returns:
        float: Convective heat transfer coefficient (W/m2/K)
    """
    Re = rho*v*Dh/mu
    Nu = 0.023*Re**(4/5)*Pr**0.4

    return Nu*k/Dh

# thermal resistance functions for a dx slice 
def R_wall(r, th, k_wall, dx):
    """Thermal resistance of a wall of inner radius r, thickness th, thermal conductivity k, slice width dx.

    Args:
        r (float): Radius (m)
        th (float): Wall thickness (m)
        k (float): Thermal conductivity (W/m/K)

    Returns:
        float: Thermal resistance for wall (K/W)
    """

    return np.log((r + th)/r)/(2*np.pi*k_wall*dx)

def R_convective(r, h, dx): # used for both coolant and exhaust gas
    """Thermal resistance of a convective boundary layer of slice length dx and convective heat transfer coefficient h.

    Args:
        r (float): Local radius of the surface at which convection occurs (m)
        h (float): Convective heat transfer coefficient (W/m2/K)
        dx (float): Boundary layer slice length (m)

    Returns:
        float: Thermal resistance for convective boundary layer (K/W)
    """

    return 1/(2*np.pi*h*r*dx)



def Dh_rectangular(w, h):
    """Hydraulic diameter of a rectangular channel.

    Args:
        w (float): Width of channel (m)
        h (float): Height of channel (m)

    Returns:
        float: Hydraulic diameter of channel (m)
    """
    return 2*w*h/(w+h)
def dp_darcy_weisbach(f_D, rho, v, Dh, drdx, dx):
    """Darcy-Weisbach equation for pressure drop in a pipe.

    Args:
        f_D (float): Darcy friction factor
        rho (float): Density of fluid (kg/m^3)
        v (float): Velocity of fluid (m/s)
        Dh (float): Hydraulic diameter of pipe (m)

    Returns:
        float: Pressure drop (Pa)
    """
    return f_D*rho*(v**2) / (2*Dh) * np.sqrt(drdx**2 + 1) * dx

def reynolds(Dh, v, rho, mu):
    """Reynolds number for a pipe.

    Args:
        Dh (float): Hydraulic diameter of pipe (m)
        v (float): Velocity of fluid (m/s)
        rho (float): Density of fluid (kg/m^3)
        mu (float): Dynamic viscosity of fluid (Pa s)

    Returns:
        float: Reynolds number
    """
    return rho*v*Dh/mu

def mach_number(gamma, v, T, MW):
    """Mach number for a gas.

    Args:
        gamma (float): Specific heat ratio
        v (float): Velocity of gas (m/s)
        T (float): Temperature of gas (K)
        MW (float): Molecular weight of gas (kg/kmol)

    Returns:
        float: Mach number
    """
    return v / np.sqrt(gamma*8.314*T/MW)


# defining isentropic flow equations for static T and P of exhaust gas given M, kc, T0, and P0
def isen_T(M, kc, T0):
    return T0*(1 + ((kc - 1)/2) * M**2)**(-1)

def isen_P(M, kc, P0):
    return P0*(1 + ((kc - 1)/2) * M**2)**((-kc)/(kc - 1))


''' Define Node Object '''
# Define Node Object
class Node:

    # constructor and setters
    def __init__(self):
        self.x = 0
        self.r = 0
        self.drdx = 0
        self.A = 0
        self.M = 0
        self.T = 0
        self.P = 0
        self.OD = 0
        self.w = 0
        self.h = 0
        self.A_c = 0
        self.Dh_c = 0
        self.T_hw = 0
        self.T_cw = 0
        self.T_c = 0
        self.P_c = 0
        self.v_c = 0
        self.dv_c = 0
        self.Re_c = 0
        self.q = 0
        self.Q = 0
        self.t_w = 0
        self.X = 0
        self.H = 0
        self.rho_c = 0
        self.k_c = 0
        self.mu_c = 0
        self.cp_c = 0
        
    
    def set_x(self, x):
        self.x = x
    
    def set_r(self, r):
        self.r = r

        # implicilty calculate area from r
        self.A = np.pi * r**2

    def set_drdx(self, drdx):
        self.drdx = drdx
        
    def set_M(self, M):
        self.M = M

    def set_T(self, T):
        self.T = T

    def set_P(self, P):
        self.P = P

    def set_OD(self, OD):
        self.OD = OD

    def set_w(self, w):
        self.w = w
    
    def set_h(self, h):
        self.h = h
        self.A_c = self.w * self.h # setting height implicitly sets area of coolant channel
        self.Dh_c = Dh_rectangular(self.w, self.h) # setting height implicitly sets hydraulic diameter of coolant channel

    def set_T_hw(self, T_hw):
        self.T_hw = T_hw
    
    def set_T_cw(self, T_cw):
        self.T_cw = T_cw

    def set_T_c(self, T_c):
        self.T_c = T_c

    def set_P_c(self, P_c):
        self.P_c = P_c

    def set_v_c(self, v_c):
        self.v_c = v_c
    
    def set_dv_c(self, dv_c):
        self.dv_c = dv_c
    
    def set_Re_c(self, Re_c):
        self.Re_c = Re_c

    def set_q(self, q):
        self.q = q

    def set_Q(self, Q):
        self.Q = Q

    def set_t_w(self, t_w):
        self.t_w = t_w

    def set_X(self, X):
        self.X = X

    def set_H(self, H):
        self.H = H
        
    def set_rho_c(self, rho_c):
        self.rho_c = rho_c
        
    def set_k_c(self, k_c):
        self.k_c = k_c
    
    def set_mu_c(self, mu_c):
        self.mu_c = mu_c
        
    def set_cp_c(self, cp_c):
        self.cp_c = cp_c


    # special methods
    def print_coolant_parameters(self):
        # Define the headers for the table
        headers = ["Parameter", "Value"]

        # Define the data for the table
        data = [
            ("Axial Position (x)", self.x),
            ("Channel Width (w)", self.w),
            ("Channel Height (h)", self.h),
            ("Coolant Bulk Temperature (T_c)", self.T_c),
            ("Coolant Pressure (P_c)", self.P_c),
            ("Coolant Velocity (v_c)", self.v_c),
            ("Coolant Velocity Change (dv_c)", self.dv_c),
            ("Heat Flow through Walls (Q)", self.Q),
            ("Temperature Drop along Wall (T_hw-T_cw)", self.T_hw - self.T_cw),
            ("Wall Thickness (t_w)", self.t_w),
            ("Vapor Quality", self.X)
        ]

        # Call the helper method to print the table
        self._print_table(headers, data)

    # Helper method to print the table
    def _print_table(self, headers, data):
        # Calculate the maximum width for each column
        max_widths = [max(len(str(item)) for item in column) for column in zip(headers, *data)]

        # Print the headers
        print(" | ".join(f"{header:<{width}}" for header, width in zip(headers, max_widths)))

        # Print the separator line
        print("-" * (sum(max_widths) + len(headers) * 3 - 1))

        # Print the data rows
        for row in data:
            print(" | ".join(f"{item:<{width}}" for item, width in zip(row, max_widths)))



''' Defining Active Functions '''

def setupDirectories():
    if not os.path.exists(plotfolder):
        os.mkdir(plotfolder)
    if not os.path.exists(solvedplotfolder):
        os.mkdir(solvedplotfolder)
    if not os.path.exists(analysisplotfolder):
        os.mkdir(analysisplotfolder)

def defineParameters(P0, Ti, P_ci, X_ci, x_off, T_hw, T_cw, N_channels, t_rib, k_wall, N):
    global def_P0, def_Ti, def_P_ci, def_X_ci, def_x_off, def_T_hw, def_T_cw, def_N_channels, def_t_rib, def_k_wall, def_N
    def_P0 = P0
    def_Ti = Ti
    def_P_ci = P_ci
    def_X_ci = X_ci
    def_x_off = x_off
    def_T_hw = T_hw
    def_T_cw = T_cw
    def_N_channels = N_channels
    def_t_rib = t_rib
    def_k_wall = k_wall
    def_N = N


USE_LINEAR_TEMP = False
def defineLinearTemperature(T_hwe, T_hwt, T_hwi, T_cwe, T_cwt, T_cwi):
    global def_T_hwe, def_T_hwt, def_T_hwi, def_T_cwe, def_T_cwt, def_T_cwi, USE_LINEAR_TEMP
    def_T_hwe = T_hwe
    def_T_hwt = T_hwt
    def_T_hwi = T_hwi
    def_T_cwe = T_cwe
    def_T_cwt = T_cwt
    def_T_cwi = T_cwi
    
    USE_LINEAR_TEMP = True

def run(verbose = False, plotExhaustGases = False, plotCoolingChannels = False, printNodes = False):
    print("Running Solver...\n")
    global gas, nodes, nodesc, x, r, th, M_coolant, x_n, r_n, o_mdot, T_hw_n, T_cw_n

    # read geometry csv file to pandas
    df = pd.read_csv("enginefiles/enginegeometry.csv", sep=",")

    # extract geometry 
    x = df["x [m]"].to_list()
    r = df["r [m]"].to_list()

    # read parameter csv file to pandas
    df = pd.read_csv("enginefiles/engineparameters.csv", sep=",")

    # extract geometry parameters
    A_t = df["A_t [m^2]"].to_list()[0]
    A_e = df["A_e [m^2]"].to_list()[0]
    A_c = df["A_c [m^2]"].to_list()[0]
    x_c = df["x_c [m]"].to_list()[0]
    x_t = df["x_t [m]"].to_list()[0]

    # read mass flow csv file to pandas
    df = pd.read_csv("enginefiles/enginemassflow.csv", sep=",")

    # extract mass flow parameters
    OF = df["OF used [-]"].to_list()[0] # OF ratio used in engine sizing
    o_mdot = df["o_mdot [kg/s]"].to_list()[0] # oxidizer mass flow rate [kg/s]
    f_mdot = df["f_mdot [kg/s]"].to_list()[0] # fuel mass flow rate [kg/s] 

    #print("x:", x, "\nr:", r, "\nA_t:", A_t, "\nx_t:", x_t, "\nA_e:", A_e, "\nA_c:", A_c, "\nx_c:", x_c, "\no_mdot:", o_mdot, "\nf_mdot:", f_mdot)
    P0 = def_P0 # First Chamber (Stagnation) Pressure Guess [Pa]
    Ti = def_Ti # First Chamber Inlet Temperature Guess [K]

    # determine coolant inlet temperature
    coolant_name = 'NITROUSOXIDE'
    def_T_ci = PropsSI('T', 'P', def_P_ci, 'Q', def_X_ci, coolant_name)
    H_ci = PropsSI('HMASS','P', def_P_ci, 'Q', def_X_ci, coolant_name)
    X_ci = def_X_ci

    P0 = def_P0 # First Chamber (Stagnation) Pressure Guess [Pa]
    Ti = def_Ti # First Chamber Inlet Temperature Guess [K]

    # Define gas
    gas = ct.Solution('gri30.yaml')          
    mixture = "C2H6:1, N2O:{}".format(OF)    # define mixture via mixture string
    gas.TPY = Ti, P0, mixture              # define state of gas before chemical balance
    gas.equilibrate("HP")                  # equilibrate keeping enthalpy and pressure constant

    # Extract Preliminary Gas Properties
    h0 = gas.h  # gas enthalpy [J]
    T0 = gas.T  # stagnation temperature [K]
    kc = gas.cp / gas.cv # specific heat ratio in chamber
    MW = gas.mean_molecular_weight # average molecular weight of gas [kg / kmol]
    mu0 = gas.viscosity # dynamic viscosity [Pa s]

    # Print Properties
    if verbose:
        print("Results of Gas calculations:")
        print("____________________________________________________")
        print("Enthalpy:", h0, "[J]\nStagnation temperature:", T0, "[K]\nSpecific heat ratio:", kc, "[-]\nMean molecular weight:", MW, "[kg/kmol]")

    T_ci = def_T_ci # [K] - coolant inlet temperature
    P_ci = def_P_ci # [Pa] - coolant inlet pressure 

    x_f = max(x) - def_x_off # [m] - coordinate of cooling jacket end

    #print("x_f:", x_f)

    # define function for thickness of chamber wall th(x)

    # Note: Currently Using Constant Thickness
    const_thickness = 0.005 # [m]
    def variable_thickness(x):
        return const_thickness
    th = np.array([variable_thickness(_x) for _x in x ]) # thickness of chamber wall [m] 


    k_wall = def_k_wall # [W / m K] - very approximate thermal conductivity of Inconel chamber wall

    # define amount of nodes
    N = def_N

    # calculate node axial increment
    dx = max(x) / N
    
    

    # finding x positons for various pairs of points
    points12x = x[:2] # points 1 and 2
    points23x = x[1:3] # points 2 and 3
    points34x = x[2:4] # points 3 and 4

    # finding r positons for various pairs of points
    points12r = r[:2] # points 1 and 2
    points23r = r[1:3] # points 2 and 3
    points34r = r[2:4] # points 3 and 4

    # finding slopes and intercepts for various pairs of points
    slope12, intercept12 = np.polyfit(points12x, points12r, 1) # points 1 and 2
    slope23, intercept23 = np.polyfit(points23x, points23r, 1) # points 2 and 3
    slope34, intercept34 = np.polyfit(points34x, points34r, 1) # points 3 and 4


    x_n = np.linspace(x[0], x[-1], N) # x positions of nodes
    r_n = np.zeros_like(x_n) # r positions of nodes
    OD_n = np.zeros_like(x_n) # outer diameter of chamber wall of nodes
    drdx_n = np.zeros_like(x_n) # slope of channel walls at nodes

    for i, el in enumerate(x_n):
        slope = 0
        # find r_n
        if el <= x[1]:
            node_r_n = slope12 * el + intercept12
            # note: for interval 1-2 slope is 0
        elif el <= x[2]:
            node_r_n = slope23 * el + intercept23
            slope = slope23
        else:
            node_r_n = slope34 * el + intercept34
            slope = slope34

        # set r_n
        r_n[i] = node_r_n

        # set wallOD
        OD_n[i] = 2 * (node_r_n + variable_thickness(el))

        # set slope
        drdx_n[i] = slope

    # print all linear equations obtained
    if verbose:
        print("\nNode Count and Geometry:")
        print("____________________________________________________")
        print("Amount of nodes:", N, "\nAxial increment:", dx * 1000, "mm")
        print("r = {}x + {}".format(slope12, intercept12))
        print("r = {}x + {}".format(slope23, intercept23))
        print("r = {}x + {}".format(slope34, intercept34))
        print("Slopes of channel walls at nodes:", drdx_n)

    # create size N array of Nodes
    nodes = np.empty(N, dtype=Node)

    # set x and r values of each node
    for i, el in enumerate(nodes):
        el = Node()
        el.set_x(x_n[i])
        el.set_r(r_n[i])
        el.set_OD(OD_n[i])
        el.set_drdx(drdx_n[i])
        nodes[i] = el


    # express r_A as a function of k and M
    def A_At(M, k):
        return 1/M * ( (2 + (k-1) * M**2) / (k + 1) )**( (k+1) / (2 * (k-1) ) )

    # define function to find M when flow is subsonic
    def get_M_subsonic(A, A_t, kc):
        # define function to give give the numerical solver
        def func_to_solve(M): # here Mach number M is the "x" variable we want scipy.optimize.root to solve for - we want to find the coordinatate of x such that the function is zero i.e. the root
            return  A/A_t - A_At(M, kc) # we want to minimize the difference / find the root (when the difference is zero)
        
        return scipy.optimize.root_scalar(func_to_solve, bracket = [1e-5, 1], x0 = 0.5).root # the bracket is the range of values to search for the root, x0 is the initial guess

    def get_M_supersonic(A, A_t, kc):
        # define function to give give the numerical solver
        def func_to_solve(M): # here Mach number M is the "x" variable we want scipy.optimize.root to solve for - we want to find the coordinatate of x such that the function is zero i.e. the root
            return  A/A_t - A_At(M, kc) # we want to minimize the difference / find the root (when the difference is zero)
        
        return scipy.optimize.root_scalar(func_to_solve, bracket = [1, 500], x0 = 1).root # the bracket is the range of values to search for the root, x0 is the initial guess


    # find Mach number at each node
    for node in nodes:
        if node.x < x_c:
            node.set_M(get_M_subsonic(A_c, A_t, kc)) # if node is in chamber cylindrical section use approximation that M cylindrical is M at nozzle inlet
        elif node.x < x_t:
            node.set_M(get_M_subsonic(node.A, A_t, kc)) # if node is in convergent section, M is subsonic
        else:
            node.set_M(get_M_supersonic(node.A, A_t, kc)) # if node is in divergent section, M is supersonic


    # plot Mach number at each node
    if plotExhaustGases:
        fig, axs = plt.subplots(figsize = (12, 6))
        fig.set_facecolor('white')
        axs.plot(x, r, color = "k") # upper innner wall
        axs.plot(x, -np.array(r), color = "k") # lower inner wall
        axs.plot(x, np.add(r, th), color = "k") # upper outer wall
        axs.plot(x, np.subtract(-np.array(r), th), color = "k") # lower outer wall
        axs.grid()
        axs.set_xlabel("Axial Position (m)")
        axs.set_ylabel("Radius (m)")
        axs.set_title("Exhaust Gas Mach Number Across Chamber")
        axs.set_aspect('equal')
        axs2 = axs.twinx()
        axs2.grid(color = "m", linestyle = "--", linewidth = 0.5)
        axs2.plot(x_n, [node.M for node in nodes], color = "m")
        axs2.set_ylabel("Mach Number [-]", color = "m")
        plt.savefig(plotfolder + "machvsenginegeometry.png", dpi=300)
        plt.show()


    # populate nodes array with T and P values - find temperature and pressue at each node
    for node in nodes:
        if node.M == 0:  # set T and P to stagnation values if node is in chamber
            node.set_T(T0)
            node.set_P(P0)
        else:
            node.set_T(isen_T(node.M, kc, T0))
            node.set_P(isen_P(node.M, kc, P0))

    # plot Temperature at each node
    if plotExhaustGases:
        fig, axs = plt.subplots(figsize = (12, 6))
        fig.set_facecolor('white')
        axs.plot(x, r, color = "k") # upper innner wall
        axs.plot(x, -np.array(r), color = "k") # lower inner wall
        axs.plot(x, np.add(r, th), color = "k") # upper outer wall
        axs.plot(x, np.subtract(-np.array(r), th), color = "k") # lower outer wall
        axs.grid()
        axs.set_xlabel("Axial Position (m)")
        axs.set_ylabel("Radius (m)")
        axs.set_title("Exhaust Gas Temperature Across Chamber")
        axs.set_aspect('equal')
        axs2 = axs.twinx()
        axs2.grid(color = "coral", linestyle = "--", linewidth = 0.5)
        axs2.plot(x_n, [node.T for node in nodes], color = "coral")
        axs2.set_ylabel("Temperature [K]", color = "coral")
        plt.savefig(plotfolder + "tempvsenginegeometry.png", dpi=300)
        plt.show()

        # plot Pressure at each node
        fig, axs = plt.subplots(figsize = (12, 6))
        fig.set_facecolor('white')
        axs.plot(x, r, color = "k") # upper innner wall
        axs.plot(x, -np.array(r), color = "k") # lower inner wall
        axs.plot(x, np.add(r, th), color = "k") # upper outer wall
        axs.plot(x, np.subtract(-np.array(r), th), color = "k") # lower outer wall
        axs.grid()
        axs.set_xlabel("Axial Position (m)")
        axs.set_ylabel("Radius (m)")
        axs.set_title("Exhaust Gas Pressure Across Chamber")
        axs.set_aspect('equal')
        axs2 = axs.twinx()
        axs2.grid(color = "g", linestyle = "--", linewidth = 0.5)
        axs2.plot(x_n, [node.P for node in nodes], color = "g")
        axs2.set_ylabel("Pressure [Pa]", color = "g")
        plt.savefig(plotfolder + "presvsenginegeometry.png", dpi=300)
        plt.show()


    # find position of last node of cooling jacket
    i_f = 0
    for i, node in enumerate(nodes):
        if node.x > x_f:
            i_f = i - 1 # if this node is larger than the final position, then the previous node is the last node of the cooling jacket
            break
        elif node.x == x_f:
            i_f = i # if this node is equal to the final position, then this node is the last node of the cooling jacket
            break

    nodesc = nodes[:i_f+1] # cooling jacket nodes (i_f + 1 is not included)

    if verbose:
        print("Result of Node Slicing:")
        print("____________________________________________________")
        print("Number of nodes:", len(nodes), "\nNumber of cooling jacket nodes:", len(nodesc), "\nLast cooling jacket node index:", i_f)

    
    # interpolate hot wall temperature T_hw(x) at each node
    # finding x positons for various pairs of points
    points12x = x[:2] # points 1 and 2
    points23x = x[1:3] # points 2 and 3
    points34x = x[2:4] # points 3 and 4

    points12T_hw = [def_T_hwi, def_T_hwi] # points 1 and 2
    points23T_hw = [def_T_hwi, def_T_hwt] # points 2 and 3
    points34T_hw = [def_T_hwt, def_T_hwe] # points 3 and 4

    # finding slopes and intercepts for various pairs of points
    slope12, intercept12 = np.polyfit(points12x, points12T_hw, 1) # points 1 and 2
    slope23, intercept23 = np.polyfit(points23x, points23T_hw, 1) # points 2 and 3
    slope34, intercept34 = np.polyfit(points34x, points34T_hw, 1) # points 3 and 4

    # print all linear equations obtained
    if verbose:
        print("Result of Hot Wall Temperature Calculation:")
        print("____________________________________________________")
        print("T_hw = {}x + {}".format(slope12, intercept12))
        print("T_hw = {}x + {}".format(slope23, intercept23))
        print("T_hw = {}x + {}".format(slope34, intercept34))

    # get all T_hw values at nodes
    T_hw_n = np.zeros_like(x_n) # hot wall temperature at nodes
    for i, el in enumerate(x_n):
        slope = 0
        # find T_hw_n
        if el <= x[1]:
            node_T_hw_n = slope12 * el + intercept12
            # note: for interval 1-2 slope is 0
        elif el <= x[2]:
            node_T_hw_n = slope23 * el + intercept23
            slope = slope23
        else:
            node_T_hw_n = slope34 * el + intercept34
            slope = slope34

        # set T_hw_n
        T_hw_n[i] = node_T_hw_n

    for i, node in enumerate(nodes):
        if USE_LINEAR_TEMP:
            node.set_T_hw(T_hw_n[i])
        else:
            node.set_T_hw(def_T_hw)

    # finding x positons for various pairs of points
    points12x = x[:2] # points 1 and 2
    points23x = x[1:3] # points 2 and 3
    points34x = x[2:4] # points 3 and 4

    points12T_cw = [def_T_cwi, def_T_cwi] # points 1 and 2
    points23T_cw = [def_T_cwi, def_T_cwt] # points 2 and 3
    points34T_cw = [def_T_cwt, def_T_cwe] # points 3 and 4

    # finding slopes and intercepts for various pairs of points
    slope12, intercept12 = np.polyfit(points12x, points12T_cw, 1) # points 1 and 2
    slope23, intercept23 = np.polyfit(points23x, points23T_cw, 1) # points 2 and 3
    slope34, intercept34 = np.polyfit(points34x, points34T_cw, 1) # points 3 and 4

    # print all linear equations obtained
    if verbose:
        print("Result of Cold Wall Temperature Calculation:")
        print("____________________________________________________")
        print("T_cw = {}x + {}".format(slope12, intercept12))
        print("T_cw = {}x + {}".format(slope23, intercept23))
        print("T_cw = {}x + {}".format(slope34, intercept34))

    # get all T_cw values at nodes
    T_cw_n = np.zeros_like(x_n) # cold wall temperature at nodes
    for i, el in enumerate(x_n):
        slope = 0
        # find T_cw_n
        if el <= x[1]:
            node_T_cw_n = slope12 * el + intercept12
            # note: for interval 1-2 slope is 0
        elif el <= x[2]:
            node_T_cw_n = slope23 * el + intercept23
            slope = slope23
        else:
            node_T_cw_n = slope34 * el + intercept34
            slope = slope34

        # set T_cw_n
        T_cw_n[i] = node_T_cw_n

    # populate nodes with target cold wall temperatures
    for i, node in enumerate(nodes):
        if USE_LINEAR_TEMP:
            node.set_T_cw(T_cw_n[i])
        else:
            node.set_T_cw(def_T_cw)


    # define cooling jacket geometry parameters
    N_channels = def_N_channels # number of cooling channels
    t_rib = def_t_rib # thickness of ribs [m]


    # calculate length of ribs at all envelope sections - for N_channels cooling channels there will be N_channels ribs (when counting going around the circumference)
    C = np.pi * OD_n # circumference of chamber outer wall at each node
    L_rt = [(N_channels) * t_rib for i in range(N)] # total circumference length occupied by ribs at each node
    w = (C - L_rt) / N_channels # width of each cooling channel

    if min(w) < 0:
        print("ERROR: Cooling channel width is negative at one or more nodes. Please check cooling jacket geometry parameters.")


    # print results
    if verbose:
        print("Result of Cooling Channel Width Calculation:")
        print("____________________________________________________")
        print("Range of cooling channel width:", min(w) * 1000.0, "to", max(w) * 1000.0, "[mm]")


     # plot results: length occupied by rib, length occupied by walls and total envelope length (circumference) at each node
    if plotCoolingChannels:
        fig, axs = plt.subplots(figsize = (12, 10))
        fig.set_facecolor('white')
        axs.plot(x_n, C, color = "r", label = "Total Envelope (Circumference) Length")
        axs.plot(x_n, C - L_rt, color = "g", label = "Total Cooling Channel Occupied Length")
        axs.plot(x_n, L_rt, color = "k", label = "Total Rib Occupied Length")
        axs.plot(x_n, w, color = "b", label = "Cooling Channel Width")
        axs.grid()
        axs.set_xlabel("Axial Position (m)")
        axs.set_ylabel("Length (m)")
        axs.set_title("Length Distributions of Cooling Jacket along Chamber")
        axs.legend()
        plt.show()


        # visualize channel geometry by plotting channel walls on top of chamber outer wall envelope
        fig, axs = plt.subplots(figsize = (12, 6))
        fig.set_facecolor('white')
        axs.grid()
        # plot chamber envelope
        axs.plot(x_n, (OD_n * np.pi) / 2, color = "k", linewidth = 2)
        axs.plot(x_n, -(OD_n * np.pi) / 2, color = "k", linewidth = 2)

        #axs.fill_between(x_n, (OD_n * np.pi) / 2, -(OD_n * np.pi) / 2, color = "gray") # color in ribs
        axs.set_xlabel("Axial Position (m)")
        axs.set_ylabel("Outer Wall Circumference (m)")
        axs.set_title("Cooling Channel Geometry over Chamber")
        for i in range(N_channels):
            axs.fill_between(x_n, - (OD_n * np.pi) / 2 + i * (t_rib + w), - (OD_n * np.pi) / 2 + i * (t_rib + w) + w, color = "b") # at the x-th position with respective outer diameter OD and channel width w plot the i-th channel
        plt.show()
        plt.savefig(plotfolder + "coolingchannelenvelope.png", dpi=300)


    # assign each node the width of the cooling channels at that node
    for i, el in enumerate(nodes):
        el.set_w(w[i])

    # reverse order of nodes array
    nodes = nodes[::-1]


    # for node 0
    Dh = 2 * nodes[0].r # diameter of gas flow at node 0
    cp = cp_exhaust(nodes[0].T, nodes[0].P) # specific heat capacity at constant pressure for the gas, in the freestream
    mu = mu_exhaust(nodes[0].T, nodes[0].P) # absolute viscosity in the freestream
    Pr = Pr_exhaust(nodes[0].T, nodes[0].P) # Prandtl number in the freestream
    rho = rho_exhaust(nodes[0].T, nodes[0].P) # density of the gas in the freestream
    v = free_vel_exhaust(nodes[0].T, nodes[0].P, nodes[0].M) # velocity of the gas in in the freestream
    rho_am = rho_exhaust((nodes[0].T + nodes[0].T_hw) / 2.0, nodes[0].P) # density of the gas, at T = (T_wall + T_freestream)/2
    mu_am = mu_exhaust((nodes[0].T + nodes[0].T_hw) / 2.0, nodes[0].P) # absolute viscosity of the gas, at T = (T_wall + T_freestream)/2

    # for each node in list calculate h_gas and heat flow rate
    for node in nodes:
        Dh = 2 * node.r # diameter of gas flow at node
        cp = cp_exhaust(node.T, node.P) # specific heat capacity at constant pressure for the gas, in the freestream
        mu = mu_exhaust(node.T, node.P) # absolute viscosity in the freestream
        Pr = Pr_exhaust(node.T, node.P) # Prandtl number in the freestream
        rho = rho_exhaust(node.T, node.P) # density of the gas in the freestream
        v = free_vel_exhaust(node.T, node.P, node.M) # velocity of the gas in in the freestream
        rho_am = rho_exhaust((node.T + node.T_hw) / 2.0, node.P) # density of the gas, at T = (T_wall + T_freestream)/2
        mu_am = mu_exhaust((node.T + node.T_hw) / 2.0, node.P) # absolute viscosity of the gas, at T = (T_wall + T_freestream)/2

        # calculate film coefficeint h_gas and heat flow rate per area q
        h_gas = h_gas_bartz(Dh, cp, mu, Pr, rho, v, rho_am, mu_am, mu0) # [W/m^2/K]
        q = h_gas * (node.T - node.T_hw) # [W/m^2]

        # calculate heat flow rate Q at slice
        R_gas = R_convective(node.r, h_gas, dx) # [K/W] - thermal resistance of gas layer slice
        Q = (node.T - node.T_hw) / R_gas # [W]

        # set node values
        node.set_q(q) # [W/m^2]
        node.set_Q(Q) # [W]


    # plot results
    # first need to reverse node list back to original order
    nodes = nodes[::-1]
    if plotExhaustGases:
        fig, axs = plt.subplots(figsize = (12, 6))
        fig.set_facecolor('white')
        axs.plot(x, r, color = "k") # upper innner wall
        axs.plot(x, -np.array(r), color = "k") # lower inner wall
        axs.plot(x, np.add(r, th), color = "k") # upper outer wall
        axs.plot(x, np.subtract(-np.array(r), th), color = "k") # lower outer wall
        axs.grid()
        axs.set_xlabel("Axial Position (m)")
        axs.set_ylabel("Radius (m)")
        axs.set_title("Heat Transfer Rate Across Chamber")
        axs.set_aspect('equal')
        axs2 = axs.twinx()
        axs2.grid(color = "r", linestyle = "--", linewidth = 0.5)
        axs2.plot(x_n, [node.q for node in nodes], color = "r")
        axs2.set_ylabel("Heat Transfer Rate [W/m^2]", color = "r")
        plt.savefig(plotfolder + "heatvsenginegeometry.png", dpi=300)
        plt.show()

    # rereverse node list
    nodes = nodes[::-1]



    # reverse cooling jacket node list
    nodesc = nodesc[::-1]

    # check that first node of cooling jacket is located at end of cooling jacket
    #print("Node 0 position:", nodesc[0].x)

    # make some inital guessees to start iteration process at first node
    h_guess = 0.005 # [m] - initial guess for coolant wall height
    t_w_guess = const_thickness # [m] - initial guess for wall thickness

    # set assumptions for first node
    nodesc[0].set_T_c(T_ci) # set coolant temperature at first node to coolant inlet temperature
    nodesc[0].set_P_c(P_ci) # set coolant pressure at first node to coolant inlet pressure

    # set enthalpy of first node to saturated liquid at 300 K
    nodesc[0].set_H(H_ci)
    nodesc[0].set_X(X_ci)

    # set initial node to node 0
    node0 = nodesc[0]
    
    
    # set properties of coolant at first node
    if node0.X < 1:
        rho_c = rho_coolant_phase(node0.P_c, node0.X) # set coolant density at first node
        k_c = k_coolant_phase(node0.P_c, node0.X) # set coolant conductivity at first node
        mu_c = mu_coolant_phase(node0.P_c, node0.X) # set coolant viscosity at first node
        cp_c = cp_coolant_phase(node0.P_c, node0.X) # set coolant specific heat at first node
    else:
        rho_c = rho_coolant_super(node0.T_c, node0.P_c)
        k_c = k_coolant_super(node0.T_c, node0.P_c)
        mu_c = mu_coolant_super(node0.T_c, node0.P_c)
        cp_c = cp_coolant_super(node0.T_c, node0.P_c)
        
    node0.set_rho_c(rho_c) # set coolant density at first node
    node0.set_k_c(k_c) # set coolant conductivity at first node
    node0.set_mu_c(mu_c) # set coolant viscosity at first node
    node0.set_cp_c(cp_c) # set coolant specific heat at first node


    # define function to find wall thickness such that necessary temperature drop across wall is achieved
    def solve_for_t_w(node):
        def func_to_solve(t_w):
            # find thermal resistance of wall
            R_w = R_wall(node.r, t_w, k_wall, dx) # [K/W]

            # find temperature drop across wall
            delta_T_w = node.Q * R_w # [K]

            # find target temperature drop
            delta_T_w_target = node.T_hw - node.T_cw # [K]

            return delta_T_w - (delta_T_w_target) # we want to minimize the difference / find the root (when the difference is zero)

        # solve for wall thickness
        return scipy.optimize.root_scalar(func_to_solve, bracket = [0.0001, 0.1], x0 = t_w_guess).root # the bracket is the range of values to search for the root, x0 is the initial guess

    t_w = solve_for_t_w(node0) # [m] - wall thickness at first node

    # add wall thickness to first node
    node0.set_t_w(t_w)

    # check result validity
    #print("t_w:", t_w, "[m]\nTemperature drop across wall:", node0.Q * R_wall(node0.r, t_w, k_wall, dx), "[K]")


    # define function to find wall height that produces target thermal resistance R_c of coolant boundary layer
    def solve_for_h(node):
        def func_to_solve(h):
            # find target thermal resistance of coolant boundary layer
            R_c_target = (node.T_cw - node.T_c) / node.Q # [K/W]

            # calculate velocity of coolant at node
            
            if node.X < 1: # determine coolant properties
                rho_c = rho_coolant_phase(node.P_c, node.X)
                Pr_c = Pr_coolant_phase(node.P_c, node.X)
                k_c = k_coolant_phase(node.P_c, node.X)
                mu_c = mu_coolant_phase(node.P_c, node.X)
            else:
                rho_c = rho_coolant_super(node.T_c, node.P_c)
                Pr_c = Pr_coolant_super(node.T_c, node.P_c)
                k_c = k_coolant_super(node.T_c, node.P_c)
                mu_c = mu_coolant_super(node.T_c, node.P_c)
            v = o_mdot / (rho_c * node.w * N_channels * h) # coolant velocity [m/s]
            
            # solve for convective heat transfer coefficient of coolant boundary layer from height h
            h_c = h_coolant_dittus_boelter(rho_c, v, Dh_rectangular(node.w, h), mu_c, Pr_c, k_c)

            # solve for thermal resistance of coolant boundary layer
            R_c = R_convective(node.r + node.t_w, h_c, dx) # [K/W]

            return R_c - R_c_target # we want to minimize the difference / find the root (when the difference is zero)

        # solve for coolant wall height
        return scipy.optimize.root_scalar(func_to_solve, bracket = [0.00001, 0.1], x0 = h_guess).root # the bracket is the range of values to search for the root, x0 is the initial guess

    # solve for channel height and velocity at node
    h = solve_for_h(node0) # [m] - coolant wall height at first node

    if node0.X < 1: # determine coolant properties
        rho_c = rho_coolant_phase(node0.P_c, node0.X)
        Pr_c = Pr_coolant_phase(node0.P_c, node0.X)
        k_c = k_coolant_phase(node0.P_c, node0.X)
        mu_c = mu_coolant_phase(node0.P_c, node0.X)
    else:
        rho_c = rho_coolant_super(node0.T_c, node0.P_c)
        Pr_c = Pr_coolant_super(node0.T_c, node0.P_c)
        k_c = k_coolant_super(node0.T_c, node0.P_c)
        mu_c = mu_coolant_super(node0.T_c, node0.P_c)


    v_c = o_mdot / (rho_c * node0.w * N_channels * h) # coolant velocity [m/s]

    # add coolant wall height and velocity to node
    node0.set_h(h)
    node0.set_v_c(v_c)
    node0.set_dv_c(0) # [m/s] - set coolant velocity change in first node to 0

    # validate result
    R_c_target = (node0.T_cw - node0.T_c) / node0.Q # [K/W]
    h_c = h_coolant_dittus_boelter(rho_c, v_c, Dh_rectangular(node0.w, h), mu_c, Pr_c, k_c)
    R_c = R_convective(node0.r + node0.t_w, h_c, dx) # [K/W]

    #print("Calculated coolant wall height in cm:", h * 100.0, "[cm]", "\nTarget thermal resistance:", R_c_target, "[K/W]", "\nCalculated thermal resistance:", R_c, "[K/W]", "\nCoolant velocity:", v_c, "[m/s]")


    # define friction factor
    f_D = 0.03002 # friction factor for turbulent flow from Moody chart

    # find pressure at next node
    def P_next(node):

        if node.X <= 1:
            rho_c = rho_coolant_phase(node.P_c, node.X)
        else:
            rho_c = rho_coolant_super(node.T_c, node.P_c)
        
        return node.P_c - dp_darcy_weisbach(f_D, rho_c, node.v_c, Dh_rectangular(node.w, node.h), node.drdx, dx) - rho_c * node.v_c * (node.dv_c) # [Pa] (currently also calculating pressure drop due to tapering in channels using velocity estimates from current channel)

    nodesc[1].set_P_c(P_next(node0)) # [Pa] - pressure at next node

    #print("Pressure at current node:", node0.P_c, "[Pa]", "\nPressure at next node:", nodesc[1].P_c, "[Pa]")



    # estimate vapor quality assuming all heat went into boiling coolant
    def H_next(node): # enthalpy at next node is enthalpy in this node + heat flow in
        return node.H + node.Q / o_mdot # [J/kg] - enthalpy at next node = enthalpy at current node + heat added to coolant / mass of coolant

    def X_next(node):
        H_n = H_next(node) # vapor quality at next node is gotten from enthalpy at next node and pressure at next node
        return X_coolant(H_n, P_next(node))

    # find temperature at next node
    def T_next(node):
        # temperature at next node is found based on vapor quality at next node
        if X_next(node) < 1:
            return T_sat(P_next(node)) # if in two phase region, temperature at next node is saturation temperature at next node pressure
        else:
            return node.T_c + node.Q / (o_mdot * cp_coolant_super(node.T_c, node.P_c)) # [K] - if in superheated state, all heat goes into raising temperature


    nodesc[1].set_H(H_next(node0)) # [J/kg] - Enthalpy at next node
    nodesc[1].set_X(X_next(node0)) # [-] - Vapor Quality at next node
    nodesc[1].set_T_c(T_next(node0)) # [K] - Temperature at next node

    #print("Temperature at current node:", node0.T_c, "[K]", "\nTemperature at next node:", nodesc[1].T_c, "[K]")
    #print("Enthalpy at current node:", node0.H, "[J/kg]", "\nEnthalpy at next node:", nodesc[1].H, "[J/kg]")
    #print("Vapor quality at current node:", node0.X, "\nVapor quality at next node:", nodesc[1].X)



    # print final result at first node
    #print("Node 0:")
    #node0.print_coolant_parameters()
    #print("\nNode 1:")
    #nodesc[1].print_coolant_parameters()



    # iterate through all nodes, starting at node 1, if error thrown, break loop and plot full chamber with all current results for all nodes
    for i, node in enumerate(nodesc[1:-1]): # iterate from N=1 to N-2 node (node before last node)
        # note: N of node is i+1 because we start at node 1 (i=0)

        # step 1 - solve for wall thickness
        try:
            t_w = solve_for_t_w(node) # [m] - wall thickness at first node
        except ValueError as error:
            print("\nValue Error at Node:", i+1)
            node.print_coolant_parameters()
            print("Error:", error)
            break
        node.set_t_w(t_w)

        # step 2 - solve for coolant wall height and velocity
        try:
            h = solve_for_h(node) # [m] - coolant wall height at first node
        except ValueError as error:
            print("\nValue Error at Node:", i+1)
            node.print_coolant_parameters()
            print("Error:", error)
            break

        if node.X < 1: # determine coolant properties
            rho_c = rho_coolant_phase(node.P_c, node.X)
            k_c = k_coolant_phase(node.P_c, node.X)
            mu_c = mu_coolant_phase(node.P_c, node.X)
            cp_c = cp_coolant_phase(node.P_c, node.X)
        else:
            rho_c = rho_coolant_super(node.T_c, node.P_c)
            k_c = k_coolant_super(node.T_c, node.P_c)
            mu_c = mu_coolant_super(node.T_c, node.P_c)
            cp_c = cp_coolant_super(node.T_c, node.P_c)

        v_c = o_mdot / (rho_c * node.w * N_channels * h) # coolant velocity [m/s]
        node.set_h(h)
        node.set_v_c(v_c)
        node.set_dv_c(v_c - nodesc[i].v_c) # [m/s] - set coolant velocity change as velocty at current node minus velocity at previous node with N = i

        # set coolant properties at node
        node.set_rho_c(rho_c) # set coolant density at first node
        node.set_k_c(k_c) # set coolant conductivity at first node
        node.set_mu_c(mu_c) # set coolant viscosity at first node
        node.set_cp_c(cp_c) # set coolant specific heat at first node
    
        # step 3 - find pressure at next node (step 3)
        nodesc[i+2].set_P_c(P_next(node))

        # step 4 - find enthalpy, quality, temperature at next node
        nodesc[i+2].set_H(H_next(node)) # [J/kg] - enthalpy at next node
        nodesc[i+2].set_X(X_next(node)) # [-] - vapor quality at next node
        nodesc[i+2].set_T_c(T_next(node)) # [K] - temperature at next node
        
        # print results
        if printNodes:
            print("\nNode", i+1, ":")
            node.print_coolant_parameters()


    # set values at final node (node N-1)
    node = nodesc[-1]
    try:
        t_w = solve_for_t_w(node) # [m] - wall thickness at first node
    except ValueError as error:
        print("\nValue Error at Node:", i+1)
        node.print_coolant_parameters()
        print("Error:", error)
    node.set_t_w(t_w)

    # step 2 - solve for coolant wall height and velocity
    try:
        h = solve_for_h(node) # [m] - coolant wall height at first node

        if node.X < 1: # determine coolant properties
            rho_c = rho_coolant_phase(node.P_c, node.X)
            k_c = k_coolant_phase(node.P_c, node.X)
            mu_c = mu_coolant_phase(node.P_c, node.X)
            cp_c = cp_coolant_phase(node.P_c, node.X)
        else:
            rho_c = rho_coolant_super(node.T_c, node.P_c)
            k_c = k_coolant_super(node.T_c, node.P_c)
            mu_c = mu_coolant_super(node.T_c, node.P_c)
            cp_c = cp_coolant_super(node.T_c, node.P_c)

        v_c = o_mdot / (rho_c * node.w * N_channels * h) # coolant velocity [m/s]
        node.set_h(h)
        node.set_v_c(v_c)
        node.set_dv_c(v_c - nodesc[-2].v_c) # [m/s] - set coolant velocity change in first node to 0

        # set properties of coolant at node
        node.set_rho_c(rho_c) # set coolant density at last
        node.set_k_c(k_c) # set coolant conductivity at last
        node.set_mu_c(mu_c) # set coolant viscosity at last
        node.set_cp_c(cp_c) # set coolant specific heat at last
    
    except ValueError as error:
        print("\nValue Error at Node:", i+1)
        node.print_coolant_parameters()
        print("Error:", error)


    # calculate and set Reynolds number at each node
    for node in nodesc:
        if node.X < 1:
            rho_c = rho_coolant_phase(node.P_c, node.X)
            mu_c = mu_coolant_phase(node.P_c, node.X)
        else:
            rho_c = rho_coolant_super(node.T_c, node.P_c)
            mu_c = mu_coolant_super(node.T_c, node.P_c)

        node.set_Re_c(reynolds(Dh_rectangular(node.w, node.h), node.v_c, rho_c, mu_c))

    # calculate mach number at each node - note holds only if gas
    M_coolant = np.zeros_like(nodesc)
    MW_coolant = 44.013e-3 # [kg/mol] - molecular weight of coolant
    for i, node in enumerate(nodesc):
        if node.X <= 1:
            rho_c = rho_coolant_phase(node.P_c, node.X)
            cp_c = cp_coolant_phase(node.P_c, node.X)
        else:
            rho_c = rho_coolant_super(node.T_c, node.P_c)
            cp_c = cp_coolant_super(node.T_c, node.P_c)

        gamma_c = (cp_c / (cp_c - 8.314 / MW_coolant)) # [-] - specific heat ratio of coolant

        M_coolant[i] = mach_number(gamma_c, node.v_c, node.T_c, MW_coolant)

    # update wall geometry lists
        
    # Some final calcualtions for plotting and saving
    
    x = [node.x for node in nodes[::-1]] # [m] - axial position at each node (in normal order)
    r = [node.r for node in nodes[::-1]] # [m] - radius at each node (in normal order)
    th = [node.t_w for node in nodesc[::-1]] # [m] - wall thickness at each node (in normal order)
    # set rest of wall thickness at chamber to constant value equal to that of cooling jacket exit wall thickness
    th = th + [nodesc[0].t_w for i in range(len(x) - len(th))]

    return nodesc[-1].T_c, nodesc[-1].P_c


''' Analyze Flow thorugh Injector '''
def runThroughInjector(C_D = 0.7, A_o = 0.0001, verbose = False, useSPI = False):
    nodef = nodesc[-1] # final node

    # define values for injector parameters
    C_D = 0.7 # [-] - discharge coefficient for coolant flow
    A_o = 0.0001 # [m^2] - area of injector

    # calculate values for injector parameters
    rho_1 = 0
    if nodef.X < 1:
        rho_1 = rho_coolant_phase(nodef.P_c, nodef.X)
    else:
        rho_1 = rho_coolant_super(nodef.T_c, nodef.P_c)
    rho_2 = PropsSI('D', 'P', def_P0, 'T', def_Ti, coolant_name) # [kg/m^3] - density of chamber nitrous estimate

    h_1 = nodef.H # [J/kg] - enthalpy of coolant before entering injector
    h_2 = PropsSI('H', 'P', def_P0, 'T', def_Ti, coolant_name) # [J/kg] - enthalpy of coolant after entering injector estimate

    # run injector calculations
    mdot_calc_HEM = 0
    mdot_calc_SPI = 0
    Del_P_SPI = 0
    
    np.seterr(all='raise')
    # HEM Calculations
    try:
        mdot_calc_HEM = C_D * A_o * np.sqrt(2 * rho_2 * (h_1 - h_2)) # [kg/s] - mass flow rate of coolant
    except:
        mdot_calc_HEM = 0
        if verbose:
            print("Warning raised: h_1 = {} less than h_2 = {}".format(h_1, h_2), "unable to retrieve HEM mass flow value.")
    np.seterr(all='warn')
    
    # SPI Calculations
    np.seterr(all='raise')
    try:
        mdot_calc_SPI = C_D * A_o * np.sqrt(2 * rho_1 * (nodef.P_c - def_P0)) # [kg/s] - mass flow rate of coolant
        Del_P_SPI = 1 / (2 * rho_1) * (o_mdot / (C_D * A_o))**2 # [Pa] - pressure drop across injector
    except:
        print("Warning raised: P_1 = {} less than P_2 = {}".format(nodef.P_c, def_P0), "unable to retrieve HEM mass flow value.")
    np.seterr(all='warn')
    
    if verbose:
        print("Enthalpy at injector:", nodes[-1].H, "[J/kg]")
        print("Enthalpy estimate inside chamber:", h_2, "[J/kg]")
        print("Difference in Enthalpy:", h_1 - h_2, "[J/kg]")
        print("Calculated mass flow rate (HEM):", mdot_calc_HEM, "[kg/s]")
        print("Calculated mass flow rate (SPI):", mdot_calc_SPI, "[kg/s]")
        print("Pressure drop across injector (SPI):", Del_P_SPI / 1e6, "[MPa]")
    
    return mdot_calc_HEM, mdot_calc_SPI, Del_P_SPI


''' Individual Plotting '''
def plot_geometry(axs):
    axs.plot(x, r, color = "k") # upper innner wall
    axs.plot(x, -np.array(r), color = "k") # lower inner wall
    axs.plot(x, np.add(r, th), color = "k") # upper outer wall
    axs.plot(x, np.subtract(-np.array(r), th), color = "k") # lower outer wall
    axs.grid()
    axs.set_xlabel("Axial Position (m)")
    axs.set_ylabel("Radius (mm)")

def plotCoolantVelocity():
    # plot velocity of coolant at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Velocity Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "darkviolet", linestyle = "--", linewidth = 0.5)
    # plot velocity
    axs2.plot([node.x for node in nodesc], [node.v_c for node in nodesc], color = "darkviolet")
    axs2.set_ylabel("Velocity [m/s]", color = "darkviolet")
    plt.savefig(solvedplotfolder + "coolantvel.png", dpi=300)
    plt.show()
    
def plotCoolantVelocityAxs(axs, ID):
    # plot velocity of coolant at each node on overlay of engine geometry
    plot_geometry(axs)
    axs.set_title("Coolant Velocity Across Chamber, ID: {}".format(ID))
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "darkviolet", linestyle = "--", linewidth = 0.5)
    # plot velocity
    axs2.plot([node.x for node in nodesc], [node.v_c for node in nodesc], color = "darkviolet")
    axs2.set_ylabel("Velocity [m/s]", color = "darkviolet")
    return axs

def plotCoolantPressure():
    # plot pressure of coolant at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Pressure Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "g", linestyle = "--", linewidth = 0.5)
    # plot pressure 
    axs2.plot([node.x for node in nodesc], [node.P_c for node in nodesc], color = "g")
    axs2.set_ylabel("Pressure [Pa]", color = "g")
    plt.savefig(solvedplotfolder + "coolantpres.png", dpi=300)
    plt.show()

def plotCoolantPressureAxs(axs, ID):
    # plot pressure of coolant at each node on overlay of engine geometry
    plot_geometry(axs)
    axs.set_title("Coolant Pressure Across Chamber, ID: {}".format(ID))
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "g", linestyle = "--", linewidth = 0.5)
    # plot pressure 
    axs2.plot([node.x for node in nodesc], [node.P_c for node in nodesc], color = "g")
    axs2.set_ylabel("Pressure [Pa]", color = "g")
    return axs   

def plotCoolantTemperature():
    # plot temperature of coolant at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Temperature Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "coral", linestyle = "--", linewidth = 0.5)
    # plot pressure 
    axs2.plot([node.x for node in nodesc], [node.T_c for node in nodesc], color = "coral")
    axs2.set_ylabel("Temperature [K]", color = "coral")
    plt.savefig(solvedplotfolder + "coolanttemp.png", dpi=300)
    plt.show()

def plotCoolantQuality():
    # plot coolant Vapor Quality at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Vapor Quality Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "r", linestyle = "--", linewidth = 0.5)
    # plot hydraulic diameter
    axs2.plot([node.x for node in nodesc], [node.X for node in nodesc], color = "r")
    axs2.set_ylabel("Vapor Quality [-]", color = "r")
    plt.show()
    
def plotCoolantQualityAxs(axs, ID):
    # plot coolant Vapor Quality at each node on overlay of engine geometry
    plot_geometry(axs)
    axs.set_title("Coolant Vapor Quality Across Chamber, ID: {}".format(ID))
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "r", linestyle = "--", linewidth = 0.5)
    # plot hydraulic diameter
    axs2.plot([node.x for node in nodesc], [node.X for node in nodesc], color = "r")
    axs2.set_ylabel("Vapor Quality [-]", color = "r")
    return axs

# Coolant Property Plotting
def plotCoolantDensity():
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Density Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "dodgerblue", linestyle = "--", linewidth = 0.5)
    # plot density
    axs2.plot([node.x for node in nodesc], [node.rho_c for node in nodesc], color = "dodgerblue")
    axs2.set_ylabel("Density [kg/m^3]", color = "dodgerblue")
    plt.show()

def plotCoolantDensityAxs(axs, ID):
    plot_geometry(axs)
    axs.set_title("Coolant Density Across Chamber, ID: {}".format(ID))
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "dodgerblue", linestyle = "--", linewidth = 0.5)
    # plot density
    axs2.plot([node.x for node in nodesc], [node.rho_c for node in nodesc], color = "dodgerblue")
    axs2.set_ylabel("Density [kg/m^3]", color = "dodgerblue")
    return axs

def plotCoolantConductivity():
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Conductivity Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "lime", linestyle = "--", linewidth = 0.5)
    # plot conductivity
    axs2.plot([node.x for node in nodesc], [node.k_c for node in nodesc], color = "lime")
    axs2.set_ylabel("Conductivity [W/m/K]", color = "lime")
    plt.show()
    
def plotCoolantConductivityAxs(axs, ID):
    plot_geometry(axs)
    axs.set_title("Coolant Conductivity Across Chamber, ID: {}".format(ID))
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "lime", linestyle = "--", linewidth = 0.5)
    # plot conductivity
    axs2.plot([node.x for node in nodesc], [node.k_c for node in nodesc], color = "lime")
    axs2.set_ylabel("Conductivity [W/m/K]", color = "lime")
    return axs

def plotCoolantViscosity():
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Viscosity Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "indigo", linestyle = "--", linewidth = 0.5)
    # plot viscosity
    axs2.plot([node.x for node in nodesc], [node.mu_c for node in nodesc], color = "indigo")
    axs2.set_ylabel("Viscosity [Pa*s]", color = "indigo")
    plt.show()
    
def plotCoolantViscosityAxs(axs, ID):
    plot_geometry(axs)
    axs.set_title("Coolant Viscosity Across Chamber, ID: {}".format(ID))
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "indigo", linestyle = "--", linewidth = 0.5)
    # plot viscosity
    axs2.plot([node.x for node in nodesc], [node.mu_c for node in nodesc], color = "indigo")
    axs2.set_ylabel("Viscosity [Pa*s]", color = "indigo")
    return axs

def plotCoolantSpecificHeat():
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Specific Heat Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "orangered", linestyle = "--", linewidth = 0.5)
    # plot specific heat
    axs2.plot([node.x for node in nodesc], [node.cp_c for node in nodesc], color = "orangered")
    axs2.set_ylabel("Specific Heat [J/kg/K]", color = "orangered")
    plt.show()
    
def plotCoolantSpecificHeatAxs(axs, ID):
    plot_geometry(axs)
    axs.set_title("Coolant Specific Heat Across Chamber, ID: {}".format(ID))
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "orangered", linestyle = "--", linewidth = 0.5)
    # plot specific heat
    axs2.plot([node.x for node in nodesc], [node.cp_c for node in nodesc], color = "orangered")
    axs2.set_ylabel("Specific Heat [J/kg/K]", color = "orangered")
    return axs

def plotCoolingJacket():
    # plot full chamber geometry visualization
    fig, axs = plt.subplots(figsize = (12, 8))
    fig.set_facecolor('white')
    # plot geometry walls by filling in between inner wall and outer wall
    axs.fill_between(x, r, np.add(r, th), color = "k")
    axs.fill_between(x, -np.array(r), np.subtract(-np.array(r), th), color = "k")
    # plot coolant walls at each node on overlay by filling in between outer wall and coolant wall
    axs.fill_between([node.x for node in nodesc], [node.r + node.t_w for node in nodesc], [node.r + node.t_w + node.h for node in nodesc], color = "b")
    axs.fill_between([node.x for node in nodesc], [-node.r - node.t_w for node in nodesc], [-node.r - node.t_w - node.h for node in nodesc], color = "b")
    axs.set_title("Cooling Jacket Geometry Across Chamber")
    axs.set_aspect('equal')
    axs.grid()
    axs2 = axs.twinx()
    axs2.grid(color = "brown", linestyle = "--", linewidth = 0.5)
    # plot coolant wall height
    axs2.plot([node.x for node in nodesc], [node.h * 1000 for node in nodesc], color = "steelblue", label="Coolant Wall Height")
    axs2.plot([node.x for node in nodesc], [node.t_w * 1000 for node in nodesc], color = "slategray", label="Chamber Wall Thickness")
    axs2.set_ylabel("Property Height [mm]", color = "brown")
    plt.legend()
    plt.savefig(solvedplotfolder + "fullchamber.png", dpi=300)
    plt.show()

def plotCoolingJacketAxs(axs, ID):
    # plot full chamber geometry visualization
    # plot geometry walls by filling in between inner wall and outer wall
    axs.fill_between(x, r, np.add(r, th), color = "k")
    axs.fill_between(x, -np.array(r), np.subtract(-np.array(r), th), color = "k")
    # plot coolant walls at each node on overlay by filling in between outer wall and coolant wall
    axs.fill_between([node.x for node in nodesc], [node.r + node.t_w for node in nodesc], [node.r + node.t_w + node.h for node in nodesc], color = "b")
    axs.fill_between([node.x for node in nodesc], [-node.r - node.t_w for node in nodesc], [-node.r - node.t_w - node.h for node in nodesc], color = "b")
    axs.set_title("Cooling Jacket Geometry Across Chamber, ID: {}".format(ID))
    axs.set_aspect('equal')
    axs.grid()
    axs2 = axs.twinx()
    axs2.grid(color = "brown", linestyle = "--", linewidth = 0.5)
    # plot coolant wall height
    axs2.plot([node.x for node in nodesc], [node.h * 1000 for node in nodesc], color = "steelblue", label="Coolant Wall Height")
    axs2.plot([node.x for node in nodesc], [node.t_w * 1000 for node in nodesc], color = "slategray", label="Chamber Wall Thickness")
    axs2.set_ylabel("Property Height [mm]", color = "brown")
    plt.legend()
    return axs

    


''''Subplot Plotting '''
subplotfig = None
axes = None
subplotcounter = 0
def plotSubplot(N):
    global subplotfig, axes, subplotcounter
    if subplotfig is None:
        subplotfig, axes = plt.subplots(N, 1, figsize = (12, 8 * N))
    
    subplotcounter += 1
    return axes[subplotcounter - 1]
    # plot 

def plotSubplot6():
    global subplotfig, axes, subplotcounter
    if subplotfig is None:
        subplotfig, axes = plt.subplots(3, 2, figsize = (18, 6 * 3))
    
    subplotcounter += 1
    if subplotcounter < 3:
        return axes[0, subplotcounter - 1]
    elif subplotcounter < 5:
        return axes[1, subplotcounter - 3]
    else:
        return axes[2, subplotcounter - 5]

def showSubplot(title = None):
    if title is not None:
        subplotfig.suptitle(title)
    plt.tight_layout()
    plt.subplots_adjust(hspace = 0.2)
    plt.show()
    
def saveSubplot(filename, title = None, show = False):
    if title is not None:
        subplotfig.title(title)
    plt.tight_layout()
    plt.subplots_adjust(hspace = 0.2)
    plt.savefig(analysisplotfolder + filename, dpi=300)
    if show:
        plt.show()
        
def resetSubplot():
    global subplotfig, axes, subplotcounter
    subplotfig = None
    axes = None
    subplotcounter = 0
    
''' Bundled Plotting '''
def plotAllCoolantProperties():
    # plot pressure of coolant at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Pressure Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "g", linestyle = "--", linewidth = 0.5)
    # plot pressure 
    axs2.plot([node.x for node in nodesc], [node.P_c for node in nodesc], color = "g")
    axs2.set_ylabel("Pressure [Pa]", color = "g")
    plt.savefig(solvedplotfolder + "coolantpres.png", dpi=300)
    plt.show()

    # plot temperature of coolant at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Temperature Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "coral", linestyle = "--", linewidth = 0.5)
    # plot pressure 
    axs2.plot([node.x for node in nodesc], [node.T_c for node in nodesc], color = "coral")
    axs2.set_ylabel("Temperature [K]", color = "coral")
    plt.savefig(solvedplotfolder + "coolanttemp.png", dpi=300)
    plt.show()

    # plot velocity of coolant at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Velocity Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "darkviolet", linestyle = "--", linewidth = 0.5)
    # plot velocity
    axs2.plot([node.x for node in nodesc], [node.v_c for node in nodesc], color = "darkviolet")
    axs2.set_ylabel("Velocity [m/s]", color = "darkviolet")
    plt.savefig(solvedplotfolder + "coolantvel.png", dpi=300)
    plt.show()

    # plot Reynolds number of coolant at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Reynolds Number Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "fuchsia", linestyle = "--", linewidth = 0.5)
    # plot Reynolds number
    axs2.plot([node.x for node in nodesc], [node.Re_c for node in nodesc], color = "fuchsia")
    axs2.set_ylabel("Reynolds Number [-]", color = "fuchsia")
    plt.savefig(solvedplotfolder + "coolantreynolds.png", dpi=300)
    plt.show()


def plotAllCoolingJacketResults():
    # plot wall thickness of chamber at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Wall Thickness Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "brown", linestyle = "--", linewidth = 0.5)
    # plot wall thickness 
    axs2.plot([node.x for node in nodesc], [node.t_w * 1000 for node in nodesc], color = "brown")
    axs2.set_ylabel("Wall Thickness [mm]", color = "brown")
    plt.savefig(solvedplotfolder + "wallthickness.png", dpi=300)
    plt.show()

    # plot coolant wall height at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Wall Height Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "b", linestyle = "--", linewidth = 0.5)
    # plot coolant wall height 
    axs2.plot([node.x for node in nodesc], [node.h * 1000 for node in nodesc], color = "b")
    axs2.set_ylabel("Coolant Wall Height [mm]", color = "b")
    plt.savefig(solvedplotfolder + "coolantwallheight.png", dpi=300)
    plt.show()

    # plot full chamber geometry visualization
    fig, axs = plt.subplots(figsize = (12, 8))
    fig.set_facecolor('white')
    # plot geometry walls by filling in between inner wall and outer wall
    axs.fill_between(x, r, np.add(r, th), color = "k")
    axs.fill_between(x, -np.array(r), np.subtract(-np.array(r), th), color = "k")
    # plot coolant walls at each node on overlay by filling in between outer wall and coolant wall
    axs.fill_between([node.x for node in nodesc], [node.r + node.t_w for node in nodesc], [node.r + node.t_w + node.h for node in nodesc], color = "b")
    axs.fill_between([node.x for node in nodesc], [-node.r - node.t_w for node in nodesc], [-node.r - node.t_w - node.h for node in nodesc], color = "b")
    axs.set_title("Cooling Jacket Geometry Across Chamber")
    axs.set_aspect('equal')
    axs.grid()
    axs2 = axs.twinx()
    axs2.grid(color = "brown", linestyle = "--", linewidth = 0.5)
    # plot coolant wall height
    axs2.plot([node.x for node in nodesc], [node.h * 1000 for node in nodesc], color = "steelblue", label="Coolant Wall Height")
    axs2.plot([node.x for node in nodesc], [node.t_w * 1000 for node in nodesc], color = "slategray", label="Chamber Wall Thickness")
    axs2.set_ylabel("Property Height [mm]", color = "brown")
    plt.legend()
    plt.savefig(solvedplotfolder + "fullchamber.png", dpi=300)
    plt.show()

    # plot coolant channel area at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Channel Area Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "aqua", linestyle = "--", linewidth = 0.5)
    # plot coolant channel area
    axs2.plot([node.x for node in nodesc], [node.A_c * 1000000 for node in nodesc], color = "aqua")
    axs2.set_ylabel("Coolant Channel Area [mm^2]", color = "aqua")
    plt.show()

    # plot coolant hydraulic diameter at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Channel Hydraulic Diameter Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "darkorange", linestyle = "--", linewidth = 0.5)
    # plot hydraulic diameter
    axs2.plot([node.x for node in nodesc], [node.Dh_c * 1000 for node in nodesc], color = "darkorange")
    axs2.set_ylabel("Coolant Hydraulic Diameter [mm]", color = "darkorange")
    plt.show()


def plotAdditionalPropertyResults():
    # plot coolant Vapor Quality at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Vapor Quality Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "r", linestyle = "--", linewidth = 0.5)
    # plot hydraulic diameter
    axs2.plot([node.x for node in nodesc], [node.X for node in nodesc], color = "r")
    axs2.set_ylabel("Vapor Quality [-]", color = "r")

    # plot coolant Mach Number (assuming gas) at each node on overlay of engine geometry
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    plot_geometry(axs)
    axs.set_title("Coolant Mach Number Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "y", linestyle = "--", linewidth = 0.5)
    # plot Mach number
    axs2.plot([node.x for node in nodesc], M_coolant, color = "y")
    axs2.set_ylabel("Mach Number [-]", color = "y")
    plt.show()

''' Other Miscaleaneous Plotting '''
def plotHotWallTargets():
    # plot extrapolated hot wall temperature T_hw(x) at each node
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    axs.plot(x, r, color = "k") # upper innner wall
    axs.plot(x, -np.array(r), color = "k") # lower inner wall
    axs.plot(x, np.add(r, th), color = "k") # upper outer wall
    axs.plot(x, np.subtract(-np.array(r), th), color = "k") # lower outer wall
    axs.grid()
    axs.set_xlabel("Axial Position (m)")
    axs.set_ylabel("Radius (m)")
    axs.set_title("Hot Wall Temperature Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "r", linestyle = "--", linewidth = 0.5)
    axs2.plot(x_n, T_hw_n, color = "r")
    axs2.set_ylabel("Temperature [K]", color = "r")
    plt.show()

def plotColdWallTargets():
    # plot extrapolated cold wall temperature T_cw(x) at each node
    fig, axs = plt.subplots(figsize = (12, 6))
    fig.set_facecolor('white')
    axs.plot(x, r, color = "k") # upper innner wall
    axs.plot(x, -np.array(r), color = "k") # lower inner wall
    axs.plot(x, np.add(r, th), color = "k") # upper outer wall
    axs.plot(x, np.subtract(-np.array(r), th), color = "k") # lower outer wall
    axs.grid()
    axs.set_xlabel("Axial Position (m)")
    axs.set_ylabel("Radius (m)")
    axs.set_title("Cold Wall Temperature Across Chamber")
    axs.set_aspect('equal')
    axs2 = axs.twinx()
    axs2.grid(color = "b", linestyle = "--", linewidth = 0.5)
    axs2.plot(x_n, T_cw_n, color = "b")
    axs2.set_ylabel("Temperature [K]", color = "b")
    plt.show()
    
    
''' Exporting Results '''
def exportResultsToCSV():
    # use pandas t export Mach number, temperature, and pressure data to csv
    nodedf = pd.DataFrame({"x [m]": x_n, "r [m]": r_n, "M [-]": [node.M for node in nodes], "T [K]": [node.T for node in nodes], "P [Pa]": [node.P for node in nodes]})

    # export to csv
    nodedf.to_csv("enginefiles/nodedata.csv")

def exportResultsToExcel():
    # get lists of node position coolant pressure, temperature, velocity, channel height, and wall thickness
    x = [node.x for node in nodesc]
    r = [node.r for node in nodesc]
    P_c = [node.P_c for node in nodesc]
    T_c = [node.T_c for node in nodesc]
    v_c = [node.v_c for node in nodesc]
    h = [node.h for node in nodesc]
    w = [node.w for node in nodesc]
    t_w = [node.t_w for node in nodesc]

    # make dataframe and save it to excel file
    df = pd.DataFrame({"Axial Position of Node [m]": x, "Radius at Node [m]": r, "Coolant Pressure [Pa]": P_c, "Coolant Temperature [K]": T_c, "Coolant Velocity [m/s]": v_c, "Channel Height [m]": h, "Channel Width [m]": w, "Wall Thickness [m]": t_w})
    df.to_excel("enginefiles/solverdata.xlsx", index = True)