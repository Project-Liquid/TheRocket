from modules import solver as sol

import matplotlib.pyplot as plt
import sys

''' Define Parameters '''
def_P0 = 2.758e+6 # [Pa] - first chamber (stagnation) pressure guess 
def_Ti = 300      # [K] - first chamber inlet temperature guess
def_P_ci = 5.5e+6 # [Pa] - coolant inlet pressure
def_X_ci = 0.4 # [-] - coolant inlet quality
def_x_off = 0.024 # [m] - cooling jacket offset
def_T_hw = 1200 # [K] - target hot wall temperature (should be below 1373 K)
def_T_cw = 600 # [K] - target cold wall temperature (should be below 623 K)
def_N_channels = 16 # [-] - number of cooling channels
def_t_rib = 0.005 # [m] thickness of ribs 
def_k_wall = 30 # [W/mK] - thermal conductivity of wall
def_N = 100 # number of nodes

# define values for injector parameters
def_C_D = 0.7 # [-] - discharge coefficient for coolant flow
def_A_o = 3.9e-5 # [m^2] - area of injector

it_P_chamber = list() # [Pa] - chamber pressure
it_T_chamber = list() # [K] - chamber temperature


''' Setup '''
sol.setupDirectories()

''' Run Solver 

sol.defineParameters(def_P0, 
                     def_Ti, 
                     def_P_ci, 
                     def_X_ci, 
                     def_x_off, 
                     def_T_hw, 
                     def_T_cw, 
                     def_N_channels, 
                     def_t_rib, 
                     def_k_wall, 
                     def_N)

T_coolant_f, P_coolant_f = sol.run()
'''
#print("Coolant Final Temperature after First Cycle: ", T_coolant_f, " K", "\nCoolant Final Pressure after First Cycle: ", P_coolant_f, " Pa")


''' Run Results Through Injector 
_, _, Del_P_SPI = sol.runThroughInjector(def_C_D, def_A_o, verbose = True)

sol.plotCoolantDensity()
#sol.plotCoolantPressure()
#sol.plotCoolantQuality()
#sol.plotCoolingJacket()
'''

'''Run Special Plots '''
#iterate over coolant quality from 0 to 0.5 in steps of 0.1
N_range = 10
def_ID = "Geometry with varying T_hw:"
for i in range(N_range):
    def_T_hw = 1200 - i*20
    sol.defineParameters(def_P0, 
                         def_Ti, 
                         def_P_ci, 
                         def_X_ci, 
                         def_x_off, 
                         def_T_hw, 
                         def_T_cw, 
                         def_N_channels, 
                         def_t_rib, 
                         def_k_wall, 
                         def_N)
    T_coolant_f, P_coolant_f = sol.run()
    #print("Coolant Final Temperature after Cycle ", i+2, ": ", T_coolant_f, " K", "\nCoolant Final Pressure after Cycle ", i+2, ": ", P_coolant_f, " Pa")
    sol.plotCoolingJacketAxs(axs = sol.plotSubplot(N_range), ID = def_ID + " {}".format(def_T_hw))

#sol.showSubplot("Velocity with varying X_ci")
sol.saveSubplot("jacketvsT_hw.png")
sol.resetSubplot()
#iterate over coolant quality from 0 to 0.5 in steps of 0.1

N_range = 10
def_ID = "Geometry with varying T_cw:"
def_T_hw = 1100
for i in range(N_range):
    def_T_cw = 600 - i*10
    sol.defineParameters(def_P0, 
                         def_Ti, 
                         def_P_ci, 
                         def_X_ci, 
                         def_x_off, 
                         def_T_hw, 
                         def_T_cw, 
                         def_N_channels, 
                         def_t_rib, 
                         def_k_wall, 
                         def_N)
    T_coolant_f, P_coolant_f = sol.run()
    #print("Coolant Final Temperature after Cycle ", i+2, ": ", T_coolant_f, " K", "\nCoolant Final Pressure after Cycle ", i+2, ": ", P_coolant_f, " Pa")
    sol.plotCoolingJacketAxs(axs = sol.plotSubplot(N_range), ID = def_ID + " {}".format(def_T_cw))

#sol.showSubplot("Velocity with varying X_ci")
sol.saveSubplot("jacketvsT_cw.png")
sol.resetSubplot()
#iterate over coolant quality from 0 to 0.5 in steps of 0.1

'''
N_range = 10
def_ID = "Pressure with varying X_ci:"
for i in range(N_range):
    def_X_ci = i*0.2
    sol.defineParameters(def_P0, 
                         def_Ti, 
                         def_P_ci, 
                         def_X_ci, 
                         def_x_off, 
                         def_T_hw, 
                         def_T_cw, 
                         def_N_channels, 
                         def_t_rib, 
                         def_k_wall, 
                         def_N)
    T_coolant_f, P_coolant_f = sol.run()
    #print("Coolant Final Temperature after Cycle ", i+2, ": ", T_coolant_f, " K", "\nCoolant Final Pressure after Cycle ", i+2, ": ", P_coolant_f, " Pa")
    sol.plotCoolantPressureAxs(axs = sol.plotSubplot(N_range), ID = def_ID + " {}".format(def_X_ci))

#sol.showSubplot("Velocity with varying X_ci")
sol.saveSubplot("pressurevsX_ci.png")
'''
sys.exit()
''' Iterate ''' 
N_iter = 20
for i in range(N_iter):
    sol.defineParameters(P_coolant_f - Del_P_SPI, 
                         T_coolant_f, 
                         def_P_ci, 
                         def_X_ci, 
                         def_x_off, 
                         def_T_hw, 
                         def_T_cw, 
                         def_N_channels, 
                         def_t_rib, 
                         def_k_wall, 
                         def_N)

    print("Iteration ", i+2, " of ", N_iter)
    try:
        T_coolant_f, P_coolant_f = sol.run()
    except:
        print("Iteration failed at cycle: {}".format(i+2))
        sys.exit()

    _, _, Del_P_SPI = sol.runThroughInjector(def_C_D, def_A_o)
    print("Coolant Final Temperature after Cycle ", i+2, ": ", T_coolant_f, " K", "\nCoolant Final Pressure after Cycle ", i+2, ": ", P_coolant_f, " Pa", "\nPressure in Chamber after Injector (SPI estimate): ",  P_coolant_f - Del_P_SPI, " Pa")

    it_P_chamber.append(P_coolant_f - Del_P_SPI)
    it_T_chamber.append(T_coolant_f)
    
''' Plot Results'''
#sol.plotCoolantVelocity()
sol.plotCoolantPressure()
sol.plotCoolantQuality()
#sol.plotCoolingJacket()

''' Additional Plots '''
# plot chamber pressure and temperature over iterations
fig = plt.figure(figsize  = (10, 5)) # plot pressure
plt.plot(range(N_iter), it_P_chamber, 'g', label = 'Chamber Pressure')
plt.xlabel('Iterations')
plt.ylabel('Pressure [Pa]')
plt.legend()
plt.show()

fig = plt.figure(figsize  = (10, 5)) # plot temperature
plt.plot(range(N_iter), it_T_chamber, 'r', label = 'Chamber Temperature')
plt.xlabel('Iterations')
plt.ylabel('Temperature [K]')
plt.legend()
plt.show()

