from modules import solver as sol

''' Define Parameters '''
def_P0 = 2.758e+6 # [Pa] - first chamber (stagnation) pressure guess 
def_Ti = 300      # [K] - first chamber inlet temperature guess
def_P_ci = 4e+6 # [Pa] - coolant inlet pressure
def_X_ci = 0 # [-] - coolant inlet quality
def_x_off = 0.024 # [m] - cooling jacket offset
def_T_hw = 1200 # [K] - target hot wall temperature (should be below 1373 K)
def_T_cw = 500 # [K] - target cold wall temperature (should be below 623 K)
def_N_channels = 16 # [-] - number of cooling channels
def_t_rib = 0.005 # [m] thickness of ribs 
def_k_wall = 20 # [W/mK] - thermal conductivity of wall
def_N = 100 # number of nodes

''' Run Solver '''
sol.setupDirectories()
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
print("Coolant Final Temperature after First Cycle: ", T_coolant_f, " K", "\nCoolant Final Pressure after First Cycle: ", P_coolant_f, " Pa")


''' Plot First Cycle Results '''	
sol.plotCoolantVelocity()
sol.plotCoolantQuality()

''' Iterate '''
for i in range(5):
    sol.defineParameters(P_coolant_f, 
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

    T_coolant_f, P_coolant_f = sol.run()
    print("Coolant Final Temperature after Cycle ", i+2, ": ", T_coolant_f, " K", "\nCoolant Final Pressure after Cycle ", i+2, ": ", P_coolant_f, " Pa")


sol.plotCoolantVelocity()
sol.plotCoolantQuality()
sol.plotCoolingJacket()

#sol.plotCoolantPropertyResults()
#sol.plotCoolingJacketGeometryResults()
#sol.plotAdditionalPropertyResults()

#sol.exportResultsToCSV()
#sol.exportResultsToExcel()

