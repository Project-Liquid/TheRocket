from rocketcea.cea_obj import CEA_Obj 
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams.update({'axes.linewidth' : 1,
                     'ytick.major.width' : 1,
                     'ytick.minor.width' : 1,
                     'xtick.major.width' : 1,
                     'xtick.minor.width' : 1,
                     'xtick.labelsize': 10, 
                     'ytick.labelsize': 10,
                     'axes.labelsize': 12,
                     'font.family': 'Serif',
                      'figure.figsize': (6.4, 4.8)
                    })

# Nitrous oxide (N2O) as our oxidizer and Ethane (C2H6) as our fuel 
cr = 5.547795942 #Contraction ratio/Subsonic area ratio (chamber area over throat area ratio)
supAr = 4.126430558 # Supersonic expansion ratio (nozzle exit area over throat area ratio)
Pinj = 400 # [psia] Pressure at the injector face 
Pamb = 14.696 # [psia] Ambient pressure 
MR = 5.0 # O/F ratio Mass ratio of oxidizer over fuel 
At = 1.73422562 # Throat Area [in^2] (0.001118853 [m^2])
Ae = 7.15615971 # Exit Area [in^2] (0.004616868 [m^2])


ispObj = CEA_Obj(oxName = 'N2O', fuelName = 'C2H6', fac_CR = cr)

# Fitting Pinj over Pcomb (Pressure in chamber) to fit from the contraction ratio 
PinjOverPcomb = 1.0 + 0.54 / cr**2.2

# Using CEA to solve for more accurate PinjOverPcomb using previous value 
PinjOverPcomb = ispObj.get_Pinj_over_Pcomb(Pc = (Pinj / PinjOverPcomb), MR=MR)

# Find Pc (Chamber Pressure) with new PinjOverPcomb
Pc = (Pinj / PinjOverPcomb) # Solved Pc [psia]

# Find Cf (Thrust Coefficient)
Cf = ispObj.get_PambCf(Pamb=Pamb, Pc=Pc, MR=MR, eps=supAr)
# Assumes optimum expansion and solves for thrust using F=Cf*At*Pc: At = Throat Area & Pc = Chamber pressure
Thrust = Cf[1] * At * Pc # [lbf]
print(f'Pressure at injector face: {Pinj:.2f}')
print(f'Thrust {Thrust:.2f}')


# Finds ambient  Specific Impulse 
# IspSL, mode = ispObj.estimate_Ambient_Isp(Pc = Pc, MR=MR, eps=supAr, Pamb=Pamb)

# Gets full report of analysis, 
# short_output = 0: Full report 
# short_output = 0: Shortened report 
s = ispObj.get_full_cea_output(Pc = (Pinj / PinjOverPcomb), MR=MR, eps=supAr, short_output=1, pc_units='psia')
print( s )


# Initializes lists 
Cfamb = []
Cf = []
Pressures = []
Pcombs = []
Pexits = []
Pinj = 400 # [psia] Pressure at injector face

# Lowers injector face pressure and finds Cf for that specific pressure 
while Pinj > 0:
    # Finds Chamber Pressure as seen above 
    PinjOverPcomb = 1.0 + 0.54 / cr**2.2
    PinjOverPcomb = ispObj.get_Pinj_over_Pcomb(Pc = (Pinj / PinjOverPcomb), MR=MR)
    Pc = (Pinj / PinjOverPcomb)

    # Adds found chamber pressure to list
    Pcombs.append(Pc)

    # Finds CFcea: Cf where ambient is equal to exit pressure and CFamb: Cf adjusted for ambient (not sure how though)
    CFcea, CFamb, mode = ispObj.get_PambCf(Pamb=Pamb, Pc=Pc, MR=MR, eps=supAr)
    
    # Adds these to the corresponding list alongside pressure at the injector 
    Cfamb.append(CFamb)
    Cf.append(CFcea)
    Pressures.append(Pinj)
    
    # Finds PcOverPe (Chamber pressure over Exit Pressure)
    PcOverPe = ispObj.get_PcOvPe(Pc=Pc, MR=MR, eps=supAr, frozen=0, frozenAtThroat=0)
    
    # Finds exit pressure 
    Pexits.append(Pc / PcOverPe)
    
    # Lowers injector pressure 
    Pinj -= 1

# Converts lists to arrays 
Cfamb = np.array(Cfamb)
Cf = np.array(Cf)
Pcombs = np.array(Pcombs)
Pexits = np.array(Pexits)

# Adjusts Cf with equation 26 in 'Project_Liquid_Rocket_Equation_Sheet.pdf' as of (02/18/25)
CfAdj = Cf + ((Pexits - Pamb) / Pcombs) * (Ae/At)

# Calculates thrust 
Thrust = Cfamb * At * Pcombs
ThrustCea = Cf * At * Pcombs
ThrustAdj = CfAdj * At * Pcombs

# Plots thrust curve against pressure 
fig, axs = plt.subplots()
fig.set_facecolor('white')
axs.plot(Pressures, Thrust, label = 'CF ambient', color='tomato') 
axs.plot(Pressures, ThrustAdj, label = 'CF adjusted', color='mediumseagreen')
axs.plot(Pressures, ThrustCea, label = 'CF Cea', color='cornflowerblue')
axs.legend() 
axs.grid()
axs.set_xlabel(r'$P_{inj} [psia]$')
axs.set_ylabel("Thrust [lbf]")
axs.set_title(r'Thrust Curve, $N_2O$ $C_2H_6$')
axs.set_aspect('auto')
axs.invert_xaxis()
plt.savefig("ThrustCurve.png", dpi=300)