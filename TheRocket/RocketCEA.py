from rocketcea.cea_obj import CEA_Obj 

# Nitrous oxide as our oxidizer and Ethane (C2H6) as our fuel 
# Contraction Ratio = 5.5509 (ratio of chamber cross sectional area and throat area) 
# with radius of the throat being 0.741303 in and the chamber radius being 1.746547
cr = 5.55
Pinj = 400 # [psia]
ispObj = CEA_Obj(oxName = 'N2O', fuelName = 'C2H6', fac_CR = cr)

PinjOverPcomb = 1.0 + 0.54 / cr**2.2

PinjOverPcomb = ispObj.get_Pinj_over_Pcomb( Pc = (Pinj / PinjOverPcomb), MR=5.0 )
s = ispObj.get_full_cea_output( Pc = (Pinj / PinjOverPcomb), MR=5.0, eps=4.3)
print( s )