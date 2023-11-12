function TankChangingState(tank, deltam, propTankQualThreshold, i)
   
%% Density Calculations
   
    %Change in density across finite time step
    deltaRho = -deltam/tank.V; %[kg/m^3]
    
    %New density across entire tank
    rho0 = tank.rho(i); %[kg/m^3]
    rho1 = rho0 + deltaRho; %[kg/m^3]
    
    %% Enthalpy Calculations
    
    
    %Retreive specific enthalpy
    h0 = tank.h(i); %[J/kg]
    
    h0liq = py.CoolProp.CoolProp.PropsSI("H","T", tank.T(i), "Q", 0, tank.prop);
    
    %Volume specific enthalpy across entire tank
    h_vol0 = rho0*h0; %[J/m^3]
    
    %Change in volume specific enthalpy (rho*h, from Karda)
    deltah_vol = -deltam*h0liq/tank.V; %[J/m^3]
    
    %New volumetric enthalpy
    h_vol1 = h_vol0 + deltah_vol; %[J/m^3]
    
    %New specific enthalpy from new volumetric enthalpy and new density
    h1 = h_vol1/rho1; %[J/kg]
    
    
    
    %% CoolProp Results & Other Parameters
    %UNSURE OF THE EXPLICIT SITUATIONS THIS COULD BE TRIGGERED
    
    
    try
        P1 = py.CoolProp.CoolProp.PropsSI("P", "H", h1, "D", rho1, tank.prop);
        
        T1 = py.CoolProp.CoolProp.PropsSI("T", "H", h1, "D", rho1, tank.prop);
        
        x1 = py.CoolProp.CoolProp.PropsSI("Q", "H", h1, "Dmass", rho1, tank.prop);
    catch
        P1 = py.CoolProp.CoolProp.PropsSI("P", "H", h0, "D", rho0, tank.prop);
        
        T1 = py.CoolProp.CoolProp.PropsSI("T", "H", h0, "D", rho0, tank.prop);
        
        x1 = py.CoolProp.CoolProp.PropsSI("Q", "H", h0, "Dmass", rho0, tank.prop);
      
    end
    
    
    
    %T1 = CoolProp("T", "H", h1, "D", rho1, tank.prop);
    %P1 = CoolProp("P", "H", h1, "D", rho1, tank.prop); 
    PVap1 = py.CoolProp.CoolProp.PropsSI("P", "T", T1, "Q", 1, tank.prop);
    
    %x1 = py.CoolProp.CoolProp.PropsSI("Q", "H", h1, "Dmass", rho1, tank.prop); %New quality
    
    m1 = rho1 * tank.V;
    mVap1 = x1 * m1;
    mLiq1 = m1 - mVap1;
    
    rhoLiq1 = py.CoolProp.CoolProp.PropsSI("Dmass", "T", T1, "Q", 0, tank.prop);
    rhoVap1 = py.CoolProp.CoolProp.PropsSI("Dmass", "T", T1, "Q", 1, tank.prop);
    hLiq1 = py.CoolProp.CoolProp.PropsSI("H", "T", T1, "Q", 0, tank.prop);
    sLiq1 = py.CoolProp.CoolProp.PropsSI("S", "T", T1, "Q", 0, tank.prop);
    
    ullage1 = (mVap1/rhoVap1) / tank.V;
    
    %% Saving Properties
    
    if x1 > propTankQualThreshold
        tank.qualUnderThreshold = false;
        return
    end
    
    tank.m(i+1) = m1;
    tank.mVap(i+1) = mVap1;
    tank.mLiq(i+1) = mLiq1;
    tank.P(i+1) = P1;
    tank.PVap(i+1) = PVap1;
    tank.T(i+1) = T1;
    tank.h(i+1) = h1;
    tank.hLiq(i+1) = hLiq1;
    tank.rho(i+1) = rho1;
    tank.rhoLiq(i+1) = rhoLiq1;
    tank.rhoVap(i+1) = rhoVap1;
    tank.ullage(i+1) = ullage1;
    tank.x(i+1) = x1;
    tank.sLiq(i+1) = sLiq1;
    
   
end