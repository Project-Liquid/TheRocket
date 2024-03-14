function [tank] = tankBlowdown(tank, deltam)
%% [tank] = tankBlowdown(tank, deltam);
%
% tank:     Function outputs object of class "tank.m" with new properties
%
% tank:     Function takes input object of class "tank.m" with initial properties
% deltam:   Finite change in mass across finite time interval
%
% Outputs the change in tank properties across some indeterminate finite
% time step with a change in mass delta-m provided an initial specific
% enthalpy across the tank (liquid and vapor phases included) and the
% density throughout the tank.

%% Liquid Property Extraction

    h0liq = CoolProp("H", "T", tank.T, "Q", 0, tank.prop);

%% Density Calculations (From Karda) 
    
    %Change in density across finite time step
    deltaRho = -deltam/tank.V; %[kg/m^3]
    
    %New density across entire tank
    rho0 = tank.rho; %[kg/m^3]
    rho1 = rho0 + deltaRho; %[kg/m^3]
    
%% Enthalpy Calculations {Also from Karda)

    %Retreive specific enthalpy
    h0 = tank.h; %[J/kg]

    %Volume specific enthalpy across entire tank
    h_vol0 = rho0*h0; %[J/m^3]
    
    %Change in volume specific enthalpy (rho*h, from Karda)
    deltah_vol = -deltam*h0liq/tank.V; %[J/m^3]
    
    %New volumetric enthalpy
    h_vol1 = h_vol0 + deltah_vol; %[J/m^3]
    
    %New specific enthalpy from new volumetric enthalpy and new density
    h1 = h_vol1/rho1; %[J/kg]
    
%% Outputs from CoolProp Tables
%Include try-catch to account for wonky draining terminal dynamics
    try
        p1 = CoolProp("P", "H", h1, "D", rho1, tank.prop);
        
        T1 = CoolProp("T", "H", h1, "D", rho1, tank.prop);
    catch
        p1 = CoolProp("P", "H", h0, "D", rho0, tank.prop);
        
        T1 = CoolProp("T", "H", h0, "D", rho0, tank.prop);
    end
    
%% Save to object properties

    tank.p = p1;
    tank.T = T1;
    tank.h = h1;
    tank.rho = rho1;
    tank.m = rho1*tank.V;
end

