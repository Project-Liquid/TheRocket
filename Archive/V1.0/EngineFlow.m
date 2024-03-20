function [Pcnew, Tc, R, gamma] = EngineFlow(oxTank, mdot_ox, fuelTank, mdot_f, mdot, Tinj, engine, Pcguess)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Sets the placeholder gas to the propellant OF ratio

Comb_Gas = GRI30;

gas_combustion_vector = strcat(CanteraConversion(oxTank.prop),':',num2str(mdot_ox),',',CanteraConversion(fuelTank.prop),':',num2str(mdot_f));
set(Comb_Gas,'T',Tinj,'P',Pcguess,'MassFractions',gas_combustion_vector);

% Simulates the steady-state combustion of the propellants
equilibrate(Comb_Gas,'HP');

% Finds adiabatic flame temperature, gas constant, and ratio of
% specific heats for post-combustions gasses
Tc = temperature(Comb_Gas);
R = gasconstant()/meanMolecularWeight(Comb_Gas);
gamma = cp_mass(Comb_Gas)/cv_mass(Comb_Gas);

% Take a guess at chamber pressure based on gas properties and flow
% rates. Based on Sutton (3-24) on pg. 59. CHECK IF TRUE
a1 = gamma*R*Tc;
ratio = (2/(gamma+1))^((gamma+1)/(gamma-1));

Pcnew = mdot/(engine.areaThroat * gamma*sqrt(ratio/a1));

end

