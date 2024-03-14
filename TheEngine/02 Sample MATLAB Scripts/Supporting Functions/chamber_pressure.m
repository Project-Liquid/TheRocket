function [pc] = chamber_pressure(Tc, mdot, R, gamma, At)
%[pc] =chamber_pressure(Tc, mdot, R, gamma)
%
%This function returns a new iteration of the chamber pressure using an
%initial guess at the pressure p_i. The new pc can be used to find a new
%mix ratio to find a new Tc, R, and gamma to iterate with.
%
%pc: New guess at chamber pressure, or p of i+1, in Pa
%
%Tc: Chamber temperature in K
%T0: Ambient/tank temperature
%R: Gas constant of the exhaust mixture in J/kg-K 
%gamma: The ratio of specific heats of the exhaust mixture

%% Throat Mdot for Stagnation/Chamber Pressure
%Using Sutton (3-24) on pg. 59 to get:
a1 = gamma*R*Tc;
ratio = (2/(gamma+1))^((gamma+1)/(gamma-1));

pc = mdot/(At*gamma*sqrt(ratio/a1));
end