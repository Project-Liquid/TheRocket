function [mdot,T1] = InjectorFlow(tank, P1, Cd, A, i)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% PARAMETERS

%rho = py.CoolProp.CoolProp.PropsSI("rho", "T", tank.T, "Q", 0, tank.prop);
h0 = tank.hLiq(i);
P0 = tank.P(i); %Take tnak pressure from index i
PVap = tank.PVap(i);
rho0 = tank.rhoLiq(i);
s0 = tank.sLiq(i);
s1 = s0; %Assuming isentropic expansion across injector
%h1 = py.CoolProp.CoolProp.PropsSI("H", "P", P1, "S", s1, tank.prop);
%rho1 = py.CoolProp.CoolProp.PropsSI("D", "P", P1, "S", s1, tank.prop);
T1 = py.CoolProp.CoolProp.PropsSI("T", "P", P1, "S", s1, tank.prop);



%% SPI: Single Phase Incompressible Model

mdot_SPI = Cd * A * sqrt(2*rho0*(P0 - P1));

%% HEM: Homogenous Equilibrium Model

%mdot_HEM = Cd * A * rho1 * sqrt(2*(h0 - h1));

%% Dyer: Non-Homogenous Non-Equilibrium Model

K = sqrt((P0 - P1) / (PVap - P1));
K = 9999;
mdot = (K/(K+1))*mdot_SPI; %+ (1/(1+K))*mdot_HEM; %Not true, just want to make sure things are working
%mdot = mdot_HEM;
mdot = mdot_SPI;
end

