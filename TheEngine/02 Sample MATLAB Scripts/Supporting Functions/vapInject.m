function [mdot, T1] = vapInject(p0, p1, T0, Cd, A, prop)
% [mdot, T1] = vapInject(p0, p1, T0, Cd, A, prop);
%
% mdot: Mass flow rate               [kg/s]
% T1:   Injector Outlet Temperature  [K]
%
% p0:   Upstream/Tank Pressure       [Pa]
% p1:   Downstream Pressure          [Pa]
% T0:   Manifold/Tank Pressure       [K]
% Cd:   Orifice Discharge Coeff.     [Unitless]
% A:    Orifice OUTLET Area          [m^2]
% prop: Propellant String Name       [String]
%
% Outputs the mass flow rate for one orifice with an arbitrary CdA for
% given propellant "prop", given chamber conditions and tank/manifold
% conditions. Employs the Dyer model for two-phase flow, meaning that this
% model can account for two-phase injection dynamics, though for low-vapor 
% pressure propellants these effects may be negligible.


%% Get Fluid Properties

pvap = CoolProp("P", "T", T0, "Q", 1, prop); %Upstream vapor pressure [Pa]

%Get saturated liquid density for initial density since that's what's
%coming out of the tank.
rho0 = CoolProp("D", "T", T0, "Q", 0, prop); %Upstream liquid density [kg/m^3]

%Get saturated liquid specific enthalpy too.
h0 = CoolProp("H", "T", T0, "Q", 0, prop); %Upstream enthalpy [J/kg]

%Get initial entropy of the liquid. Because flow time is short flow can be 
%assumed to be isentropic, i.e. s1 = s0 = const.
s0 = CoolProp("S", "T", T0, "Q", 0, prop); %Specific entropy [J/kg]

%Use const. entropy to find density and enthalpy at outlet, i.e. at
%downstream pressure w/ same entropy
rho1 = CoolProp("D", "S", s0, "P", p1, prop); %Outlet density [kg/m^3]
h1 = CoolProp("H", "S", s0, "P", p1, prop); %Outlet sp. enthalpy [J/kg]
T1 = CoolProp("T", "S", s0, "H", h1, prop); %Outlet temperature [K]

%% Dyer Model for Two-Phase Flow

k = sqrt((p0-p1)/(pvap-p1)); %Dyer two-phase parameter [Unitless]

if p1 > pvap
    SPIFlag = 1; %If chamber pressure exceeds vapor pressure, standard SPI applies
else
    SPIFlag = 0;
end
    
%Homogenous Equilibrium Model (Enthalpy-based)
if h0-h1 > 0
    mdot_HEM = Cd*A*rho1*sqrt(2*(h0-h1)); %Enthalpy-driven mass flow rate [kg/s]
else
    mdot_HEM = 0; %Provide an exception in case we get imaginary outputs
end

if p0-p1 > 0
    mdot_SPI = Cd*A*sqrt(2*rho0*(p0-p1)); %Pressure-driven mass flow rate [kg/s]
else
    mdot_SPI = 0; %Provide an exception in case we get imaginary outputs
end

%Combine the two models for a two-phase flow rate
if SPIFlag == 0
    mdot = (k/(1+k))*mdot_SPI + (1/(1+k))*mdot_HEM; %Actual mass flow rate [kg/s]
else
    mdot = mdot_SPI; %If chamber pressure isn't below vapor pressure, SPI model applies
end

end

