function [rocket,exhaust] = engine_param_univ(props,injector,engine_geometry,cstarEff,p_amb,p_guess,damp)
%engine_param_univ: Finds the performance of a rocket engine given injector
%and engine geometries and propellant choice.

%% Inputs
%   props: Array of tank.m objects for fuel, ox, and BLC. If no BLC is
%   specified, the fuel object will be counted as BLC
%   engine_geometry: Object of engine geometry parameters
%   ambient_conditions: Conditions of outside atmosphere:
%       [Ambient Temp, Ambient Pressure]
%   p_guess: Initial guess of pressure [Pa]
%   cStarEff: C* Efficiency, measure of efficiency of chamber and nozzle in
%   extracting useful energy from propellants
%   damp: Damping coefficient for Gauss-Siedel Method

%% Outputs
% rocket: rocketProperties.m object containing all the necessary
% performance and chamber properties
% exhaust: Cantera object representing the exhaust gasses

if nargin < 4
    cstarEff=0.85;
elseif nargin < 5
    ambient_conditions = [290,101325];
elseif nargin < 6
    p_guess = 101325;
elseif nargin < 7
    damp = 0;
end

%% Conversion Factors
psi2pa = 6894.76; % Conversion from psi to Pascals
in2m = 0.0254; % Conversion from inches to meters
N2lbf = 0.224809; % Conversion from Newtons to pound-force

%% Engine Geometry
% Extracts necessary engine geometries
Exp_Ratio = engine_geometry.expRatio;
At = engine_geometry.At;
Ae = engine_geometry.Ae;

%% Propellant Properties
% Stores propellant name strings
ox = props(1);
fuel = props(2);

if length(props) == 2
    BLC = fuel;
else
    BLC = props(3);
end

Fuel_prop = fuel.prop;
Ox_prop = ox.prop;

CD_f = injector.fuelCd; % Coefficient of discharge of fuel orifices
CD_ox = injector.oxCd; % Coefficient of discharge of oxidizer orifices

A_f = injector.fuelArea; % Total area of fuel injector orifices (m^3)
A_ox = injector.oxArea; % Total area of ox injector orifices (m^3)


T_f = fuel.T; % Fuel temperature (K)
T_ox = ox.T; % Oxidizer temperature (K)

p_f = fuel.p; % Fuel Feed Pressure (Pa)
p_ox = ox.p; % Oxidizer Feed Pressure (Pa)

BLC_prop = BLC.prop;
CD_BLC = injector.BLCCd;
A_BLC = injector.BLCArea;


%% Iteration Parameters
p_c = p_guess;
T_c = 290;
iter = 0;
maxiter = 100;
MaxRelErr = 1E-5;
RelErr = 1;

%If initial pressure guess is above ullage pressure, get new, lower
%guess at pressure.
if p_c > min(ox.p, fuel.p)
    p_c = 0.75*min(ox.p, fuel.p);
end

Comb_Gas = GRI30;

%% Iteration
while RelErr>MaxRelErr
    p_old = p_c;
    T_old = T_c;
    
    [mdot_f, T_1_f] = vapInject(p_f,p_c,T_f,CD_f,A_f,Fuel_prop); % Fuel mass flow rate and injection temperature (kg/s, K)
    
    [mdot_ox, T_1_ox] = vapInject(p_ox,p_c,T_ox,CD_ox,A_ox,Ox_prop); % Ox mass flow rate and injection temperature (kg/s, K)
    
    [mdot_BLC, T_1_BLC] = vapInject(BLC.p,p_c,BLC.T,CD_BLC,A_BLC,BLC_prop);
    
    
    % Calculates O/F ratio an total mass flow rate
    OF = mdot_ox/mdot_f;
    mdot = mdot_ox+mdot_f+mdot_BLC;
    % Estimates combined propellant temperature at injection
    T_inject = (mdot_ox*T_1_ox+mdot_f*T_1_f/(mdot_ox+mdot_f));
    
    % Sets the placeholder gas to the propellant OF ratio
    gas_combustion_vector = strcat(cantera_conversion(Ox_prop),':',num2str(mdot_ox),',',cantera_conversion(Fuel_prop),':',num2str(mdot_f));
    set(Comb_Gas,'T',T_inject,'P',p_c,'MassFractions',gas_combustion_vector);
    
    % Simulates the steady-state combustion of the propellants
    equilibrate(Comb_Gas,'HP');
    
    % Finds adiabatic flame temperature, gas constant, and ratio of
    % specific heats for post-combustions gasses
    T_c = temperature(Comb_Gas);
    R = gasconstant()/meanMolecularWeight(Comb_Gas);
    gamma = cp_mass(Comb_Gas)/cv_mass(Comb_Gas);
    
    % Take a guess at chamber pressure based on gas properties and flow
    % rates
    p_c = chamber_pressure(T_c,mdot_f+mdot_ox,R,gamma,At);
    
    % If pressure guess is above ullage pressure, get new, lower
    % guess at pressure.
    if p_c > min(ox.p, fuel.p)
        p_c = 0.75*min(ox.p, fuel.p);
    end
    
    p_Err = abs((p_c-p_old)/p_c);
    T_Err = abs((T_c-T_old)/T_c);
    RelErr = max([p_Err,T_Err]);
    p_c = p_c-damp*(p_c-p_old);
    T_c = T_c-damp*(T_c-T_old);
    
    iter = iter+1;
    
    if iter >= maxiter
        break
    end
end

gamma_ratio = (gamma+1)/(gamma-1); % Useful gamma ratio
% Creates an implicit function for exhaust mach number, and finds the mach
% number of exhaust at nozzle exit
M = @(m) (1/Exp_Ratio)*((2/(gamma+1))+m^2/gamma_ratio)^(gamma_ratio/2)-m;
mach = bisect_method(M,1,100,1E-8);
% Uses exit mach to calculate exit pressure
pe = p_c/((1+0.5*(gamma-1)*mach^2)^(gamma/(gamma-1)));
% Calculates characteristic velocity
char_v = (p_c*At/mdot)*cstarEff;

% Calculates coefficient of thrust
Cf = sqrt((2*gamma^2/(gamma-1))*((2/(gamma+1))^((gamma+1)/(gamma-1)))*(1-(pe/p_c)^((gamma-1)/gamma)))+(Ae/At)*(pe-p_amb)/p_c;

% Uses previous values and mass flow rate to calculate total thrust
Thrust = Cf*mdot*char_v;

% Uses relation between thrust and mass flow rate to calculate the
% specific impulse
Isp = Thrust/(9.81*mdot);

%% Output

% Outputs Cantera gas object of the exhaust
exhaust = Comb_Gas;
% Creates output object
rocket = rocketProperties();
%Perform boundary layer cooling calculations with final data
if length(props) < 3
    rocket = BLC_Calc(engine_geometry, mdot_BLC, mdot, exhaust, props(2), T_1_BLC, rocket);
else
    rocket = BLC_Calc(engine_geometry, mdot_BLC, mdot, exhaust, props(3), T_1_BLC, rocket);
end
    
    
% Assigns output properties
rocket.pc = p_c;
rocket.Tc = T_c;
rocket.gamma = gamma;
rocket.R = R;

rocket.mdot = mdot;
rocket.mdot_ox = mdot_ox;
rocket.mdot_fuel = mdot_f;
rocket.mdot_BLC = mdot_BLC;

rocket.thrust = Thrust;
rocket.ISP = Isp;
rocket.cstar = char_v;
rocket.Cf = Cf;
rocket.pe = pe;

% Prints table of values if not output is specified
if nargout == 0
    fprintf("Thrust: \t %.1f lbf\n",Thrust*N2lbf)
    fprintf("ISP: \t %.1f seconds\n",Isp)
    fprintf("Chamber Pressure: \t %.1f psi\n",p_c/psi2pa)
    fprintf("Chamber Temperature: \t %.1f Kelvin\n",T_c)
    fprintf("Mixture Ratio: \t %.1f \n",OF)
end

end

