clear
clc
%% General Comments

%This script is intended to demonstrate the functionality of the
%transientModel.m function for some arbitrary injector, nozzle, and tank
%conditions

%% Dependencies

%vapInject.m
%transientModel.m
%tankBlowdown.m
%H2O_emmissions_coeff.m
%engine_param_univ.m
%emittance.m
%CoolProp.m
%CO2_emissions_coeff.m
%chamber_pressure.m
%cantera_conversion.m
%BLC_Calc.m
%bisect_method.m

%tankProperties.m
%rocketProperties.m
%nozzleGeometry.m
%injectorGeometry.m

%Cantera Installation
%CoolProp Installation

addpath("Supporting Functions") %This adds the supporting functions from the folder
addpath("Supporting Classes") %This adds the supporting object classes from the folder

%% Conversion Factors
%Conversion factors from imperial to metric
in2m = 0.0254;
psi2pa = 6894.76;
lbf2N = 4.44822;

%% Initial Condition Demonstrator
%=========================================================================
%Normally an input argument, but hardcoded for reference
pBlow = 0*psi2pa; %Constant blowdown pressure [Pa]

pSCOx = 800*psi2pa; %One-time initial ox supercharge pressure [Pa]
VSCpcOx = 0.49; %Ox ullage volume [%]

pSCFu = 0*psi2pa; %One-time initial fuel supercharge pressure [Pa]
VSCpcFu = 0.05; %Fuel ullage volume [%]


pAmbient = 101325; %Ambient pressure [Pa]

T0Ox = 266.5; %Initial ox temperature [K]
T0Fu = 299.8; %Initial fuel temperature[K]

timeStep = 0.025; %Time Step [s]

damping = 0.3; %Damping coefficient of iterative solvers

cstarEff = 0.85; %Assumed cstar efficiency
%=========================================================================

%% Injector Geometry

%Zone 1 injector geometry - 16 count
AZone1 = pi/4*0.00097^2; %[m^2] Areas Using Chamber geometry and Zones specified in Slide 18 of Delta CDR
CdZone1 = 0.85;
Count1 = 16;

%Zone 2 injector geometry - 24 count
AZone2 = in2m^2*pi/4*0.018^2; %[m^2]
CdZone2 = 0.7;
Count2 = 24;

%Zone 3 injector geometry - 52 count
AZone3 = in2m^2*pi/4*0.018^2; %[m^2]
CdZone3 = 0.7;
Count3 = 52;

% Venturi Geometry
CdVenturi = 0.995;
AVenturiExit = pi*0.00275^2;
CountVenturi = 1;

%Initialize injector as an object.
injector = injectorGeometry([CdVenturi, AVenturiExit, CountVenturi],...
                            [CdZone1, AZone1, Count1; CdZone2, AZone2, Count2],...
                            [CdZone3, AZone3, Count3]);

%% Engine Geometry

Dc = 3.834*in2m; %Combustion chamber diameter [m]
Dt = 1.335*in2m; %Throat diameter [m]
Lc = 4.25*in2m; %Combustion chamber length [m]
At = pi/4*(1.335*in2m)^2; %Throat Area [m^2]
r1 = 0; %Converging curvature radius [m]
r2 = 0; %Diverging curvature radius [m]
theta1 = 30; %Contraction angle [deg]
theta2 = 15; %Expansion angle [deg]
Cont_R = (Dc/Dt)^2; % Chamber contraction ratio [Unitless]
Exp_R = 5; %Expansion ratio [Unitless]

%Initialize nozzle as an object
nozzle = nozzleGeometry(Exp_R, Cont_R, Dt, r1, theta1, r2, theta2, Lc);


%% Tank Properties

%Static Fire 2023 Tanks (Static Fire 1-Gal)

%Ox tank
oxidizer = 'N2O';
VOx = 0.00378541; %[m^3]
VSC_Ox = VOx*VSCpcOx; %[m^3] Vapor/supercharge volume above props in ox tank

oxTank = tankProperties(pBlow, T0Ox, VOx, oxidizer);
oxTank.VSC = VSC_Ox;

%Fuel tank
fuel = 'Ethane';
VFu = 0.00378541; %[m^3]
VSC_Fu = VFu*VSCpcFu; %[m^3] Vapor/supercharge volume above props in fuel tank

fuelTank = tankProperties(pBlow, T0Fu, VFu, fuel);
fuelTank.VSC = VSC_Fu;

%Place in array of objects for input
tankList = [oxTank, fuelTank];

pSC = [pSCOx, pSCFu];

%% Run Transient Model

[modelData, timevec] = transientModel(tankList, injector, nozzle, pBlow, pSC, pAmbient, pAmbient*2, cstarEff, damping, timeStep);

%% Print Outputs

%%Print returns

%Return burn time and flameout characteristics
fprintf('Burn time: %.2f seconds\n', timevec(end))
if modelData.QOx(end)-0.5 > modelData.QFuel(end)
    fprintf('Fuel-rich flameout\n')
elseif modelData.QFuel(end)-0.5 > modelData.QOx(end)
    fprintf('Ox-rich flameout\n')
else
    fprintf('Approx. simultaneous flameout\n')
end

%Return maximum temperature
fprintf('\nMax Wall Temperature: %.1f K or %.1f F\n', max(modelData.Tmax), max(modelData.Tmax)*9/5-460)
fprintf('Max Throat Temperature: %.1f K or %.1f F\n', max(modelData.Tthroat), max(modelData.Tmax)*9/5-460)

%Return maximum performance
fprintf('\nMax Chamber Temp: %.1f K or %.1f R\n', max(modelData.Tc), max(modelData.Tc)*9/5);
fprintf('Max Chamber Pressure: %.2f bar or %.1f psi\n', max(modelData.pc)/1E5, max(modelData.pc)/psi2pa)
fprintf('Max Thrust Coefficient: %.2f\n', max(modelData.Cf))
fprintf('Max Isp: %.1f s \n', max(modelData.Cf.*modelData.cstar./9.81))
fprintf('Max Thrust: %.1f N or %.1f lbf\n', max(modelData.thrust), max(modelData.thrust)/lbf2N)

%% Plot Outputs

%Tank Conditions
figure(1)
clf
subplot(2, 4, 1)
plot(timevec(1:end-1), modelData.T0Ox(1:end-1), 'g')
title('Ox Tank Temperature vs Time')
ylabel('Ox Temperature, T, [K]')
xlabel('Time, t, [s]')

subplot(2, 4, 2)
plot(timevec(1:end-1), modelData.p0Ox(1:end-1)./psi2pa, 'g')
title('Ox Tank Pressure vs Time')
ylabel('Ox Pressre, p, [psi]')
xlabel('Time, t, [s]')

subplot(2, 4, 3)
plot(timevec, modelData.mOx, 'g')
title('Ox Tank Mass vs Time')
ylabel('Ox Mass Remaining, m, [kg]')
xlabel('Time, t, [s]')
ylim([0, max(modelData.mOx)*1.1])

subplot(2, 4, 4)
plot(timevec, 100*(modelData.QOx), 'g')
title('Ox Tank Quality vs Time')
ylabel('Ox Quality, Q, [%]')
xlabel('Time, t, [s]')
ylim([0, 100])

subplot(2, 4, 5)
plot(timevec(1:end-1), modelData.T0Fuel(1:end-1), 'r')
title('Fuel Tank Temperature vs Time')
ylabel('Fuel Temperature, T, [K]')
xlabel('Time, t, [s]')

subplot(2, 4, 6)
plot(timevec(1:end-1), modelData.p0Fuel(1:end-1)./psi2pa, 'r')
title('Fuel Tank Pressure vs Time')
ylabel('Fuel Pressre, p, [psi]')
xlabel('Time, t, [s]')

subplot(2, 4, 7)
plot(timevec, modelData.mFuel, 'r')
title('Fuel Tank Mass vs Time')
ylabel('Fuel Mass Remaining, m, [kg]')
xlabel('Time, t, [s]')
ylim([0, max(modelData.mFuel)*1.1])

subplot(2, 4, 8)
plot(timevec, 100*(modelData.QFuel), 'r')
title('Fuel Tank Quality vs Time')
ylabel('Fuel Quality, Q, [%]')
xlabel('Time, t, [s]')
ylim([0, 100])

%Chamber Conditions
figure(2)
clf

subplot(2, 2, 1)
plot(timevec, modelData.pc./psi2pa)
title('Chamber Pressure vs Time')
ylabel('Chamber Pressre, pc, [psi]')
xlabel('Time, t, [s]')

subplot(2, 2, 2)
plot(timevec, modelData.Tc)
title('Chamber Temperature vs Time')
ylabel('Chamber Temperature, Tc, [K]')
xlabel('Time, t, [s]')

subplot(2, 2, 3)
hold on
plot(timevec, modelData.mdot, 'k')
plot(timevec, modelData.mdot_ox, 'g')
plot(timevec, modelData.mdot_fuel, 'r')
title('Mass Flow Rates vs Time')
ylabel('Mass Flow Rate, mdot, [kg/s]')
xlabel('Time, t, [s]')
legend('Total', 'Ox', 'Fuel')

subplot(2, 2, 4)
hold on
plot(timevec, modelData.mdot_ox./modelData.mdot_fuel, 'r')
plot(timevec, modelData.mdot_ox./(modelData.mdot_fuel + modelData.mdot_BLC), 'k')
title('Mix Ratio vs Time')
ylabel('Mixture Ratio, r, O/F')
xlabel('Time, t, [s]')
legend('Combustion O/F', 'Total O/F')

%Performance Parameters
figure(3)
clf

subplot(2, 2, 1)
plot(timevec, modelData.thrust./lbf2N)
title('Thrust vs Time')
ylabel('Thrust, F, [lbf]')
xlabel('Time, t, [s]')

subplot(2, 2, 2)
plot(timevec, modelData.cstar.*modelData.Cf./9.81)
title('Specific Impulse vs Time')
ylabel('Specific Impulse, Isp, [s]')
xlabel('Time, t, [s]')

subplot(2, 2, 3)
plot(timevec, modelData.cstar)
title('Characteristic Velocity vs Time')
ylabel('Characteristic Velocity, c*, [m/s]')
xlabel('Time, t, [s]')

subplot(2, 2, 4)
plot(timevec, modelData.Cf)
title('Thrust Coefficient vs Time')
ylabel('Thrust Coefficient, Cf')
xlabel('Time, t, [s]')

% Thermal Parameters (Woohoo - Austin)
figure(4)
clf

subplot(2, 2, 1)
plot(timevec,modelData.Tthroat)
title('Max Steady-State Wall Temperature')
ylabel('Max Wall Temperature, Tthroat [K]')
xlabel('Time, t, [s]')

subplot(2, 2, 2)
plot(modelData.L_array(1, :)./in2m, modelData.h_array(1, :))
title('Burnout Film Coefficient Contour')
ylabel('Film Coefficient, h, [W/m^2-K]')
xlabel('Axial Length, L, [in]')

subplot(2, 2, 3)
plot(modelData.L_array(1, :)./in2m, modelData.Tw_array(1, :))
title('Initial Temperature Contour')
ylabel('Wall Temperature, Tw, [K]')
xlabel('Axial Length, L, [in]')

subplot(2, 2, 4)
plot(modelData.L_array(end, :)./in2m, modelData.Tw_array(end, :))
title('Burnout Temperature Contour')
ylabel('Wall Temperature, Tw, [K]')
xlabel('Axial Length, L, [in]')