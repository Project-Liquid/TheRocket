clear
clc
clearvars

tic
%% Target

% 1. To output thrust profile of vehicle


%NOTES
%Intermediate calculations are mostly in metric

i = 1; %Start at i = 1, since this is the default starter index for matlab
preallocation = 10000; %Preallocation of memory for any arrays
timeVec = zeros(1,preallocation);

%Estimate Tank Sizing Ratio
OFtarget = 5;
oxDensityStart = 786;
fuelDensityStart = 339;
heightRatioOF = OFtarget * fuelDensityStart / oxDensityStart;


%All in inches below
heightNoseCone = 24;
heightMainBay = 22.5;
heightAvBay = 5;
heightDrogueBay = 10;
heightCoupler3 = 5;
heightOx = 25;
heightFuel = heightOx / heightRatioOF; %originally 11.5
heightCoupler2 = 5;
heightCoupler1 = 5;
heightEngine = 10; 
OD = 5;

heightVehicle = (heightEngine + heightCoupler1 + heightOx + heightCoupler2 + heightFuel + heightCoupler3 + heightDrogueBay + heightAvBay + heightMainBay + heightNoseCone)/12;

%TESTING EXCEL READ
%type ./MasterParameters.xlsx
T = readtable("../../MasterParameters.xlsx");

oxDensityStart = excelReader("oxDensityStart", T);

%% Conversion Factors
%Conversion factors from imperial to metric
in2m = 0.0254;
psi2pa = 6894.76;
lbf2N = 4.44822;
lbm2kg = 0.45359237; %for converting pounds to kg

%% %% GEOMETRIC POSITIONING OF ELEMENTS
oxPosition = heightEngine + heightCoupler1 + (heightOx/2);
fuelPosition = heightEngine + heightCoupler1 + heightOx + heightCoupler2 + (heightFuel/2);

oxPositionSI = oxPosition * in2m;
fuelPositionSI = fuelPosition * in2m;


%% Key User Inputs

... Tank temperature, tank volume, tank ullage volume
    
%% TANKS

%Tank General: Inputs
tankThick = 0.25; %in
temperature = 293; %Kelvin. Tempearture of each prop tank at t=0


%Tank General: Intermediate Calculations
OD_SI = OD * in2m; %SI
ID_SI = (OD - (2*tankThick)) * in2m; %SI
radiusSI = OD_SI / 2;
tankAreaInternal = (pi/4) * ID_SI^2; %SI
heightOx = heightOx * in2m; %SI
heightFuel = heightFuel * in2m; %SI

%OX: Inputs
oxidizer = 'N2O';
ullageOx = 0.15; %Between 0 to 1. Indicates volume fraction of ox
volumeOx = tankAreaInternal * heightOx;
%massOx = 8; %lbs of prop. Includes liquid and gas phase

%OX: Intermediate Calculations
%massOxSI = massOx * lbm2kg;

%Fuel: Inputs
fuel = 'Ethane';
ullageFuel = 0.15; %Between 0 to 1. Indicates volume fraction of ox
volumeFuel = tankAreaInternal * heightFuel;
%massFuel = 2; %lbs of prop. Includes liquid and gas phase

%Fuel: Intermediate Calculations
%massFuelSI = massFuel * lbm2kg;

%TANK OBJECTS
oxTank = TheTanks(oxidizer, ullageOx, temperature, tankAreaInternal, i, OD_SI, heightOx, volumeOx);
fuelTank = TheTanks(fuel, ullageFuel, temperature, tankAreaInternal, i, OD_SI, heightFuel, volumeFuel);

%% ENGINE PARAMETERS

diamThroat = 1.65; %in... original 1.48
diamThroat = diamThroat * in2m; %Converting to metric
cstarEff = 0.85;
Pamb = 101000;

engine = EngineParameters(diamThroat, i);
heightEngineSI = heightEngine * in2m;
nozzle_radiusSI = OD_SI / 2;


%% INJECTOR PARAMETERS

CdOx = 0.75;
CdFuel = 0.75;

elementDiamOx = 0.2637;%0.3259; %Inches, originally 0.44
elementDiamOx = elementDiamOx * in2m; %Metric

elementDiamFuel = 0.1797; %Inches, originally 0.17
elementDiamFuel = elementDiamFuel * in2m; %Metric

elementCountOx = 1;
elementCountFuel = 1;

injectorOx = InjectorParameters(CdOx, elementDiamOx, elementCountOx);
injectorFuel = InjectorParameters(CdFuel, elementDiamFuel, elementCountFuel);


%% Changing State of the tank

maxRelErr = 0.01;
Pcguess = 0.99 * min(oxTank.P(i), fuelTank.P(i)); %Setting min Pcguess for first iteration
%interim = 0.99 * min(oxTank.P(i), fuelTank.P(i)); %Setting min Pcguess for first iteration

time = 0;
iter = 1;
itermax = 100;
propTankQualThreshold = 0.99;
damping = 0.05;
timeDelta = 0.025;

relErrCompile = zeros(2,400);

%Reseting to 0 until I get clever
i = 0;
%ITERATE UNTIL ONE TANK IS COMPLETELY VAPOR
while fuelTank.qualUnderThreshold && oxTank.qualUnderThreshold
    
    relErr = 1; %Initialize to start mass flow iterative solver
    iter = 1;
    
    %ITERATE ON i FOR EVERY STEP EXCEPT THE FIRST
    
        i = i + 1;
    
    pressureDebugger = zeros(11, 100);
    
    %% TIME ITERATION 
    timeVec(i) = time;
    
    %Return progress through time step
    fprintf('Time Step: \t%0.2f s\n', timeVec(i))
    
    %% Determine SS flow rate at time step using InjectorFlow and EngineFlow
    while relErr > maxRelErr
         
        %% FLOW RATE: INJECTOR.... use Pcguess
        [mdot_f, Tc_f] = InjectorFlow(fuelTank, Pcguess, injectorFuel.Cd, injectorFuel.cumulativeArea, i);
        [mdot_ox, Tc_ox] = InjectorFlow(oxTank, Pcguess, injectorOx.Cd, injectorOx.cumulativeArea, i);
        m_dot = mdot_f + mdot_ox;
        OF = mdot_ox / mdot_f;
        Tinj = (mdot_ox * Tc_ox + mdot_f * Tc_f) / (mdot_ox + mdot_f); %Weighted average of input temperature in chamber before combustion
        
        %% FLOW RATE: ENGINE.... use m_dot from Pcguess to calculated Pcnew through throat
        [Pcnew, Tc, R, gamma] = EngineFlow(oxTank, mdot_ox, fuelTank, mdot_f, m_dot, Tinj, engine, Pcguess);
        
        %% COMPARE FLOW RATE GUESSES
        relErr = (Pcguess - Pcnew) / (Pcguess);
       
        
        
        %% SET MAXIMUM CYCLES
        
        if iter < itermax
            iter = iter + 1;
        else
            break
        end
        
        %% DEBUGGER FOR PRESSURE COMPARISON
        pressureDebugger(1, iter - 1) = Pcguess;
        pressureDebugger(2, iter - 1) = Pcnew;
        pressureDebugger(3, iter - 1) = mdot_f;
        pressureDebugger(4, iter - 1) = mdot_ox;
        pressureDebugger(5, iter - 1) = Tc;
        pressureDebugger(6, iter - 1) = m_dot;
        pressureDebugger(7, iter - 1) = OF;
        pressureDebugger(8, iter - 1) = gamma;
        pressureDebugger(9, iter - 1) = R;
        pressureDebugger(10, iter - 1) = relErr;
        pressureDebugger(11, iter - 1) = engine.areaThroat;
        
        relErrCompile(i,iter-1) = relErr;
        
        %% DAMPEN PRESSURE FOR NEXT ITERATION
        Pcguess = Pcguess - (damping * (Pcguess - Pcnew));
        
    end
    
    %% ENGINE STORAGE
    % Some code storing Pcnew into engine parameter
    engine.Pc(i) = Pcnew;
    engine.m_dot(i) = m_dot;
    engine.m_dotox(i) = mdot_ox;
    engine.m_dotfuel(i) = mdot_f;
    engine.Tc(i) = Tc;
    engine.gamma(i) = gamma;
    engine.R(i) = R;
    
    %% THRUST PARAMETERS
    %dummy values for now
    
    %Assume exit pressure is ambient for simplicity
    Pe = Pamb;
    
    % Calculates characteristic velocity
    engine.charV(i) = (engine.Pc(i) * engine.areaThroat / engine.m_dot(i)) * cstarEff;
    
    % Calculates coefficient of thrust
    engine.Cf(i) = sqrt((2*engine.gamma(i)^2/(engine.gamma(i)-1))*((2/(engine.gamma(i)+1))^((engine.gamma(i)+1)/(engine.gamma(i)-1)))*(1-(Pe/engine.Pc(i))^((engine.gamma(i)-1)/engine.gamma(i))));
    
    % Uses previous values and mass flow rate to calculate total thrust
    engine.thrust(i) = engine.Cf(i) * engine.m_dot(i) * engine.charV(i);
    
    % Uses relation between thrust and mass flow rate to calculate the
    % specific impulse
    engine.Isp(i) = engine.thrust(i)/(9.81*engine.m_dot(i));
   
    
    %% TANK ITERATION
    deltamOx = mdot_ox * timeDelta;
    deltamFuel = mdot_f * timeDelta;
    
    %Store new tank state in i+1 as long as new qual doesn't exceed
    %threshold
    TankChangingState(fuelTank, deltamFuel, propTankQualThreshold, i);
    TankChangingState(oxTank, deltamOx, propTankQualThreshold, i);
    
    
    %% TIME ITERATION
    time = time + timeDelta;
    Pcguess = engine.Pc(i); %Set new chamber pressure guess for next round
    
    
    %clc
end

%% DELETING ZERO VALUES OF MATRICES
%starting with just some of the import ones for CSV files
timeVec(i+1:end) = [];
engine.thrust(i+1:end) = [];
oxTank.mLiq(i+1:end) = [];
oxTank.mVap(i+1:end) = [];
fuelTank.mLiq(i+1:end) = [];
fuelTank.mVap(i+1:end) = [];

%Prepping CSV files
enginePerformance = [timeVec; engine.thrust];
oxTankMassLiq = [timeVec; oxTank.mLiq];
oxTankMassVap = [timeVec; oxTank.mVap];
fuelTankMassLiq = [timeVec; fuelTank.mLiq];
fuelTankMassVap = [timeVec; fuelTank.mVap];

%Transpose for time - data format
enginePerformance = transpose(enginePerformance);
oxTankMassLiq = transpose(oxTankMassLiq);
oxTankMassVap = transpose(oxTankMassVap);
fuelTankMassLiq = transpose(fuelTankMassLiq);
fuelTankMassVap = transpose(fuelTankMassVap);

%% EXPORTING TO CSV FILES

csvwrite('enginePerformance', enginePerformance)
csvwrite('oxTankMassLiq', oxTankMassLiq) 
csvwrite('oxTankMassVap', oxTankMassVap)
csvwrite('fuelTankMassLiq', fuelTankMassLiq) 
csvwrite('fuelTankMassVap', fuelTankMassVap)

toc

%% PLOTS



%% Plot Outputs


%Tank Conditions
figure(1)
clf

subplot(2, 5, 1)
plot(timeVec(1:i), oxTank.T(1:i), 'g')
title('Ox Tank Temperature vs Time')
ylabel('Ox Temperature, T, [K]')
xlabel('Time, t, [s]')

subplot(2, 5, 2)
plot(timeVec(1:i), oxTank.P(1:i)./psi2pa, 'g')
title('Ox Tank Pressure vs Time')
ylabel('Ox Pressre, p, [psi]')
xlabel('Time, t, [s]')

subplot(2, 5, 3)
plot(timeVec(1:i), oxTank.m(1:i), 'g')
title('Ox Tank Mass vs Time')
ylabel('Ox Mass Remaining, m, [kg]')
xlabel('Time, t, [s]')
%ylim([0, max(modelData.mOx)*1.1])

subplot(2, 5, 4)
plot(timeVec(1:i), 100*(oxTank.x(1:i)), 'g')
title('Ox Tank Quality vs Time')
ylabel('Ox Quality, Q, [%]')
xlabel('Time, t, [s]')
%ylim([0, 100])

subplot(2,5,5)
hold on
plot(timeVec(1:i), oxTank.rhoLiq(1:i), '--g', 'LineWidth', 2.5)
plot(timeVec(1:i), oxTank.rhoVap(1:i), ':g', 'LineWidth', 2.5)
ylabel('Ox Tank Density vs Time')
xlabel('Time, t, [s]')
legend('Liquid Phase', 'Vapor Phase')



subplot(2, 5, 6)
plot(timeVec(1:i), fuelTank.T(1:i), 'r')
title('Fuel Tank Temperature vs Time')
ylabel('Fuel Temperature, T, [K]')
xlabel('Time, t, [s]')

subplot(2, 5, 7)
plot(timeVec(1:i), fuelTank.P(1:i)./psi2pa, 'r')
title('Fuel Tank Pressure vs Time')
ylabel('Fuel Pressre, p, [psi]')
xlabel('Time, t, [s]')

subplot(2, 5, 8)
plot(timeVec(1:i), fuelTank.m(1:i), 'r')
title('Fuel Tank Mass vs Time')
ylabel('Fuel Mass Remaining, m, [kg]')
xlabel('Time, t, [s]')
%ylim([0, max(modelData.mFuel)*1.1])

subplot(2, 5, 9)
plot(timeVec(1:i), 100*(fuelTank.x(1:i)), 'r')
title('Fuel Tank Quality vs Time')
ylabel('Fuel Quality, Q, [%]')
xlabel('Time, t, [s]')
%ylim([0, 100])

subplot(2,5,10)
hold on
plot(timeVec(1:i), fuelTank.rhoLiq(1:i), '--r', 'LineWidth', 2.5)
plot(timeVec(1:i), fuelTank.rhoVap(1:i), ':r', 'LineWidth', 2.5)
ylabel('Fuel Tank Density vs Time')
xlabel('Time, t, [s]')
legend('Liquid Phase', 'Vapor Phase')


%Chamber Conditions
figure(2)
clf

subplot(2, 2, 1)
plot(timeVec(1:i), engine.Pc(1:i)./psi2pa)
title('Chamber Pressure vs Time')
ylabel('Chamber Pressre, pc, [psi]')
xlabel('Time, t, [s]')

subplot(2, 2, 2)
plot(timeVec(1:i), engine.Tc(1:i))
title('Chamber Temperature vs Time')
ylabel('Chamber Temperature, Tc, [K]')
xlabel('Time, t, [s]')

subplot(2, 2, 3)
hold on
plot(timeVec(1:i), engine.m_dot(1:i), 'k')
plot(timeVec(1:i), engine.m_dotox(1:i), 'g')
plot(timeVec(1:i), engine.m_dotfuel(1:i), 'r')
title('Mass Flow Rates vs Time')
ylabel('Mass Flow Rate, mdot, [kg/s]')
xlabel('Time, t, [s]')
legend('Total', 'Ox', 'Fuel')

subplot(2, 2, 4)
hold on
plot(timeVec(1:i), engine.m_dotox(1:i) ./ engine.m_dotfuel(1:i), 'r')
title('Mix Ratio vs Time')
ylabel('Mixture Ratio, r, O/F')
xlabel('Time, t, [s]')

%Performance Parameters
figure(3)
clf

subplot(2, 2, 1)
plot(timeVec(1:i), engine.thrust(1:i) ./ lbf2N)
title('Thrust vs Time')
ylabel('Thrust, F, [lbf]')
xlabel('Time, t, [s]')

subplot(2, 2, 2)
plot(timeVec(1:i), engine.Isp(1:i))
title('Specific Impulse vs Time')
ylabel('Specific Impulse, Isp, [s]')
xlabel('Time, t, [s]')

subplot(2, 2, 3)
plot(timeVec(1:i), engine.charV(1:i))
title('Characteristic Velocity vs Time')
ylabel('Characteristic Velocity, c*, [m/s]')
xlabel('Time, t, [s]')

subplot(2, 2, 4)
plot(timeVec(1:i), engine.Cf(1:i))
title('Thrust Coefficient vs Time')
ylabel('Thrust Coefficient, Cf')
xlabel('Time, t, [s]')





%{
% Thermal Parameters (Woohoo - Austin)
figure(4)
clf


subplot(2, 2, 1)
plot(timeVec(1:i),modelData.Tthroat)
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

%}

%%
%%What's next: First, mimic tankBlowdown to show change in tank properties
%%with given mass flow rate. Then go back and use SPI/HEM model  to solve
%%for actual expected mass flow rates
%%Reference roughyl 252 to 281 in transientModel, which calls on the two
