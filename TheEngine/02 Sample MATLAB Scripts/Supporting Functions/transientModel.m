function [rocketData, t_array] = transientModel(tanksArray, injector, nozzle, pBlow, pSC, p3, p1guess, cstarEff, damp, timeStep)
%[rocketData, t_array] = transientModel(tanksArray, injector, nozzle, p0, p3, p1gues, cstarEff, damp, timeStep);
%
%Outputs:
%
%   rocketData:     Object of class transientData containing transient
%                   data points of relevant characteristics
%
%   t_array:        Array of time values transient data points were
%                   collected against
%
%Inputs:
%
%   tanksArray:     Object array of class(es) tankProperties, of syntax:
%                   [oxTank, fuelTank, BLCTank] or [oxTank, fuelTank]
%
%   injector:       Object of class injectorGeometry containing injector
%                   design characteristics
%
%   nozzle:         Object of class nozzleGeometry containing conical
%                   nozzle & chamber design characteristics
%
%   pBlow:          Blowdown ullage gas pressure applied        [Pa]
%                   Set to zero for self-pressurizing
%
%   pSC:            List of supercharge pressures applied       [Pa]
%
%   p3:             Ambient pressure                            [Pa]
%
%   p1guess:        Initial chamber pressure guess              [Pa]
%
%   cstarEff:       Decimal estimate of cstar-efficiency
%
%   damp:           Gauss-Seidel damping applied to iteration
%
%   timeStep:       Desired time step interval                  [s]
%
%This function utilizes the engine_param_univ.m function to generate
%transient models of some arbitrarily defined rocket (tanks, propellants,
%injector, and nozzle) until one propellant tank drains. This function can
%support any non-halogen propellant combination including the list of 
%hydrocarbons below, provided CoolProp and Cantera support them. Separate
%boundary layer coolant tanks are supported, but not required. Tank design,
%injector geometry, and any conical nozzle shape are parameterized to allow
%for the modeling and simulation of a wide variety of designs.
%
%Supported Hydrocarbons:
%
%   Input hydrocarbons are to be designated by their names as specified by
%   CoolProp supported fluids list in input arguments. cantera_conversion.m
%   will then retrieve the chemical formula expression for the hydrocarbon
%   internally.
%
%   Methane:    CH4
%   Methanol:   CH3OH
%   Ethane:     C2H6
%   Propane:    C3H8
%
%   Support for additional propellants can be installed provided that
%   these propellants (dodecane, for example) are supported by CoolProp and
%   are installed into the GRI30 object of Cantera. Following this, adding
%   the conversion from proper name to chemical formula to
%   cantera_conversion.m

%% Iteration Parameters
    %Assume tank drain level/vapor quality starts at 0
    Qmax = 0; %Vapor quality [unitless]
    i = 0; %Iteration counter
    
    %Initialize vectors
    dt = timeStep; %[s]
    t_array = 0; %[s]
    
    %% Supercharge Initialization
    
    pVapOx = CoolProp('P', 'Q', 1, 'T', tanksArray(1).T, tanksArray(1).prop); %[Pa]
    pVapFuel = CoolProp('P', 'Q', 1, 'T', tanksArray(2).T, tanksArray(2).prop); %[Pa]
    
    pressures = [pBlow, pVapOx, pVapOx];
    
    if length(tanksArray) == 3
    pVapBLC = CoolProp('P', 'Q', 1, 'T', tanksArray(3).T, tanksArray(3).prop); %[Pa]
    pressures(4) = pVapBLC;
    end
    
    if pSC(1) > max([pBlow, pVapOx]) && pSC(2) > max([pBlow, pVapFuel])
        %Apply supercharge to both if both exceed blowdown and vapor
        %pressure
        if length(pSC) == 1
            pSC = pSC*ones(1, 3);
        end
        
        %Ox tank
        
        %Calculate the change in density due to application of supercharge
        %pressure using the formula rho2/rho1 = V1/V2. Inverse ratios and all
        %that.

        %Supercharge tank
        VoidOx = tanksArray(1).VSC; %[m^3] Assumed initial supercharge volume
        tanksArray(1) = tankProperties(pSC(1), tanksArray(1).T, tanksArray(1).V, tanksArray(1).prop);
        rhoVapOx = CoolProp('D', 'T', tanksArray(1).T, 'Q', 1, tanksArray(1).prop);
        rhoLiqOx = CoolProp('D', 'T', tanksArray(1).T, 'Q', 0, tanksArray(1).prop);
        mOx = rhoVapOx*VoidOx + rhoLiqOx*(tanksArray(1).V - VoidOx);
        
        tanksArray(1).rho = mOx./tanksArray(1).V; %Set new propellant density
        tanksArray(1).h = CoolProp('H', 'T', tanksArray(1).T, 'D', tanksArray(1).rho, tanksArray(1).prop);

        %Get number of moles of blowdown gas from ideal gas law (R = 8.31J/mol-K)
        %& partial pressure, where partial pressure is total pressure minus vapor
        %pressure
        
        nSCOx = (pSC(1) - pVapOx)*VoidOx/(8.31*tanksArray(1).T); %[mol]

        %Fuel tank

        %Supercharge tank
        VoidFuel = tanksArray(2).VSC; %[m^3] Assumed initial supercharge volume
        tanksArray(2) = tankProperties(pSC(2), tanksArray(2).T, tanksArray(2).V, tanksArray(2).prop);
        rhoVapFuel = CoolProp('D', 'T', tanksArray(2).T, 'Q', 1, tanksArray(2).prop);
        rhoLiqFuel = CoolProp('D', 'T', tanksArray(2).T, 'Q', 0, tanksArray(2).prop);
        mFuel = rhoVapFuel*VoidFuel + rhoLiqFuel*(tanksArray(2).V - VoidFuel);
        
        tanksArray(2).rho = mFuel./tanksArray(2).V; %Set new propellant density
        tanksArray(2).h = CoolProp('H', 'T', tanksArray(2).T, 'D', tanksArray(2).rho, tanksArray(2).prop);

        %Get number of moles of blowdown gas from ideal gas law (R = 8.31J/mol-K)
        %& partial pressure, where partial pressure is total pressure minus vapor
        %pressure

        nSCFuel = (pSC(2) - pVapFuel)*VoidFuel/(8.31*tanksArray(2).T); %[mol]

    elseif pSC(1) > max([pBlow, pVapOx])
        %Otherwise, only apply supercharge to ox tank if you can
        %Supercharge tank
        VoidOx = tanksArray(1).VSC; %[m^3] Assumed initial supercharge volume
        tanksArray(1) = tankProperties(pSC(1), tanksArray(1).T, tanksArray(1).V, tanksArray(1).prop);
        rhoVapOx = CoolProp('D', 'T', tanksArray(1).T, 'Q', 1, tanksArray(1).prop);
        rhoLiqOx = CoolProp('D', 'T', tanksArray(1).T, 'Q', 0, tanksArray(1).prop);
        mOx = rhoVapOx*VoidOx + rhoLiqOx*(tanksArray(1).V - VoidOx);
        
        tanksArray(1).rho = mOx./tanksArray(1).V; %Set new propellant density
        tanksArray(1).h = CoolProp('H', 'T', tanksArray(1).T, 'D', tanksArray(1).rho, tanksArray(1).prop);

        %Get number of moles of blowdown gas from ideal gas law (R = 8.31J/mol-K)
        %& partial pressure, where partial pressure is total pressure minus vapor
        %pressure
        
        nSCOx = (pSC(1) - pVapOx)*VoidOx/(8.31*tanksArray(1).T); %[mol]
        
        %Same for fuel, just forget the supercharge
        VoidFuel = tanksArray(2).VSC; %[m^3] Assumed initial supercharge volume
        tanksArray(2) = tankProperties(pBlow, tanksArray(2).T, tanksArray(2).V, tanksArray(2).prop);
        rhoVapFu = CoolProp('D', 'T', tanksArray(2).T, 'Q', 1, tanksArray(2).prop);
        rhoLiqFu = CoolProp('D', 'T', tanksArray(2).T, 'Q', 0, tanksArray(2).prop);
        mFuel = rhoVapFu*VoidFuel + rhoLiqFu*(tanksArray(2).V - VoidFuel);
        
        tanksArray(2).rho = mFuel./tanksArray(2).V; %Set new propellant density
        tanksArray(2).h = CoolProp('H', 'T', tanksArray(2).T, 'D', tanksArray(2).rho, tanksArray(2).prop);
            
    elseif pSC(2) > max([pBlow, pVapFuel])
            
        %Supercharge fuel tank
        VoidFuel = tanksArray(2).VSC; %[m^3] Assumed initial supercharge volume
        tanksArray(2) = tankProperties(pBlow, tanksArray(2).T, tanksArray(2).V, tanksArray(2).prop);
        rhoVapFu = CoolProp('D', 'T', tanksArray(2).T, 'Q', 1, tanksArray(2).prop);
        rhoLiqFu = CoolProp('D', 'T', tanksArray(2).T, 'Q', 0, tanksArray(2).prop);
        mFuel = rhoVapFu*VoidFuel + rhoLiqFu*(tanksArray(2).V - VoidFuel);

        tanksArray(2).rho = mFuel./tanksArray(2).V; %Set new propellant density
        tanksArray(2).h = CoolProp('H', 'T', tanksArray(2).T, 'D', tanksArray(2).rho, tanksArray(2).prop);

        %Unsupercharged ox initialization
        VoidOx = tanksArray(1).VSC; %[m^3] Assumed initial supercharge volume
        tanksArray(1) = tankProperties(pBlow, tanksArray(1).T, tanksArray(1).V, tanksArray(1).prop);
        rhoVapOx = CoolProp('D', 'T', tanksArray(1).T, 'Q', 1, tanksArray(1).prop);
        rhoLiqOx = CoolProp('D', 'T', tanksArray(1).T, 'Q', 0, tanksArray(1).prop);
        mOx = rhoVapOx*VoidOx + rhoLiqOx*(tanksArray(1).V - VoidOx);

        tanksArray(1).rho = mOx./tanksArray(1).V; %Set new propellant density
        tanksArray(1).h = CoolProp('H', 'T', tanksArray(1).T, 'D', tanksArray(1).rho, tanksArray(1).prop);

    else
        %Otherwise, just allow for plain vapor pressure ullage volume
        
            VoidOx = tanksArray(1).VSC; %[m^3] Assumed initial supercharge volume
            tanksArray(1) = tankProperties(pBlow, tanksArray(1).T, tanksArray(1).V, tanksArray(1).prop);
            rhoVapOx = CoolProp('D', 'T', tanksArray(1).T, 'Q', 1, tanksArray(1).prop);
            rhoLiqOx = CoolProp('D', 'T', tanksArray(1).T, 'Q', 0, tanksArray(1).prop);
            mOx = rhoVapOx*VoidOx + rhoLiqOx*(tanksArray(1).V - VoidOx);
            
            tanksArray(1).rho = mOx./tanksArray(1).V; %Set new propellant density
            tanksArray(1).h = CoolProp('H', 'T', tanksArray(1).T, 'D', tanksArray(1).rho, tanksArray(1).prop);
            
            
            VoidFuel = tanksArray(2).VSC; %[m^3] Assumed initial supercharge volume
            tanksArray(2) = tankProperties(pBlow, tanksArray(2).T, tanksArray(2).V, tanksArray(2).prop);
            rhoVapFu = CoolProp('D', 'T', tanksArray(2).T, 'Q', 1, tanksArray(2).prop);
            rhoLiqFu = CoolProp('D', 'T', tanksArray(2).T, 'Q', 0, tanksArray(2).prop);
            mFuel = rhoVapFu*VoidFuel + rhoLiqFu*(tanksArray(2).V - VoidFuel);
            
            tanksArray(2).rho = mFuel./tanksArray(2).V; %Set new propellant density
            tanksArray(2).h = CoolProp('H', 'T', tanksArray(2).T, 'D', tanksArray(2).rho, tanksArray(2).prop);
            
        if length(tanksArray) == 3
            %Supercharge tank
            VoidBLC = tanksArray(3).VSC; %[m^3] Assumed initial supercharge volume
            tanksArray(3) = tankProperties(pBlow, tanksArray(3).T, tanksArray(3).V, tanksArray(3).prop);
            rhoVapBLC = CoolProp('D', 'T', tanksArray(3).T, 'Q', 1, tanksArray(3).prop);
            rhoLiqBLC = CoolProp('D', 'T', tanksArray(3).T, 'Q', 0, tanksArray(3).prop);
            mBLC = rhoVapBLC*VoidBLC + rhoLiqBLC*(tanksArray(3).V - VoidBLC);

            tanksArray(3).rho = mBLC./tanksArray(3).V; %Set new propellant density
            tanksArray(3).h = CoolProp('H', 'T', tanksArray(3).T, 'D', tanksArray(3).rho, tanksArray(3).prop);
        end
    end
    
    %Coolant tank, if applicable
    if length(tanksArray) == 3
        if pSC(3) > max([pBlow, pBLC])
            %Supercharge tank
            VoidBLC = tanksArray(3).VSC; %[m^3] Assumed initial supercharge volume
            tanksArray(3) = tankProperties(pSC(3), tanksArray(3).T, tanksArray(3).V, tanksArray(3).prop);
            rhoVapBLC = CoolProp('D', 'T', tanksArray(3).T, 'Q', 1, tanksArray(3).prop);
            rhoLiqBLC = CoolProp('D', 'T', tanksArray(3).T, 'Q', 0, tanksArray(3).prop);
            mBLC = rhoVapBLC*VoidBLC + rhoLiqBLC*(tanksArray(3).V - VoidBLC);

            tanksArray(3).rho = mBLC./tanksArray(3).V; %Set new propellant density
            tanksArray(3).h = CoolProp('H', 'T', tanksArray(3).T, 'D', tanksArray(3).rho, tanksArray(3).prop);

            %Get number of moles of blowdown gas from ideal gas law (R = 8.31J/mol-K)
            %& partial pressure, where partial pressure is total pressure minus vapor
            %pressure

            nSCBLC = (pSC(3) - pVapBLC)*VoidBLC/(8.31*tanksArray(3).T); %[mol]
        else
            tanksArray(3) = tankProperties(pBlow, tanksArray(3).T, tanksArray(3).V, tanksArray(3).prop);
            rhoVapBLC = CoolProp('D', 'T', tanksArray(3).T, 'Q', 1, tanksArray(3).prop);
            rhoLiqBLC = CoolProp('D', 'T', tanksArray(3).T, 'Q', 0, tanksArray(3).prop);
            mBLC = rhoVapBLC*VoidBLC + rhoLiqBLC*(tanksArray(3).V - VoidBLC);

            tanksArray(3).rho = mBLC./tanksArray(3).V; %Set new propellant density
            tanksArray(3).h = CoolProp('H', 'T', tanksArray(3).T, 'D', tanksArray(3).rho, tanksArray(3).prop);
        end
    end
        

%% Iteration
    
%While the vapor pressure is less than 99%, or while the tank is more than
%1 percent full of liquid...
while Qmax < 0.99
    %% Iteration & Data Collection
    
    i = i + 1; %Iteration count
    
    %Time of this step is latest time plus time step
    t_array(i) = t_array(end) + dt;
    
    %Return progress through time step
    fprintf('Time Step: \t%0.2f s\n', t_array(end))
    
    %Model the engine at this time step
    [rocket,exhaust] = engine_param_univ(tanksArray, injector, nozzle, cstarEff, p3, p1guess, damp);
    
    %Get mass flow rates from rocket object
    mdot_ox   = rocket.mdot_ox;
    mdot_fuel = rocket.mdot_fuel;
    mdot_BLC  = rocket.mdot_BLC;
    
    %% Tank Dynamics
    
    %Get finite change in mass, delta-m, across this time step
    if length(tanksArray) == 2
        %If only two tanks, BLC draws from fuel tank
        deltam_ox   = mdot_ox*dt; %Oxidizer delta-m 
        deltam_fuel = (mdot_fuel + mdot_BLC)*dt; %Fuel delta-m
        
        %Apply blowdown model to tank objects
        tanksArray(1) = tankBlowdown(tanksArray(1), deltam_ox);
        tanksArray(2) = tankBlowdown(tanksArray(2), deltam_fuel);
        
        %% Supercharge Gas Dynamics
        
        %If supercharge is applicable, get supercharged pressure
        if pSC(1) > max([pBlow, pVapOx])
            %Ox Tank
            try %Empty try-catch to account for funny terminal dynamics
            pVapOx = CoolProp('P', 'Q', 1, 'T', tanksArray(1).T, tanksArray(1).prop); %[Pa]
            rhoLiqOx = CoolProp('D', 'Q', 0, 'T', tanksArray(1).T, tanksArray(1).prop); %[kg/m^3]
            end %If CoolProp throws an error, just leave these properties be
            
            VoidOx = VoidOx + deltam_ox/rhoLiqOx; %[m^3]
            
            %Sum partial pressures for ullage pressure
            pGasOx = nSCOx*(8.31)*tanksArray(1).T/VoidOx;
            tanksArray(1).p = pVapOx + pGasOx; %[Pa]
        end
        
        if pSC(2) > max([pBlow, pVapOx])
            %Fuel Tank
            try %Empty try-catch to account for funny terminal dynamics
            pVapFuel = CoolProp('P', 'Q', 1, 'T', tanksArray(2).T, tanksArray(2).prop); %[Pa]
            rhoLiqFuel = CoolProp('D', 'Q', 0, 'T', tanksArray(2).T, tanksArray(2).prop); %[kg/m^3]
            end %If CoolProp throws an error, just leave these properties be
            
            VoidFuel = VoidFuel + deltam_fuel/rhoLiqFuel; %[m^3]
            
            %Sum partial pressures for ullage pressure
            pGasFuel = nSCFuel*(8.31)*tanksArray(2).T/VoidFuel;
            tanksArray(2).p = pVapFuel + pGasFuel; %[Pa]
        end
        
        %If blowdown gas is applied, set pressure to blowdown
        if pBlow > tanksArray(1).p
            tanksArray(1).p = pBlow;
        end
        if pBlow > tanksArray(2).p
            tanksArray(2).p = pBlow;
        end
        
    else
        %Otherwise, assume separate BLC feed
        deltam_ox   = mdot_ox*dt;   %Oxidizer delta-m
        deltam_fuel = mdot_fuel*dt; %Fuel delta-m
        deltam_BLC  = mdot_BLC*dt;  %BLC delta-m
        
        %Apply blowdown model to tank objects
        tanksArray(1) = tankBlowdown(tanksArray(1), deltam_ox);
        tanksArray(2) = tankBlowdown(tanksArray(2), deltam_fuel);
        tanksArray(3) = tankBlowdown(tanksArray(3), deltam_BLC);
        
        %If supercharge is applicable, get supercharged pressure
        %If supercharge is applicable, get supercharged pressure
        if pSC(1) > max([pBlow, pVapOx])
            %Ox Tank
            try %Empty try-catch to account for funny terminal dynamics
            pVapOx = CoolProp('P', 'Q', 1, 'T', tanksArray(1).T, tanksArray(1).prop); %[Pa]
            rhoLiqOx = CoolProp('D', 'Q', 0, 'T', tanksArray(1).T, tanksArray(1).prop); %[kg/m^3]
            end %If CoolProp throws an error, just leave these properties be
            
            VoidOx = VoidOx + deltam_ox/rhoLiqOx; %[m^3]
            
            %Sum partial pressures for ullage pressure
            pGasOx = nSCOx*(8.31)*tanksArray(1).T/VoidOx;
            tanksArray(1).p = pVapOx + pGasOx; %[Pa]
        end
        
        if pSC(2) > max([pBlow, pVapOx])
            %Fuel Tank
            try %Empty try-catch to account for funny terminal dynamics
            pVapFuel = CoolProp('P', 'Q', 1, 'T', tanksArray(2).T, tanksArray(2).prop); %[Pa]
            rhoLiqFuel = CoolProp('D', 'Q', 0, 'T', tanksArray(2).T, tanksArray(2).prop); %[kg/m^3]
            end %If CoolProp throws an error, just leave these properties be
            
            VoidFuel = VoidFuel + deltam_fuel/rhoLiqFuel; %[m^3]
            
            %Sum partial pressures for ullage pressure
            pGasFuel = nSCFuel*(8.31)*tanksArray(2).T/VoidFuel;
            tanksArray(2).p = pVapFuel + pGasFuel; %[Pa]
        end
        
        if pSC(3) > max([pBlow, pVapBLC])
            %Fuel Tank
            try %Empty try-catch to account for funny terminal dynamics
            pVapBLC = CoolProp('P', 'Q', 1, 'T', tanksArray(3).T, tanksArray(3).prop); %[Pa]
            rhoLiqBLC = CoolProp('D', 'Q', 0, 'T', tanksArray(3).T, tanksArray(3).prop); %[kg/m^3]
            end %If CoolProp throws an error, just leave these properties be
            
            VoidBLC = VoidBLC + deltam_BLC/rhoLiqBLC; %[m^3]
            
            %Sum partial pressures for ullage pressure
            pGasBLC = nSCBLC*(8.31)*tanksArray(3).T/VoidBLC;
            tanksArray(3).p = pVapBLC + pGasBLC; %[Pa]
        end
        
        %If blowdown gas is applied, set pressure to blowdown if tank
        %pressure drops below blowdown pressure
        if pBlow > tanksArray(1).p
            tanksArray(1).p = pBlow;
        end
        if pBlow > tanksArray(2).p
            tanksArray(2).p = pBlow;
        end
        if pBlow > tanksArray(3).p
            tanksArray(3).p = pBlow;
        end
    end
    
    %First guess for new pressure is last chamber pressure
    p1guess = rocket.pc;
    
    %% Calculate & check vapor quality criteria
    if length(tanksArray) == 2
        %Include try-catch to account for wonky draining terminal dynamics
        try
            QOx   = CoolProp("Q", "H", tanksArray(1).h, "D", tanksArray(1).rho, tanksArray(1).prop);
        catch
            QOx = 1;
        end
        
        try
            QFuel = CoolProp("Q", "H", tanksArray(2).h, "D", tanksArray(2).rho, tanksArray(2).prop);
        catch
            QFuel = 1;
        end
        
        try
            QBLC  = CoolProp("Q", "H", tanksArray(2).h, "D", tanksArray(2).rho, tanksArray(2).prop);
        catch
            QBLC  = 1;
        end
        
        %Qmax, or essentially the highest drain %, is what will guide
        %while-loop condition
        Qmax = max([QOx, QFuel]);
    else
        %Include try-catch to account for wonky draining terminal dynamics
        try
            QOx = CoolProp("Q", "H", tanksArray(1).h, "D", tanksArray(1).rho, tanksArray(1).prop);
        catch
            QOx = 1;
        end
        
        try
            QFuel = CoolProp("Q", "H", tanksArray(2).h, "D", tanksArray(2).rho, tanksArray(2).prop);
        catch
            QOx = 1;
        end
        
        try
            QBLC  = CoolProp("Q", "H", tanksArray(3).h, "D", tanksArray(3).rho, tanksArray(3).prop);
        catch
            QBLC = 1;
        end
        %Qmax, or essentially the highest drain %, is what will guide
        %while-loop condition
        Qmax = max([QOx, QFuel, QBLC]); %Don't wanna run out of coolant, either
        
    end
    
    %% Data Output
    
    oxQuality(i) = QOx;
    fuelQuality(i) = QFuel;
    BLCQuality(i) = QBLC;
    
    %Store data in arrays as a function of time
    if length(tanksArray) == 2
        %Ox Tank Data
        oxpvec(i) = tanksArray(1).p;
        oxTvec(i) = tanksArray(1).T;
        oxmvec(i) = tanksArray(1).m;
        %Fuel Tank Data
        fuelpvec(i) = tanksArray(2).p;
        fuelTvec(i) = tanksArray(2).T;
        fuelmvec(i) = tanksArray(2).m;
        %BLC Tank Data
        BLCpvec(i) = tanksArray(2).p;
        BLCTvec(i) = tanksArray(2).T;
        BLCmvec(i) = tanksArray(2).m;
    else
        %Ox Tank Data
        oxpvec(i) = tanksArray(1).p;
        oxTvec(i) = tanksArray(1).T;
        oxmvec(i) = tanksArray(1).m;
        %Fuel Tank Data
        fuelpvec(i) = tanksArray(2).p;
        fuelTvec(i) = tanksArray(2).T;
        fuelmvec(i) = tanksArray(2).m;
        %BLC Tank Data
        BLCpvec(i) = tanksArray(3).p;
        BLCTvec(i) = tanksArray(3).T;
        BLCmvec(i) = tanksArray(3).m;
    end
    
    %Data-dump transient data into objects
    
    pcvec(i) = rocket.pc;
    Tcvec(i) = rocket.Tc;
    gammavec(i) = rocket.gamma;
    Rvec(i) = rocket.R;

    mdotvec(i) = rocket.mdot;
    mdot_oxvec(i) = rocket.mdot_ox;
    mdot_fuelvec(i) = rocket.mdot_fuel;
    mdot_BLCvec(i) = rocket.mdot_BLC;

    Lmat(i,:) = rocket.L_array;
    pmat(i,:) = rocket.p_array;
    Trmat(i,:) = rocket.Tr_array;
    Twmat(i,:) = rocket.Tw_array;
    machmat(i,:) = rocket.Mach_array;
    hmat(i,:) = rocket.h_array;

    thrustvec(i) = rocket.thrust;
    cstarvec(i) = rocket.cstar;
    Cfvec(i) = rocket.Cf;
    pevec(i) = rocket.pe;

    Tmaxvec(i) = rocket.Tmax;
    Tthroatvec(i) = rocket.Tthroat;
    
    clc %Clear time-step print
end

rocketData = transientData;

rocketData.p0Ox = oxpvec;
rocketData.T0Ox = oxTvec;
rocketData.QOx  = oxQuality;
rocketData.mOx  = oxmvec;

rocketData.p0Fuel = fuelpvec;
rocketData.T0Fuel = fuelTvec;
rocketData.QFuel  = fuelQuality;
rocketData.mFuel  = fuelmvec;

rocketData.p0BLC = BLCpvec;
rocketData.T0BLC = BLCTvec;
rocketData.QBLC  = BLCQuality;
rocketData.mBLC  = BLCmvec;

rocketData.pc = pcvec;
rocketData.Tc = Tcvec;
rocketData.gamma = gammavec;
rocketData.R = Rvec;

rocketData.mdot = mdotvec;
rocketData.mdot_ox = mdot_oxvec;
rocketData.mdot_fuel = mdot_fuelvec;
rocketData.mdot_BLC = mdot_BLCvec;

rocketData.L_array = Lmat;
rocketData.p_array = pmat;
rocketData.Tr_array = Trmat;
rocketData.Tw_array = Twmat;
rocketData.Mach_array = machmat;
rocketData.h_array = hmat;

rocketData.thrust = thrustvec;
rocketData.cstar = cstarvec;
rocketData.Cf = Cfvec;
rocketData.pe = pevec;

rocketData.Tmax = Tmaxvec;
rocketData.Tthroat = Tthroatvec;
end