classdef TheTanks < handle
    % Important properties of a propellant tank
    % tankProperties(oxidizer, massOx, ullageOx, temperature, tankAreaInternal);;
    % Public, tunable properties
    properties
        prop; %Propellant name
        m; %Tank mass
        mVap;
        mLiq;
        P; %Tank pressure: Pa
        PVap; %Vapor pressure: Pa
        T; %Tank temperature: K
        h; %Tank specific enthalpy:J/kg 
        hLiq; %Liquid phase specific enthalpy
        rho; %Tank average density
        rhoLiq;
        rhoVap;
        V; %Tank volume:
        height; %Tank height
        ullage; %Tank ullage
        x; %Tank quality
        sLiq; %Liquid phase entropy
        qualUnderThreshold;
        diameter;
        length;
        thickness;
        area;
    end

    methods
        function obj = TheTanks(prop, ullage, T, area, i, OD_SI, height, volume)
            preallocation = 10000;
            
            %Pre-allocation of memory
            obj.T = zeros(1,preallocation);
            obj.ullage = zeros(1,preallocation);
            obj.m = zeros(1,preallocation);
            obj.P = zeros(1,preallocation);
            obj.PVap = zeros(1,preallocation);
            obj.rhoLiq = zeros(1,preallocation);
            obj.rhoVap = zeros(1,preallocation);
            obj.mLiq = zeros(1,preallocation);
            obj.mVap = zeros(1,preallocation);
            obj.x = zeros(1, preallocation);
            obj.h = zeros(1,preallocation);
            obj.hLiq = zeros(1,preallocation);
            obj.rho = zeros(1,preallocation);
            obj.sLiq = zeros(1,preallocation);
            
            
            
            % Perform one-time calculations, such as computing constants
            obj.T(i) = T;
            obj.prop = prop;
            obj.ullage(i) = ullage;
            obj.area = area;
            obj.height = height;
            obj.V = volume;
            obj.diameter = OD_SI;
            %obj.m(i) = m;
            
            obj.P(i) = py.CoolProp.CoolProp.PropsSI("P", "T", T, "Q", 0, prop); %Quality doesn't matter for pressure... come back if matters later for supercharge
            obj.PVap(i) = py.CoolProp.CoolProp.PropsSI("P","T",T,"Q",1,prop);
            
            rhoLiq = py.CoolProp.CoolProp.PropsSI("Dmass", "T", T, "Q", 0, prop); %Initial density of liquid phase ox
            obj.rhoLiq(i) = rhoLiq;
            
            rhoVap = py.CoolProp.CoolProp.PropsSI("Dmass", "T", T, "Q", 1, prop); %Initial density of vapor phase ox
            obj.rhoVap(i) = rhoVap;
            
            m = volume * ((rhoLiq*(1-ullage)) + (rhoVap * ullage));
            %volume = m / ((rhoLiq*(1-ullage)) + (rhoVap * ullage)); %req volume of ox tank
            %obj.V(i) = volume;
            obj.m(i) = m;
            
            massLiq = rhoLiq * volume * (1 - ullage);
            obj.mLiq(i) = massLiq;
            
            massVap = rhoVap * volume * ullage;
            obj.mVap(i) = massVap;
            
            quality = massVap / m;
            obj.x(i) = quality;
            
            obj.h(i) = py.CoolProp.CoolProp.PropsSI("H", "T", T, "Q", quality, prop);
            obj.hLiq(i) = py.CoolProp.CoolProp.PropsSI("H", "T", T, "Q", 0, prop);
            
            obj.rho(i) = m / volume; %average density of the tank

            %height = volume / area; %req usable height of ox tank given dimensions
            
            
            obj.sLiq(i) = py.CoolProp.CoolProp.PropsSI("S", "T", T, "Q", 0, prop);
            
            obj.qualUnderThreshold = true;
            
            
            
        end
        
        
    end
end
