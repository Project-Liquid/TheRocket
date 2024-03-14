classdef tankProperties < handle
    %TANK Tank class with relevant properties
    %Constructor format:
    %
    %   object = tank(p, T, V, prop)
    %
    %Constructor inputs:
    %   
    %   p:      Tank ullage pressure                    [Pa]
    %           Input 'sat' for vapor pressure at T     [String]
    %
    %   T:      Tank average fluid temperature          [K]
    %
    %   h:      Tank specific enthalpy                  [J/kg]
    %
    %   rho:    Tank propellant density                 [kg/m^3]
    %
    %   V:      Tank volume                             [m^3]
    %
    %   prop:   Propellant composition in tank          [String]
    %
    %Methods:
    %
    %   obj.setTemp(T)
    %       Sets propellant tank temperature to T       [K]
    %
    %   obj.setPress(p)
    %       Sets propellant ullage pressure to p        [Pa]
    %
    %   obj.setVolume(V)
    %       Sets propellant tank volume to V            [m^3]
    %
    %   obj.setProp(prop)
    %       Sets tank propellant composition to prop    [String]
    %
    %   obj.saturate()
    %       Sets tank to saturation pressure at already
    %       defined tank temperature
    
    %% PROPERTY INITIALIZATION
    properties
        %Tank pressure property
        p
        %Tank temperature property
        T
        %Tank enthalpy property
        h
        %Tank average density property
        rho
        %Tank volume property
        V
        %Tank mass property
        m
        %Propellant composition property
        prop
        %Initial supercharge volume property (manually defined only)
        VSC
    end
    
    methods
        function obj = tankProperties(p, T, V, prop)
            %TANK Construct an instance of this class
            
            %Collect inputs and add to objects - these are taken raw
            obj.T    = T;
            obj.prop = prop;
            obj.V    = V;
 
            
            %Check for string 'sat' input or low pressure
            pStr = lower(p);
            if strcmp(pStr, 'sat')
                obj.p = CoolProp('P', 'T', obj.T, 'Q', 1, obj.prop);
            elseif p < CoolProp('P', 'T', obj.T, 'Q', 1, obj.prop)
                obj.p = CoolProp('P', 'T', obj.T, 'Q', 1, obj.prop);
            else
                %Otherwise just set pressure to input value
                obj.p = p;
            end
            
            obj.h    = CoolProp("H", "Q", 0, "T", T, prop);
            obj.rho  = CoolProp("D", "Q", 0, "T", T, prop);
            obj.m = obj.rho*obj.V;
        end
        
        %% MANUAL PROPERTY DEFINITION
        function obj = setTemp(obj, T)
            %SETTEMP Sets temperature of prop tank to input T
            obj.T = T;
        end
        
        function obj = setVolume(obj, V)
            %SETVOLUME Sets volume of prop tank to input V
            obj.V = V;
        end
        
        function obj = setPress(obj, p)
            %SETPRESS Sets pressure of prop tank to input p
            obj.p = p;
        end
        
        function obj = setProp(obj, prop)
            %SETPROP Sets composition of propellant to string name 
            %designated by inputprop. If you want to do that. For some
            %reason.
            obj.prop = prop;
        end
        
        
        %% SATURATION OPERATIONS
        function obj = saturate(obj)
            %SATURATE Sets tank to saturation pressure as in self
            %pressurizing blowdown cycles
            obj.p = CoolProp('P', 'T', obj.T, 'Q', 1, obj.prop);
        end
    end
end

